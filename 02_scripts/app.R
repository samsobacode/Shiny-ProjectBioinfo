#############################################################################
# app.R - Analizador bioinformático
#############################################################################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22")

BiocManager::install(c("GenomeInfoDb", "Rhtslib"), force = TRUE)
BiocManager::install("Biostrings", force = TRUE)
BiocManager::install("ShortRead", force = TRUE)
BiocManager::install("msa", force = TRUE)
BiocManager::install("ggmsa", force = TRUE)

install.packages("shiny")
install.packages("ape")
install.packages("ggplot2")
install.packages("seqinr")

# Instala dependencias fuera de la app antes de ejecutarla (ver sección "Instalación" abajo).

library(shiny)
library(ggplot2)
library(ShortRead)
library(shiny)
library(Biostrings)
library(msa)
library(ggmsa)
library(ape)
library(seqinr)
library(ggplot2)

#############################################################################
#                                UI
#############################################################################

ui <- fluidPage(
  titlePanel("Análisis bioinformático (FASTA → MSA → Árbol)"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("archivo", "Sube archivo FASTA",
                accept = c(".fa", ".fasta", ".fna")),
      checkboxInput("showConsensus", "Mostrar consenso", TRUE),
      hr(),
      h4("MSA"),
      selectInput("method", "Algoritmo:",
                  choices = c("ClustalW", "ClustalOmega", "Muscle")),
      actionButton("runMSA", "Ejecutar MSA"),
      downloadButton("downloadAlignedFasta", "Descargar MSA (FASTA)")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Resumen",
                 verbatimTextOutput("tipoDetected"),
                 tableOutput("seqSummary")),
        
        tabPanel("MSA (texto)",
                 verbatimTextOutput("msaText")),
        
        tabPanel("Matriz de distancias",
                 tableOutput("distMatrix")),
        
        tabPanel("Árbol filogenético",
                 plotOutput("treePlot", height = "600px"))
      )
    )
  )
)

#############################################################################
#                                SERVER
#############################################################################

server <- function(input, output, session) {
  
  ###########################################################################
  #           Lectura FASTA (DNA o RNA)
  ###########################################################################
  
  seqs <- reactive({
    req(input$archivo)
    path <- input$archivo$datapath
    
    # Intentar DNA
    tryDNA <- try(readDNAStringSet(path), silent = TRUE)
    
    if(!inherits(tryDNA, "try-error")) {
      merged <- paste(as.character(unlist(tryDNA)), collapse = "")
      
      if(grepl("U", merged) & !grepl("T", merged)) {
        return(readRNAStringSet(path))   # Es RNA
      } else {
        return(tryDNA)                   # Es DNA
      }
    }
    
    # Intentar RNA si DNA falló
    tryRNA <- try(readRNAStringSet(path), silent = TRUE)
    if(!inherits(tryRNA, "try-error"))
      return(tryRNA)
    
    stop("No se pudo leer como DNA ni RNA.")
  })
  
  output$tipoDetected <- renderText({
    req(seqs())
    tipo <- ifelse(is(seqs(), "DNAStringSet"), "DNA", "RNA")
    paste("Tipo detectado:", tipo, " | Secuencias:", length(seqs()))
  })
  
  output$seqSummary <- renderTable({
    req(seqs())
    data.frame(
      ID = names(seqs()),
      Length = width(seqs()),
      stringsAsFactors = FALSE
    )
  })
  
  ###########################################################################
  #              MSA, distancias, árbol
  ###########################################################################
  
  rv <- reactiveValues(
    alignment = NULL,
    aln_seqinr = NULL,
    distmat = NULL,
    tree = NULL
  )
  
  ###########################################################################
  #               Ejecutar MSA
  ###########################################################################
  
  observeEvent(input$runMSA, {
    req(seqs())
    method <- input$method
    
    withProgress(message = paste("Ejecutando MSA con", method), value = 0, {
      incProgress(0.1)
      
      s <- seqs()
      
      alignment <- tryCatch({
        switch(method,
               "ClustalW"     = msa(s, method="ClustalW"),
               "ClustalOmega" = msa(s, method="ClustalOmega"),
               "Muscle"       = msa(s, method="Muscle"))
      }, error = function(e){
        showNotification(paste("Error en msa():", e$message), type="error")
        return(NULL)
      })
      
      if(is.null(alignment)) return()
      rv$alignment <- alignment
      
      incProgress(0.5)
      
      # Convertir MSA a formato seqinr
      rv$aln_seqinr <- tryCatch({
        msaConvert(alignment, type = "seqinr::alignment")
      }, error = function(e){
        showNotification("Error en msaConvert()", "error")
        return(NULL)
      })
      
      incProgress(0.75)
      
      # Matriz de distancias
      if(!is.null(rv$aln_seqinr)) {
        d <- tryCatch({
          dist.alignment(rv$aln_seqinr, "identity")
        }, error=function(e){
          showNotification("Error en dist.alignment()", "error")
          return(NULL)
        })
        
        rv$distmat <- as.matrix(d)
        
        # Árbol NJ
        rv$tree <- tryCatch({
          nj(rv$distmat)
        }, error=function(e){
          showNotification("Error en nj()", "error")
          return(NULL)
        })
      }
      
      incProgress(1)
    })
  })
  
  ###########################################################################
  #                       Mostrar MSA en texto
  ###########################################################################
  
  output$msaText <- renderPrint({
    req(rv$alignment)
    
    if(input$showConsensus){
      print(rv$alignment, show="alignment")
    } else {
      print(rv$alignment, show="alignment", showConsensus=FALSE)
    }
  })
  
  ###########################################################################
  #                        Matriz de distancias
  ###########################################################################
  
  output$distMatrix <- renderTable({
    req(rv$distmat)
    round(rv$distmat, 4)
  }, rownames = TRUE)
  
  ###########################################################################
  #                        Árbol filogenético
  ###########################################################################
  
  output$treePlot <- renderPlot({
    req(rv$tree)
    plot(rv$tree, main = paste("Árbol NJ -", input$method))
  })
  
  ###########################################################################
  #                     Descargar MSA en FASTA
  ###########################################################################
  
  output$downloadAlignedFasta <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$archivo$name), "_aligned.fasta")
    },
    content = function(file) {
      req(rv$aln_seqinr)
      aln <- rv$aln_seqinr
      
      seqs_list <- as.list(aln$seq)
      names(seqs_list) <- aln$nam
      
      write.fasta(sequences=seqs_list, names=names(seqs_list), file.out=file)
    }
  )
}

#############################################################################
#                              LAUNCH APP
#############################################################################

shinyApp(ui, server)