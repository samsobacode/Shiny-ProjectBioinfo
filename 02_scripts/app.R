# app.R - Analizador de secuencias: FASTA/FASTQ -> MSA & Árbol (versión "pro")
# NOTA IMPORTANTE: NO incluimos install.packages ni BiocManager::install aquí.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22")

BiocManager::install(c("GenomeInfoDb", "Rhtslib"), force = TRUE)
BiocManager::install("Biostrings", force = TRUE)
BiocManager::install("msa", force = TRUE)
BiocManager::install("ggmsa", force = TRUE)

install.packages("shiny")
install.packages("ape")
install.packages("ggplot2")
install.packages("seqinr")

# Instala dependencias fuera de la app antes de ejecutarla (ver sección "Instalación" abajo).

library(shiny)
library(Biostrings)   # manejo de XStringSet
library(msa)          # función msa()
library(ggmsa)        # visualización con ggplot2
library(ape)          # árboles (nj)
library(seqinr)       # write.fasta, dist.alignment helpers
library(ggplot2)

#############################################################################
#############             Interfaz de usuario básica            #############
#############################################################################


ui <- fluidPage(
  titlePanel("Análisis bioinformático"),
  sidebarLayout(
    sidebarPanel(
      fileInput("archivo", "Sube archivo FASTA (.fasta, .fa, .fna)",
                accept = c(".fa", ".fasta", ".fna")),
      hr(),
      checkboxInput("showConsensus", "Mostrar consenso al imprimir MSA", value = TRUE)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Resumen",
                 verbatimTextOutput("tipoDetected"),
                 tableOutput("seqSummary")
        ),
        
        # -------------------- MSA ---------------------
        tabPanel("MSA (texto)",
                 
                 # Controles movidos desde el sidebar:
                 selectInput("method", "Algoritmo MSA:",
                             choices = c("ClustalW", "ClustalOmega", "Muscle"),
                             selected = "ClustalW"),
                 
                 actionButton("runMSA", "Ejecutar MSA"),
                 
                 downloadButton("downloadAlignedFasta", "Descargar MSA (FASTA)"),
                 
                 hr(),
                 verbatimTextOutput("msaText")
        ),
        
        # ------------------------------------------------
        tabPanel("Matriz de distancias",
                 tableOutput("distMatrix")
        ),
        tabPanel("Árbol filogenético (NJ)",
                 plotOutput("treePlot", height = "600px")
        )
      )
    )
  )
)

#############################################################################
#############             Interfaz de usuario básica            #############
#############################################################################

server <- function(input, output, session) {
  
  # Reactive: leer secuencias con Biostrings (intenta DNA, si detecta U -> RNA)
  seqs <- reactive({
    req(input$archivo)
    path <- input$archivo$datapath
    
    # Intentar leer como DNA; si detecta 'U' y no 'T' convertimos a RNA
    tryDNA <- try(readDNAStringSet(path), silent = TRUE)
    if(!inherits(tryDNA, "try-error")) {
      merged <- paste(as.character(unlist(tryDNA)), collapse = "")
      if(grepl("U", merged) & !grepl("T", merged)) {
        # Es probable que sea RNA
        rna <- readRNAStringSet(path)
        return(rna)
      } else {
        return(tryDNA)
      }
    } else {
      # Intentar RNA si DNA fail
      tryRNA <- try(readRNAStringSet(path), silent = TRUE)
      if(!inherits(tryRNA, "try-error")) return(tryRNA)
      stop("No se pudo leer el archivo como DNA ni RNA. Revisa el formato FASTA.")
    }
  })
  
  output$tipoDetected <- renderText({
    req(seqs())
    tipo <- ifelse(is(seqs(), "DNAStringSet"), "DNA", "RNA")
    paste0("Tipo detectado (según Biostrings): ", tipo, " — secuencias: ", length(seqs()))
  })
  
  output$seqSummary <- renderTable({
    req(seqs())
    df <- data.frame(
      ID = names(seqs()),
      Length = width(seqs()),
      stringsAsFactors = FALSE
    )
    df
  })

--------------------------------------------------------------------------------
#                                          MSA
--------------------------------------------------------------------------------
  
  # ReactiveValues para guardar alignment y objetos derivados
  rv <- reactiveValues(alignment = NULL, aln_seqinr = NULL, distmat = NULL, tree = NULL)
  
  observeEvent(input$runMSA, {
    req(seqs())
    method <- input$method
    
# Crea el cuadro para ejecutar el análisis con msa.
# Función de Shiny withProgress:
# La clase de referencia proporciona una API orientada a objetos.
#
# message: Un vector de caracteres de un solo elemento; el mensaje que se 
# mostrará al usuario, o NULL para ocultar el mensaje actual (si lo hay).
#
# Link: https://shiny.posit.co/r/reference/shiny/1.0.1/withprogress.html
#
    withProgress(message = paste("Ejecutando MSA con", method, "..."), value = 0, {
# Aumentas el progreso al 10%.
            incProgress(0.1)
# 
# El objeto es un reactive() que lee las secuencias FASTA.
#
      s <- seqs()
      
      # tryCatch para manejar errores si faltan ejecutables externos
      alignment <- tryCatch({
        if(method == "ClustalW"){
          msa(s, method = "ClustalW")        # usa ClustalW (si está disponible)
        } else if(method == "ClustalOmega"){
          msa(s, method = "ClustalOmega")    # Clustal Omega
        } else {
          msa(s, method = "Muscle")          # MUSCLE
        }
      }, error = function(e){
        showNotification(paste("Error en msa():", e$message), type = "error")
        return(NULL)
      })
      
      incProgress(0.6)
      
      if(is.null(alignment)) return()
      
      rv$alignment <- alignment
      
      # Convertir a objeto seqinr::alignment para distancia y escritura fasta
      aln_seqinr <- tryCatch({
        msaConvert(alignment, type = "seqinr::alignment")
      }, error = function(e){
        showNotification(paste("Error en msaConvert():", e$message), type = "error")
        return(NULL)
      })
      
      rv$aln_seqinr <- aln_seqinr
      
      incProgress(0.85)
      
      # Calcular matriz de distancias (identity -> proporción de diferencias)
      if(!is.null(aln_seqinr)){
        d <- tryCatch({
          dist.alignment(aln_seqinr, "identity")
        }, error = function(e){
          showNotification(paste("Error en dist.alignment():", e$message), type = "error")
          return(NULL)
        })
        rv$distmat <- as.matrix(d)
        # Construir árbol NJ
        if(!is.null(d)){
          rv$tree <- tryCatch({
            nj(d)
          }, error = function(e){
            showNotification(paste("Error en nj():", e$message), type = "error")
            return(NULL)
          })
        }
      }
      
      incProgress(1)
    }) # withProgress
  }) # observeEvent runMSA
  
  # Mostrar MSA como texto
  output$msaText <- renderPrint({
    req(rv$alignment)
    if(input$showConsensus){
      print(rv$alignment, show = "alignment")
    } else {
      print(rv$alignment, show = "alignment", showConsensus = FALSE)
    }
  })
  
  # Mostrar matriz de distancias
  output$distMatrix <- renderTable({
    req(rv$distmat)
    # formateo legible
    round(rv$distmat, 4)
  }, rownames = TRUE)
  
  # Dibuja el árbol NJ
  output$treePlot <- renderPlot({
    req(rv$tree)
    plot(rv$tree, main = paste("Árbol NJ -", input$method))
    # opcional: añadir tip labels rotadas si muchas secuencias
    tiplabels()
  })
  
  # Download handler: descarga el MSA en FASTA (secuencias alineadas)
  output$downloadAlignedFasta <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$fastaFile$name), "_aligned.fasta")
    },
    content = function(file) {
      req(rv$aln_seqinr)
      aln <- rv$aln_seqinr
      # seqinr::write.fasta requiere lista de secuencias (character) y nombres
      seqs_list <- as.list(aln$seq)
      names(seqs_list) <- aln$nam
      write.fasta(sequences = seqs_list, names = names(seqs_list), file.out = file)
    }
  )

#############################################################################
#############             Interfaz de usuario básica            #############
#############################################################################
  
} # server

shinyApp(ui, server)