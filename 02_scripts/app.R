# app.R - Analizador de secuencias: FASTA/FASTQ -> MSA & Árbol (versión "pro")
# NOTA IMPORTANTE: NO incluimos install.packages ni BiocManager::install aquí.

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
library(Biostrings)   # manejo de XStringSet
library(msa)          # función msa()
library(ggmsa)        # visualización con ggplot2
library(ape)          # árboles (nj)
library(seqinr)       # write.fasta, dist.alignment helpers
library(ggplot2)
library(ShortRead)

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
# 
# Aumentas el progreso al 10%.
# Link: https://rstudio-pubs-static.s3.amazonaws.com/28353_bf4353b1c63f40f08082d4f91009edef.html
#
            incProgress(0.1)
# 
# El objeto es un reactive() que contiene las secuencias FASTA.
#      
# Link: https://rdrr.io/bioc/Biostrings/man/XStringSet-io.html
#
      s <- seqs()
#
# Alineamiento de secuencias múltiple mediatne el método ClustalW,
# ClustalOmega y MUSCLE, con el paquete msa.
# El if ayuda para elegir el tipo de alineamiento a ecoger y
# tryCatch evita errores en Shiny.
#
# Link: https://stackoverflow.com/questions/30038676/r-trycatch-in-place-with-err-and-warn-handlers-but-shiny-still-crashes   
#
      # tryCatch para manejar errores si faltan ejecutables externos
      alignment <- tryCatch({
        if(method == "ClustalW"){
          msa(s, method = "ClustalW")        
        } else if(method == "ClustalOmega"){
          msa(s, method = "ClustalOmega")    
        } else {
          msa(s, method = "Muscle")          
        }
      }, error = function(e){
#
# El showNotification() muestra NULL, si hay un error.
# Link: https://shiny.posit.co/r/reference/shiny/0.14/shownotification.html
#  
        showNotification(paste("Error en msa():", e$message), type = "error")
        return(NULL)
      })
# 
# Aumentas el progreso al 10%.
# Link: https://rstudio-pubs-static.s3.amazonaws.com/28353_bf4353b1c63f40f08082d4f91009edef.html
#
      incProgress(0.6)
#
# El if ayuda para retornar, si no hay alineamiento. Ayuda de ChatGTP,
# porque sali problemas de error.
#
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

#############################################################################
#############             Interfaz de usuario básica            #############
#############################################################################

library(Biostrings)   # Para manipulación de secuencias biológicas
library(msa)          # Para alineamiento múltiple (MUSCLE)
library(ape)          # Para análisis filogenético y construcción de árboles
library(ggmsa)        # Para visualización del alineamiento

# Leer múltiples archivos FASTA de secuencias ITS y combinarlos
# (Ajustar el patrón o lista de archivos según su directorio de trabajo)
archivos <- list.files(pattern = "^sequence.*\\.fasta$")
secuencias_lista <- lapply(archivos, readDNAStringSet)
todas_secuencias <- do.call(c, secuencias_lista)
# Asignar nombres basados en el nombre de archivo (sin extensión)
names(todas_secuencias) <- sub("\\.fasta$", "", archivos)

# Alinear las secuencias combinadas usando MUSCLE
alineamiento <- msa(todas_secuencias, method = "Muscle")
# Convertir el alineamiento a formato compatible con ape (seqinr)
aline_seqinr <- msaConvert(alineamiento, type = "seqinr::alignment")
aline_dnabin <- as.DNAbin(aline_seqinr)


# Construir árbol filogenético (vecino más cercano - NJ) usando distancia genética (modelo K80)
distancias <- dist.dna(aline_dnabin, model = "K80")
arbol_filogenetico <- nj(distancias)
# Opcional: graficar el árbol filogenético
plot(arbol_filogenetico, main = "Árbol filogenético de Aspergillus (ITS)")

# Función para agregar una nueva secuencia y actualizar el análisis
analizar_nueva_secuencia <- function(nueva_secuencia) {
  # Convertir la nueva secuencia a DNAStringSet y asignarle un nombre
  nueva_set <- DNAStringSet(nueva_secuencia)
  names(nueva_set) <- "Secuencia_Nueva"
  
  # Agregar la nueva secuencia al conjunto existente
  secuencias_actualizadas <- c(todas_secuencias, nueva_set)
  
  # Realizar nuevo alineamiento con todas las secuencias
  aline2 <- msa(secuencias_actualizadas, method = "Muscle")
  # Convertir a formato seqinr y DNAbin para ape
  aline2_seqinr <- msaConvert(aline2, type = "seqinr::alignment")
  aline2_dnabin <- as.DNAbin(aline2_seqinr)
  
  # Construir nuevo árbol filogenético con la secuencia nueva
  dist2 <- dist.dna(aline2_dnabin, model = "K80")
  arbol2 <- nj(dist2)
  plot(arbol2, main = "Árbol filogenético actualizado")
  
  # Evaluar relación filogenética: encontrar la especie más cercana en el árbol
  dmat <- cophenetic.phylo(arbol2)
  nombre_nueva <- names(secuencias_actualizadas)[length(secuencias_actualizadas)]
  otros_nombres <- names(secuencias_actualizadas)[-length(secuencias_actualizadas)]
  cercano <- otros_nombres[which.min(dmat[nombre_nueva, otros_nombres])]
  cat("La secuencia nueva es filogenéticamente más cercana a:", cercano, "\n")
  
  # Calcular contenido GC y conteo de cada nucleótido de la nueva secuencia
  frec <- alphabetFrequency(DNAString(nueva_secuencia), baseOnly = TRUE)
  gc <- frec["G"] + frec["C"]
  total <- sum(frec[c("A", "C", "G", "T")])
  contenido_gc <- round( (gc/total) * 100, 2 )
  cat("Contenido GC (%):", contenido_gc, "\n")
  cat("Conteo nucleótidos (A, C, G, T):", frec["A"], frec["C"], frec["G"], frec["T"], "\n")
  
  # Visualizar -las primeras 100 posiciones alineadas usando ggmsa
  # Convertir el alineamiento actualizado a un objeto DNAStringSet para exportar
  aline2_set <- as(aline2, "DNAStringSet")
  writeXStringSet(aline2_set, file = "alineamiento_actualizado.fasta")
  print(ggmsa("alineamiento_actualizado.fasta", start = 1, end = 100, color = "Chemistry_NT"))
}

# Ejemplo de uso de la función:
# nueva_seq <- "ATGCGTAACGTAGCTAGCTAGCTAGCATCGATCG..."
# analizar_nueva_secuencia(nueva_seq)
Nueva1 <- "AACGACCCCCCAGAGCCGGAAAGTTGGTCAAACCCGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTACCGAGTGCGGGTCTTTATGGCCCAACCTCCCACCCGTGACTATTGTACCTTGTTGCTTCGGCGGGCCCGCCAGCGTTGCTGGCCGCCGGGGGGCGACTCGCCCCCGGGCCCGTGCCCGCCGGAGACCCCAACATGAACCCTGTTCTGAAAGCTTGCAGTCTGAGTTGTGATTCTTTGCAATCAGTTAAAACTTTCAACAATGGATCTCTTGGTTCCGGCATCGATGAAGAACGCAGCGAAATGCGATAACTAATGTGAATTGCAGAATTCAGTGAATCATCGAGTCTTTGAACGCACATTGCGCCCCCTGGTATTCCGGGGGGCATGCCTGTCCGAGCGTCATTGCTGCCCTCAAGCCCGGCTTGTGTGTTGGGCCCTCGTCCCCCGGCTCCCGGGGGACGGGCCCGAAAGGCAGCGGCGGCACCGCGTCCGGTCCTCGAGCGTATGGGGCTTCGTCTTCCGCTCCGTAGGCCCGGCCGGCGCCCGCCGACGCATT"


analizar_nueva_secuencia(Nueva1)






