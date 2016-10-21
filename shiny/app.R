library(DegeneratePrimerTools)
library(Biostrings)
library(msaR)

library(dplyr)
library(shiny)

server <- function(input, output) {
    output$msa <- renderMsaR({
      inFile <- input$file1
      
      if (is.null(inFile)) return(NULL)
      
      dp <- degeprimer(readDNAStringSet(inFile$datapath)) %>%
        run_alignment() %>%
        build_tree() %>%
        design_primers( maxdegeneracies=c(1, 10, 20), number_iterations=10)
      
      
      msaR(dp@msa) 
    })
}

ui <- fluidPage(
    titlePanel("Uploading A Fasta File"),
    sidebarLayout(
      sidebarPanel(
        fileInput('file1', 'Choose a Fasta File',
                  accept=c('.fasta',
                           '.fa',
                           '.fna'))
      ),
      mainPanel(
        msaROutput("msa", width="100%")
      )
    )
  )

shinyApp(ui = ui, server = server)
