library(DegeneratePrimerTools)
library(Biostrings)
library(msaR)

library(dplyr)
library(shiny)

#' Server Function
#'
server <- function(input, output) {
  
  # capture the input in a reactive component
  inFile <- reactive(
    if (is.null(input$file1)){
      NULL
    }else {
      readDNAStringSet(input$file1$datapath)
    }
  )
  
  
  # return the display values
  output$mainpanel <- renderUI({
    if(is.null(inFile())) { 
      return ()
    } else {
     dp <- degeprimer(inFile()) %>%
        run_alignment() %>%
        build_tree() %>%
        design_primers( maxdegeneracies=c(1, 10, 20), number_iterations=10)
      
      tabsetPanel(id = "mpanel", type = "pill"
                  ,tabPanel("Full Alignment", renderMsaR({msaR(dp@msa)}))
                  ,tabPanel("Region 1", renderMsaR({msaR(dp@msa)}))
                  ,tabPanel("Region 2", renderMsaR({msaR(dp@msa)}))
                  ,tabPanel("Region 3", renderMsaR({msaR(dp@msa)})))
    }
  })
}

#' UI Function
#'
#' Handles file input and proveides a generic "Main Panel" output
#' The content of that output is controlled by the server function
ui <- fluidPage(
    titlePanel("Uploading A Fasta File"),
    sidebarLayout(
      sidebarPanel(
        fileInput('file1', 'Choose a Fasta File',
                  accept=c('.fasta',
                           '.fa',
                           '.fna'))
      ),
      mainPanel(uiOutput("mainpanel"))
    )
  )


#' Call the Widget
#'
shinyApp(ui = ui, server = server)
