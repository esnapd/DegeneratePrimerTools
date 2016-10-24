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
        design_primers( maxdegeneracies=c(1,10), number_iterations=12)
     
     #obtain the locations of peaks 
     primerlocs <- autofind_primers(dp, keepprimers = 4)
     
     # Region1
     msa1 <- add_primers_to_MSA(dp, position = primerlocs[[1]])
     
     msa2 <- add_primers_to_MSA(dp, position = primerlocs[[2]])
     msa3 <- add_primers_to_MSA(dp, position = primerlocs[[3]])
     msa4 <- add_primers_to_MSA(dp, position = primerlocs[[4]])
     
     # Region 4
     tabsetPanel(
       id = "mpanel", 
       type = "pill",
       tabPanel("Full Alignment", renderMsaR({msaR(dp@msa)})),
       tabPanel("Region 1", renderMsaR({msaR(msa1)})),
       tabPanel("Region 2", renderMsaR({msaR(msa2)})),
       tabPanel("Region 3", renderMsaR({msaR(msa3)})),
       tabPanel("Region 4", renderMsaR({msaR(msa4)})))
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
