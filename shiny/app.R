library(DegeneratePrimerTools)
library(Biostrings)
library(msaR)
library(DT)
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
        design_primers(maxdegeneracies=as.numeric(input$checkGroup), number_iterations=10, ncpus = 4)
     
     #obtain the locations of peaks 
     primerlocs <- autofind_primers(dp, keepprimers = 5)
     
     # Region1
     msa1 <- add_primers_to_MSA(dp, position = primerlocs[[1]])
     
     msa2 <- add_primers_to_MSA(dp, position = primerlocs[[2]], strict = FALSE)
     msa3 <- add_primers_to_MSA(dp, position = primerlocs[[3]], strict = FALSE)
     msa4 <- add_primers_to_MSA(dp, position = primerlocs[[4]], strict = FALSE)
     
     t1 <- data.frame(dp@primerdata) %>% filter(Pos==primerlocs[[1]])
     t2 <- data.frame(dp@primerdata) %>% filter(Pos==primerlocs[[2]])
     t3 <- data.frame(dp@primerdata) %>% filter(Pos==primerlocs[[3]])
     t4 <- data.frame(dp@primerdata) %>% filter(Pos==primerlocs[[4]])
       
     DT::renderDataTable(t1)
     # Region 4
     tabsetPanel(
       id = "mpanel", 
       type = "pill",
       tabPanel("Full Alignment", 
                renderMsaR({msaR(dp@msa, menu = F, alignmentHeight = 300, labelNameLength = 250, seqlogo=F)})),
       tabPanel("Region 1",
                renderMsaR({msaR(msa1, menu = F, alignmentHeight = 300, labelNameLength = 250, seqlogo=F)}),
                DT::renderDataTable(t1)),
       tabPanel("Region 2",
                renderMsaR({msaR(msa2, menu = F, alignmentHeight = 300, labelNameLength = 250, seqlogo=F)}),
                DT::renderDataTable(t2)),
       tabPanel("Region 3",
                renderMsaR({msaR(msa3, menu = F, alignmentHeight = 300, labelNameLength = 250, seqlogo=F)}),
                DT::renderDataTable(t3)),
       tabPanel("Region 4",
                renderMsaR({msaR(msa4, menu = F, alignmentHeight = 300, labelNameLength = 250, seqlogo=F)}),
                DT::renderDataTable(t4)))
    }
  })
}

#' UI Function
#'
#' Handles file input and proveides a generic "Main Panel" output
#' The content of that output is controlled by the server function
ui <- fluidPage(
    titlePanel("Upload A Fasta File to Design Degenerate Primers"),
    sidebarLayout(
      sidebarPanel(
        fileInput('file1', 'Choose a Fasta File',
                  accept=c('.fasta',
                           '.fa',
                           '.fna')),
      checkboxGroupInput(
      "checkGroup", 
       label = h6("Choose Degeneracy Values"), 
       choices = list("Degeneracy 10" = 10, 
                      "Degeneracy 50" = 50, 
                      "Degeneracy 100" = 100,
                      "Degeneracy 200" = 200,
                      "Degeneracy 500" = 500,
                      "Degeneracy 1000" = 1000,
                      "Degeneracy 2000" = 2000),
       selected = c(50, 100, 200, 500))
      
      ),
      mainPanel(uiOutput("mainpanel"))
    )
  )


#' Call the Widget
#'
shinyApp(ui = ui, server = server)
