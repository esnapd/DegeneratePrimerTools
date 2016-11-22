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
     primerlocs <- autofind_primers(dp, keepprimers = 4)
     
     # Add all Peak Info to the MSA
     msa1 <- add_primers_to_MSA(dp, positions = primerlocs)
     t1 <- data.frame(dp@primerdata) %>% 
       filter(Pos %in% primerlocs) %>% 
       select(Pos, PrimerSeq, PrimerDeg, degeneracy, coverage) %>%
       arrange(Pos)

     # Create the Return MSA
     rowheight <- 15
     mainPanel(
       id = "mpanel", 
       renderMsaR({ msaR(msa1, menu = F, alignmentHeight = nrow(msa1)*rowheight, 
                        leftheader = FALSE, labelNameLength = 160, seqlogo=F)}),
       
       DT::renderDataTable(t1, rownames = FALSE, options = list(dom = 't'))
     )

       #          renderTable(data.frame(code = names(Biostrings::IUPAC_CODE_MAP),
       #                                 nucleotides = Biostrings::IUPAC_CODE_MAP))))
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
