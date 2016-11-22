library(DegeneratePrimerTools)
library(Biostrings)
library(msaR)
library(DT)
library(dplyr)
library(shiny)
library(xtable)

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
  
  # this determines file upload status and is used for our UI to determine what to display
  output$fileUploaded <- reactive({return(is.null(inFile()))})
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  
  output$IUPAC <- renderTable(data.frame(
    code = names(Biostrings::IUPAC_CODE_MAP),
    nucleotides = Biostrings::IUPAC_CODE_MAP))
  
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
       
       DT::renderDataTable(t1, 
                           rownames = FALSE, 
                           options = list(
                             dom = 't',
                             pageLength = 5))
     )

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
        
        # Include this Conditonal Panel BEFORE File Upload
        conditionalPanel(
          condition="output.fileUploaded",
          
          p("Choose a file to upload and the degeneracy values you would like
            to pass to DEGEPRIME. This program will desing primers against your sequences 
            useing each of the specified values."),
          
          p("Higher values will generate more degenerate primers at the expense of speed."),
          
          checkboxGroupInput(
            "checkGroup", 
            label = h4("Choose Degeneracy Values"), 
            choices = list("Degeneracy 10" = 10, 
                           "Degeneracy 50" = 50, 
                           "Degeneracy 100" = 100,
                           "Degeneracy 200" = 200,
                           "Degeneracy 500" = 500,
                           "Degeneracy 1000" = 1000,
                           "Degeneracy 2000" = 2000),
            selected = c(50, 100, 200, 500)),
          
          hr(),
          
          fileInput('file1', 
                    label=h4('Choose a Fasta File'),
                    accept=c('.fasta', '.fa','.fna')
          )
        ),
      
        conditionalPanel(
          condition="!output.fileUploaded",
          
          h4("Primer Design Results"),
          
          p("To the right are two panels. The top displays a multiple seqeuence alignment
            of the input sequences along with the primers designed using the degeneracy values you have specified.
            Highlighting the schematic will navigate around the multiple sequence alignment.
            Below is a table of primers."),
          
          br(),
          
          h4("IUPAC Degenerate Nucleotide Codes"),
          
          tableOutput("IUPAC")
          )
      ),
      
      mainPanel(uiOutput("mainpanel"))
    )
  )


#' Call the Widget
#'
shinyApp(ui = ui, server = server)
