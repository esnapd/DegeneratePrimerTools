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
      dnas <- readDNAStringSet(input$file1$datapath)
      if (length(dnas) < 2) {
        stop("This fastafile does not have enough sequences to align. 
             If you believe this message is in error check your fasta file for 
             the use of carriage returns used to demarcate line ends.")
      }
      dnas
    }
  )
  
  
  # this determines file upload status and is used for our UI to determine what to display
  output$fileUploaded <- reactive({return(is.null(inFile()))})
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  
  output$fastaName <- reactive(input$file1$name)
  
  # create the IUPAC table
  output$IUPAC <- renderTable(data.frame(
    code = names(Biostrings::IUPAC_CODE_MAP),
    nucleotides = Biostrings::IUPAC_CODE_MAP))
  
  # obtain data from file and return the main panel output
  output$mainpanel <- renderUI({
    if(is.null(inFile())) { 
      return ()
    } else {
      
     sliderLength <- as.numeric(input$sliderLength)
      
     dp <- degeprimer(inFile()) %>%
        run_alignment() %>%
        build_tree() %>%
        design_primers(oligolengths = sliderLength, maxdegeneracies=as.numeric(input$checkGroup), number_iterations=10, ncpus = 4)
     
     # obtain the locations of peaks and check for the presence of NAs
      primerlocs   <- autofind_primers(dp, keepprimers = input$numberofsites, minsequences = input$minseqs)
      primerlocNAs <- any(is.na(primerlocs))
     
     #primerlocNAs <- 1
     # if autofind has difficulties, provide one output,
     # otherwise provide a main output.
     rowheight <- 15
     
     if (primerlocNAs) {
       mainPanel(
         id = "mpanel", 
         h2("Warning! There is a problem with the autodetection of primers."),
         
         br(),
         
         h3("MSA of the input sequences."),
         renderMsaR({ msaR(dp@msa, 
                           menu = F, 
                           height= nrow(dp@msa)*rowheight, 
                           alignmentHeight = nrow(dp@msa)*rowheight, 
                           leftheader = FALSE,
                           labelNameLength = 160, 
                           labelid = FALSE,
                           seqlogo=F)}),
         br(),
         h3("Primers Designed against this MSA"),
         
         p("This plot shows the coverage(y-axis) of primers designed at each location (x-axis)
           across your ungapped alignment using a maximum degeneracy value of (color). If you are seeing this
           plot there was an error in choosing primers. The auto-detection of primers looks
           for peaks in this graph. If you see very few peaks, or peaks with very low coverage use the
           MSA and then plot to troubleshoot. Most likely you will need to increase your maximum degeneracy."),
         
         shiny::renderPlot(plot_degeprimer(dp@primerdata)
         )
       )
     } else {
       
       # create the MSA and Table and Pass them to the Main Panel
       msa1 <- add_primers_to_MSA(dp, positions = primerlocs, mode = "consensus")
       
       t1 <- data.frame(dp@primerdata) %>% 
         filter(Pos %in% primerlocs) %>% 
         mutate(RevComp= as.character(reverseComplement(DNAStringSet(PrimerSeq)))) %>%
         select(Pos, PrimerSeq, RevComp, PrimerDeg, degeneracy, coverage) %>%
         arrange(Pos)
       
       mainPanel(
         id = "mpanel", 
           
         renderMsaR({ msaR(msa1, 
                           menu = F, 
                           height= nrow(msa1)*rowheight, 
                           alignmentHeight = nrow(msa1)*rowheight, 
                           leftheader = FALSE, 
                           labelNameLength = 160, 
                           labelid = FALSE,
                           seqlogo=F)}),
         
         DT::renderDataTable(t1, 
                             rownames = FALSE, 
                             options = list(
                               dom = 't',
                               pageLength = 50))
         )
     }
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
        
        # Include this Conditional Panel BEFORE File Upload
        conditionalPanel(
          condition="output.fileUploaded",
          
          p("Choose a file to upload and the degeneracy values you would like
            to pass to DEGEPRIME. This program will design primers against your sequences 
            using each of the specified values.Higher values will generate more degenerate primers at the expense of speed."),
          
          fileInput('file1', 
                    label=h4('Choose a Fasta File'),
                    accept=c('.fasta', '.fa','.fna')
          ),
          
          hr(),
          
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
          
          
          sliderInput("sliderLength", "Length of primers:",
                      min = 10, max = 30,
                      value = 21),
          
          numericInput("numberofsites", "Number of Primer Locations to Return:", 6, min = 1, max = 1000),
          numericInput("minseqs", "Minimum Sequence to Be Considered for Primer Design:", 3, min = 2, max = 10)
        ),
        
        # Include this Conditional Panel AFTER File Upload
        conditionalPanel(
          condition="!output.fileUploaded",
          
          h4("Primer Design Results"),
          
          p("To the right are two panels. The top displays a multiple seqeuence alignment
            of the input sequences along with the primers designed using the degeneracy values you have specified.
            Highlighting the schematic will navigate around the multiple sequence alignment.
            Below is a table of primers."),
          
          br(),
          
          p("Fasta file:"),
          htmlOutput("fastaName"),
          htmlOutput("input$length"),
          
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
