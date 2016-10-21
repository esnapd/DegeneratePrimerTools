library(shiny)
library(Biostrings)
library(msaR)


server <- function(input, output) {
    output$msa <- renderMsaR({
      inFile <- input$file1
      
      if (is.null(inFile))
        return(NULL)
      
      msaR(inFile$datapath) 
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
