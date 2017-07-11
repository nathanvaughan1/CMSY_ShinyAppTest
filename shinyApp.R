#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
setwd("/home/enrico/Work/BlueBridge/CMSY/NathanVersion/")
source(paste0(getwd(),"/CMSY_Parametrized.R"))



# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("DLMTools CMSY Test page"),
   
   sidebarLayout(
      sidebarPanel(
        fileInput("file1", "Choose Stock CSV File",
                  accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv")
        ),
        fileInput("file2", "Choose Stock ID CSV File",
                  accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv")
        ),
        uiOutput("fill"),
        numericInput("kv_pairs", "Select K V Pairs", 100, min = 100, max = 10000, step=100),
        numericInput("uncert", "Uncertainty in catch", 0.1, min = 0.1, max = 2, step=0.1),
        numericInput("sigmaR", "Process error for CMSY", 0.1, min = 0.1, max = 2, step=0.1),
        numericInput("ni", "Iterations for r-k-startbiomass combinations", 3, min = 0.1, max = 4, step=0.1),
        numericInput("nab", "Minimum number of years with abundance data to run BSM", 5, min = 2, max = 10, step=1),
        actionButton("go", "Go")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        plotlyOutput("vectorized"),
        HTML("<br><br><br>"),
        plotlyOutput("original"),
        textOutput("renderInfo")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   cmsy <- reactiveValues()
   
   output$fill <- renderUI({
     inFile1 <- input$file1
     inFile2 <- input$file2
     
     if (is.null(inFile1)) {
       return(NULL)
     }
     if (is.null(inFile2)) {
       return(NULL)
     }
     
     a <- read.csv(inFile1$datapath)
     

     selectInput("stocks", "Select a stock", sort(unique(a$Stock)))
   })
   
   output$renderInfo <- renderText({
     paste0(cmsy$res$comments)
   })

   output$vectorized <- renderPlotly({ 
     if ("res" %in% names(cmsy)) {
       vectorizedPlot <- plot_ly(as.data.frame(cmsy$res$Vect), x = ~X1, y = ~X2,type='scatter',mode='markers')
       vectorizedPlot <- layout(vectorizedPlot, title="Vectorized plot")
       vectorizedPlot
     } else {
       return (NULL)
     }
     
  })
   
   output$original <- renderPlotly({ 
     if ("res" %in% names(cmsy)) {
       originalPlot <- plot_ly(as.data.frame(cmsy$res$Orig), x = ~X1, y = ~X2,type='scatter',mode='markers')
       originalPlot <- layout(originalPlot, title="Original plot")
       originalPlot
     } else {
       return (NULL)
     }
   })
  
   observeEvent(input$go, {
     inFile1 <- input$file1
     inFile2 <- input$file2
     
     if (is.null(inFile1)) {
       return(NULL)
     }
     if (is.null(inFile2)) {
       return(NULL)
     }
     
     ret = runCMSY(inFile1$datapath, inFile2$datapath, input$stocks, input$uncert, input$sigmaR, input$kv_pairs, input$ni, input$nab)
     cmsy$res <- ret
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

