#appR
library(shiny)
library(shinythemes)

# Define UI ----
ui <- fluidPage(theme = shinytheme("yeti"),
                navbarPage(
                  "FCSimple transformer",
                  tabPanel("transform",
                           sidebarPanel(
                             uiOutput("raw_data"),
                             uiOutput("transform_type"),
                             uiOutput("hyperparameters")
                           ),
                           mainPanel(
                             plotOutput("densPlot")
                           )),
                  tabPanel("translate")
                  )
)

# Define server logic ----
server <- function(input, output) {
  Data <- read.csv("E:/sample_FCS/test_data.csv", check.names = FALSE)
  output$raw_data <- renderUI({
    selectInput("raw_data", h4("channel"), choices = colnames(Data))
  })
  output$transform_type <- renderUI({
    selectInput("transform_type", h4("transform function"), choices = c("asinh","biexponential","hyperlog"))
  })
  output$hyperparameters <- renderUI({
    if(!is.null(input$transform_type)) {
      if(input$transform_type=="asinh") {
        fluidRow(column(12,numericInput("cofactor", h5("cofactor"), value = 5)))
      } else if(input$transform_type=="biexponential") {
        fluidRow(column(12,sliderInput("biexp_pos", h5("biexponential pos"), min = 0.1, max = 10, value = 4.5, step = 0.1),
                        sliderInput("biexp_neg", h5("biexponential neg"), min = 0, max = 10, value = 0, step = 0.1),
                        sliderInput("biexp_width", h5("biexponential width"), min = -20, max = 20, value = -10, step = 0.1)
        ))
      } else if(input$transform_type=="hyperlog") {
        fluidRow(column(12,numericInput("hyperlog_T", h5("hyperlog T"), value = 100000),
                        sliderInput("hyperlog_M", h5("hyperlog M"), min = 0.1, max = 10, value = 5),
                        sliderInput("hyperlog_W", h5("hyperlog W"), min = 0.001, max = 1, value = 0.01),
                        sliderInput("hyperlog_A", h5("hyperlog A"), min = 0.1, max = 10, value = 2)
        ))
      }
    }
  })
  output$densPlot <- renderPlot({
    if(!is.null(input$raw_data)) {
      plot(density(Data[,input$raw_data]), main = input$raw_data, xlab = "intensity", ylab = "density")
    }
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)

