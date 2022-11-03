#appR
library(shiny)
library(shinythemes)
library(flowWorkspace)
library(flowCore)

# Define UI ----
ui <- fluidPage(theme = shinytheme("yeti"),
                navbarPage(
                  "FCSimple transformer",
                  tabPanel("transform",
                           sidebarPanel(
                             uiOutput("channel"),
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
  output$channel <- renderUI({
    selectInput("channel", h4("channel"), choices = colnames(Data))
  })
  output$transform_type <- renderUI({
    selectInput("transform_type", h4("transform function"), choices = c("linear","asinh","biexponential","hyperlog"), selected = "linear")
  })
  output$hyperparameters <- renderUI({
    if(!is.null(input$transform_type)) {
      if(input$transform_type=="asinh") {
        fluidRow(column(12,numericInput("cofactor", h5("cofactor"), value = 5)))
      } else if(input$transform_type=="biexponential") {
        fluidRow(column(12,sliderInput("biexp_pos", h5("biexponential pos"), min = 2, max = 7, value = 4, step = 0.05),
                        sliderInput("biexp_neg", h5("biexponential neg"), min = 0, max = 1, value = 0, step = 0.1),
                        sliderInput("biexp_width", h5("biexponential width"), min = 0, max = 3, value = 0, step = 0.05)
        ))
      } else if(input$transform_type=="hyperlog") {
        fluidRow(column(12,numericInput("hyperlog_T", h5("hyperlog T"), value = 100000),
                        sliderInput("hyperlog_M", h5("hyperlog M"), min = 2, max = 7, value = 4, step = 0.1),
                        sliderInput("hyperlog_W", h5("hyperlog W"), min = 0.01, max = 3, value = 0.01, step = 0.1),
                        sliderInput("hyperlog_A", h5("hyperlog A"), min = 0, max = 3, value = 2, step = 0.1)
        ))
      }
    }
  })
  output$densPlot <- renderPlot({
    if(!is.null(input$channel)) {
      in_data <- Data[,input$channel]
      in_hyperlog <- matrix(Data[,input$channel]); colnames(in_hyperlog) <- "tmp"
      if(input$transform_type=="linear") {
        plot(density(in_data), main = input$channel, xlab = "intensity", ylab = "density")
      } else if(input$transform_type=="asinh") {
        get_cofactor <- input$cofactor
        plot(density(base::asinh(in_data/get_cofactor)), main = input$channel, xlab = "intensity", ylab = "density")
      } else if(input$transform_type=="biexponential") {
        biexp_fun <- flowWorkspace::flowjo_biexp(pos = input$biexp_pos,
                                                 neg = input$biexp_neg,
                                                 widthBasis = (10^input$biexp_width)*-1)
        plot(density(biexp_fun(in_data)), main = input$channel, xlab = "intensity", ylab = "density")
      } else if(input$transform_type=="hyperlog") {
        hyperlog_fun <- flowCore::hyperlogtGml2(parameters = "tmp", T = input$hyperlog_T, M = input$hyperlog_M,
                                                W = input$hyperlog_W, A = input$hyperlog_A)
        plot(density(eval(hyperlog_fun)(in_hyperlog)), main = input$channel, xlab = "intensity", ylab = "density")
      }
    }
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)

