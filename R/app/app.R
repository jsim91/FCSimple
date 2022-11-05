#appR
library(shiny)
library(shinythemes)
library(flowWorkspace)
library(flowCore)
library(scales)
library(flowCore)
library(parallel)
library(uwot)

# to do: store final transform values, save out transform parameters;
# allow user to update these values every time the button is pressed; finally, apply transform to data set
# and produce an output with the transform values used, allow an option to parse this output for future
# transform setting on files using identical panels, ensures consistency and streamlines workflow consider
# a json format or a csv file with columns channel names, and rows parameters, perhaps a row for alternate
# names, such as fluorochrome + marker or metal + marker to allow for more robust matching of future panels

# Define UI ----
ui <- fluidPage(theme = shinytheme("yeti"),
                navbarPage(
                  "FCSimple transformer",
                  tabPanel("transform",
                           sidebarPanel(
                             uiOutput("channel"),
                             uiOutput("transform_type"),
                             uiOutput("hyperparameters"),
                             hr(),
                             actionButton(inputId = "apply_transform", label = "apply transform"),
                             hr(),
                             uiOutput("channel_x"),
                             uiOutput("channel_y"),
                             actionButton(inputId = "apply_plot", label = "update 2d plot")
                           ),
                           mainPanel(
                             plotOutput("densPlot"),
                             textOutput("dens_hyperlog_error"),
                             textOutput("dens_asinh_error"),
                             hr(),
                             plotOutput("biax"),
                             hr(),
                             uiOutput("transform_text")
                           )),
                  tabPanel("selected transforms",
                           tableOutput("parameter_df")),
                  tabPanel("UMAP preview",
                           sidebarPanel(
                             actionButton(inputId = "umap_button", label = "calculate sample UMAP"),
                           ),
                           mainPanel(
                             plotOutput("sample_umap")
                           ))
                           # actionButton(inputId = "umap_button", label = "calculate sample UMAP"))
                )
)

# Define server logic ----
server <- function(input, output) {
  Data <- read.csv("E:/sample_FCS/test_data.csv", check.names = FALSE)
  Data_dynamic <- Data
  param_df <- as.data.frame(matrix(data = NA,nrow=9,ncol=ncol(Data)))
  colnames(param_df) <- colnames(Data)
  row.names(param_df) <- c("algo","asinh_cofactor","biexp_pos","biexp_neg","biexp_width",
                           "hyperlog_T","hyperlog_M","hyperlog_W","hyperlog_A")
  param_settings <- reactiveValues(reactive_data = NULL)
  param_settings$reactive_data <- param_df
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
  output$channel_x <- renderUI({
    selectInput("channel_x", h5("x-axis channel"), choices = colnames(Data))
  })
  output$channel_y <- renderUI({
    selectInput("channel_y", h5("y-axis channel"), choices = colnames(Data))
  })
  output$densPlot <- renderPlot({
    if(!is.null(input$channel)) {
      in_data <- Data[,input$channel]
      in_hyperlog <- matrix(Data[,input$channel]); colnames(in_hyperlog) <- "tmp"
      if(input$transform_type=="linear") {
        plot(density(in_data), main = input$channel, xlab = "intensity", ylab = "density")
      } else if(input$transform_type=="asinh") {
        if(all(!is.null(input$cofactor),input$cofactor!=0)) {
          get_cof <- input$cofactor
          plot(density(base::asinh(in_data/get_cof )), main = input$channel, xlab = "intensity", ylab = "density")
        }
      } else if(input$transform_type=="biexponential") {
        if(all(!is.null(input$biexp_pos),!is.null(input$biexp_neg),!is.null((10^input$biexp_width)*-1))) {
          biexp_fun <- flowWorkspace::flowjo_biexp(pos = input$biexp_pos,
                                                   neg = input$biexp_neg,
                                                   widthBasis = (10^input$biexp_width)*-1)
          plot(density(biexp_fun(in_data)), main = input$channel, xlab = "intensity", ylab = "density")
        }
      } else if(input$transform_type=="hyperlog") {
        if(all(!is.null(input$hyperlog_T),!is.null(input$hyperlog_M),!is.null(input$hyperlog_W),!is.null(input$hyperlog_A))) {
          # in the future, make apply transform button unclickable if hyperlogtGml2 function is invalid, see error messages when certain values are selected
          # consider https://stackoverflow.com/questions/40621393/disabling-buttons-in-shiny
          if(all(input$hyperlog_W<(input$hyperlog_M/2), input$hyperlog_A<=(input$hyperlog_M - (input$hyperlog_W * 2)))) {
            hyperlog_fun <- flowCore::hyperlogtGml2(parameters = "tmp", T = input$hyperlog_T, M = input$hyperlog_M,
                                                    W = input$hyperlog_W, A = input$hyperlog_A)
            plot(density(eval(hyperlog_fun)(in_hyperlog)), main = input$channel, xlab = "intensity", ylab = "density")
          }
        }
      }
    }
  })
  output$dens_hyperlog_error <- renderText({
    if(all(!is.null(input$channel),input$transform_type=="hyperlog")) {
      if(any(input$hyperlog_W>=(input$hyperlog_M/2), input$hyperlog_A>(input$hyperlog_M - (input$hyperlog_W * 2)))) {
        paste0("W must be less than M/2 .. AND .. A must be less than or equal to M - 2W")
      }
    }
  })
  output$dens_asinh_error <- renderText({
    if(all(!is.null(input$channel),input$transform_type=="asinh")) {
      if(all(!is.null(input$cofactor),input$cofactor==0)) {
        "cofactor must not equal 0"
      }
    }
  })
  chosen_transf <- eventReactive(input$apply_transform, {
    if(is.null(input$transform_type)) {
      return("")
    } else {
      tmp_transform_choice <- input$transform_type
      if(input$transform_type=="linear") {
        return(paste0(input$channel,":",tmp_transform_choice))
      } else if(input$transform_type=="asinh") {
        tmp_cof <- input$cofactor
        return(paste0(input$channel,":",tmp_transform_choice,",",tmp_cof))
      } else if(input$transform_type=="biexponential") {
        tmp_pos <- input$biexp_pos
        tmp_neg <- input$biexp_neg
        tmp_wid <- (10^input$biexp_width)*-1
        return(paste0(input$channel,":",tmp_transform_choice,",",tmp_pos,",",tmp_neg,",",tmp_wid))
      } else if(input$transform_type=="hyperlog") {
        tmp_t <- input$hyperlog_T
        tmp_m <- input$hyperlog_M
        tmp_w <- input$hyperlog_W
        tmp_a <- input$hyperlog_A
        return(paste0(input$channel,":",tmp_transform_choice,",",tmp_t,",",tmp_m,",",tmp_w,",",tmp_a))
      }
    }
  })
  output$transform_text <- renderText({
    chosen_transf()
  })
  observeEvent(input$apply_transform, {
    insert_channel <- gsub("\\:((linear|asinh|biexponential|hyperlog)|(linear|asinh|biexponential|hyperlog).+$)","",chosen_transf())
    col_index <- which(colnames(param_settings$reactive_data)==insert_channel)
    if(input$transform_type=="linear") { # if transform isn't edited currently, values default to linear, table should reflect this, consider changing default value to asinh
      param_settings$reactive_data[,col_index] <- c("linear",rep(NA,8))
    } else if(input$transform_type=="asinh") {
      param_settings$reactive_data[,col_index] <- c("asinh",input$cofactor,rep(NA,7))
    } else if(input$transform_type=="biexponential") {
      param_settings$reactive_data[,col_index] <- c("biexponential",NA,input$biexp_pos,input$biexp_neg,input$biexp_width,rep(NA,4))
    } else if(input$transform_type=="hyperlog") {
      param_settings$reactive_data[,col_index] <- c("hyperlog",rep(NA,4),input$hyperlog_T,input$hyperlog_M,input$hyperlog_W,input$hyperlog_A)
    }
  })

  # update_2d should be able to draw from param_settings$reactive_data
  # update_2d <- eventReactive(input$apply_plot, {
  #   if(is.null(input$transform_type)) {
  #     return(NULL) # should return the initial plot based on initial values
  #   } else {
  #     tmp_transform_choice <- input$transform_type
  #     if(input$transform_type=="linear") {
  #       return(paste0(input$channel,": ",tmp_transform_choice))
  #     } else if(input$transform_type=="asinh") {
  #       tmp_cof <- input$cofactor
  #       return(paste0(input$channel,": ",tmp_transform_choice," - ",tmp_cof))
  #     } else if(input$transform_type=="biexponential") {
  #       tmp_pos <- input$biexp_pos
  #       tmp_neg <- input$biexp_neg
  #       tmp_wid <- (10^input$biexp_width)*-1
  #       return(paste0(input$channel,": ",tmp_transform_choice," - ",tmp_pos," - ",tmp_neg," - ",tmp_wid))
  #     } else if(input$transform_type=="hyperlog") {
  #       tmp_t <- input$hyperlog_T
  #       tmp_m <- input$hyperlog_M
  #       tmp_w <- input$hyperlog_W
  #       tmp_a <- input$hyperlog_A
  #       return(paste0(input$channel,": ",tmp_transform_choice," - ",tmp_t," - ",tmp_m," - ",tmp_w," - ",tmp_a))
  #     }
  #   }
  # })
  # output$transform_text <- renderPlot({
  #   update_2d()
  # })

  output$biax <- renderPlot({
    if(!is.null(input$transform_type)) {
      plot(x = Data[,input$channel_x], y = Data[,input$channel_y], pch = 19, col = scales::alpha("black", 0.2),
           xlab = input$channel_x, ylab = input$channel_y, cex = 0.7) # inherit transform chosen and parameter choice/alter the Data directly when transform is applied with button, when button is selected, first take Data back to pre-transform then apply transform
    }
  })
  output$transform_text <- renderText({
    chosen_transf()
  })
  calc_umap <- eventReactive(input$umap_button, {
    if(is.null(input$transform_type)) {
      return(NULL)
    } else {
      set.seed(123)
      Data_dynamic <- Data[sample(1:nrow(Data),size=10000,replace=F),]
      for(i in 1:ncol(param_settings$reactive_data)) {
        if(any(is.na(param_settings$reactive_data[1,i]),param_settings$reactive_data[1,i]=="linear")) {
          Data_dynamic[,i] <- Data_dynamic[,i]
        } else if(param_settings$reactive_data[1,i]=="asinh") {
          Data_dynamic[,i] <- asinh(Data_dynamic[,i]/as.numeric(param_settings$reactive_data[2,i]))
        } else if(param_settings$reactive_data[1,i]=="biexponential") {
          Data_dynamic[,i] <- flowWorkspace::flowjo_biexp(pos = as.numeric(param_settings$reactive_data[3,i]),
                                                          neg = as.numeric(param_settings$reactive_data[4,i]),
                                                          widthBasis = as.numeric(param_settings$reactive_data[5,i]))(Data_dynamic[,i])
        } else if(param_settings$reactive_data[1,i]=="hyperlog") {
          fs_data <- as.data.frame(Data_dynamic[,i]); colnames(fs_data) <- colnames(Data_dynamic)[i]
          tmp_fs <- new("flowFrame",exprs=as.matrix(fs_data))
          transform_FUN <- flowCore::hyperlogtGml2(parameters = flowCore::colnames(tmp_fs)[1],
                                                   T = as.numeric(param_settings$reactive_data[6,i]),
                                                   M = as.numeric(param_settings$reactive_data[7,i]),
                                                   W = as.numeric(param_settings$reactive_data[8,i]),
                                                   A = as.numeric(param_settings$reactive_data[9,i]))
          Data_dynamic[,i] <- eval(transform_FUN)(exprs(tmp_fs))
        }
      }
      return(map <- uwot::umap(X = Data_dynamic, init = "spca", min_dist = 0.1, n_threads = parallel::detectCores(), verbose = TRUE))
    }
  })
  output$sample_umap <- renderPlot({ # allow color by channel
    plot_data <- calc_umap()
    if(is.null(plot_data)) {
      plot.new()
    } else {
      plot(plot_data)
    }
  })
  output$parameter_df <- renderTable(t(param_settings$reactive_data), rownames = TRUE)
}

# Run the app ---- # this will be changed later so that the app is called from fcs_join with runApp
shinyApp(ui = ui, server = server)

