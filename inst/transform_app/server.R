#appR

# Define server logic ----
# server <- function(input, output) {
function(input,output) {
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
  output$init_input <- renderUI({
    selectInput("init_input", h5("UMAP init"), choices = c("spectral","normlaplacian","random","lvrandom","laplacian",
                                                           "pca","spca","agspectral"), selected = "spca")
  })
  output$min_dist_input <- renderUI({
    numericInput("min_dist_input", h5("UMAP min dist (0.0001 - 3)"), value = 0.1, min = 0.0001, max = 3)
  })
  output$spread_input <- renderUI({
    numericInput("spread_input", h5("UMAP spread (0.25 - 3)"), value = 1, min = 0.25, max = 3)
  })
  output$neighbors_input <- renderUI({
    numericInput("neighbors_input", h5("UMAP neighbors (5 - 200)"), value = 30, min = 5, max = 200)
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
  # output$transform_text <- renderText({
  #   chosen_transf()
  # })
  output$warning_message <- renderText({
    "Clicking finalize will lock in the transform values. Are you sure you want to continue?"
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
  observeEvent(input$apply_transform_all, {
    # insert_channel <- gsub("\\:((linear|asinh|biexponential|hyperlog)|(linear|asinh|biexponential|hyperlog).+$)","",chosen_transf())
    # col_index <- which(colnames(param_settings$reactive_data)==insert_channel)
    if(input$transform_type=="linear") { # if transform isn't edited currently, values default to linear, table should reflect this, consider changing default value to asinh
      param_settings$reactive_data <- param_df
      param_settings$reactive_data[1,] <- rep("linear",ncol(param_settings$reactive_data))
    } else if(input$transform_type=="asinh") {
      param_settings$reactive_data <- param_df
      param_settings$reactive_data[1,] <- rep("asinh",ncol(param_settings$reactive_data))
      param_settings$reactive_data[which(row.names(param_settings$reactive_data)=="asinh_cofactor"),] <- rep(input$cofactor,ncol(param_settings$reactive_data))
    } else if(input$transform_type=="biexponential") {
      param_settings$reactive_data <- param_df
      param_settings$reactive_data[1,] <- rep("biexponential",ncol(param_settings$reactive_data))
      param_settings$reactive_data[which(row.names(param_settings$reactive_data)=="biexp_pos"),] <- rep(input$biexp_pos,ncol(param_settings$reactive_data))
      param_settings$reactive_data[which(row.names(param_settings$reactive_data)=="biexp_neg"),] <- rep(input$biexp_neg,ncol(param_settings$reactive_data))
      param_settings$reactive_data[which(row.names(param_settings$reactive_data)=="biexp_width"),] <- rep(input$biexp_width,ncol(param_settings$reactive_data))
    } else if(input$transform_type=="hyperlog") {
      param_settings$reactive_data <- param_df
      param_settings$reactive_data[1,] <- rep("hyperlog",ncol(param_settings$reactive_data))
      param_settings$reactive_data[which(row.names(param_settings$reactive_data)=="hyperlog_T"),] <- rep(input$hyperlog_T,ncol(param_settings$reactive_data))
      param_settings$reactive_data[which(row.names(param_settings$reactive_data)=="hyperlog_M"),] <- rep(input$hyperlog_M,ncol(param_settings$reactive_data))
      param_settings$reactive_data[which(row.names(param_settings$reactive_data)=="hyperlog_W"),] <- rep(input$hyperlog_W,ncol(param_settings$reactive_data))
      param_settings$reactive_data[which(row.names(param_settings$reactive_data)=="hyperlog_A"),] <- rep(input$hyperlog_A,ncol(param_settings$reactive_data))
    }
  })
  twoD_plot <- eventReactive(input$apply_plot, {
    if(!is.null(input$channel)) {
      # x_channel <- input$channel_x; y_channel <- input$channel_y
      transform_options <- c("linear","asinh","biexponential","hyperlog")
      xcol <- which(colnames(param_settings$reactive_data)==input$channel_x); ycol <- which(colnames(param_settings$reactive_data)==input$channel_y)
      in_data_x <- Data[,xcol]; in_data_y <- Data[,ycol]
      in_hyperlog_x <- matrix(Data[,xcol]); colnames(in_hyperlog_x) <- "tmp_x"
      in_hyperlog_y <- matrix(Data[,ycol]); colnames(in_hyperlog_y) <- "tmp_y"
      get_x_transf <- param_settings$reactive_data[1,xcol]; get_y_transf <- param_settings$reactive_data[1,ycol]
      transf_settings_x <- param_settings$reactive_data[-1,xcol]; transf_settings_y <- param_settings$reactive_data[-1,ycol]
      if(all(!is.na(get_x_transf),!is.na(get_y_transf))) {
        if(get_x_transf=="linear") {
          x_data <- in_data_x
        } else if(get_x_transf=="asinh") {
          x_data <- base::asinh(in_data_x/as.numeric(transf_settings_x[1]))
        } else if(get_x_transf=="biexponential") {
          biexp_fun_x <- flowWorkspace::flowjo_biexp(pos = as.numeric(transf_settings_x[2]),
                                                     neg = as.numeric(transf_settings_x[3]),
                                                     widthBasis = (10^as.numeric(transf_settings_x[4]))*-1)
          x_data <- biexp_fun_x(in_data_x)
        } else if(get_x_transf=="hyperlog") {
          # in the future, make apply transform button unclickable if hyperlogtGml2 function is invalid, see error messages when certain values are selected
          # consider https://stackoverflow.com/questions/40621393/disabling-buttons-in-shiny
          hyperlog_fun_x <- flowCore::hyperlogtGml2(parameters = "tmp_x", T = as.numeric(transf_settings_x[5]), M = as.numeric(transf_settings_x[6]),
                                                    W = as.numeric(transf_settings_x[7]), A = as.numeric(transf_settings_x[8]))
          x_data <- eval(hyperlog_fun_x)(in_hyperlog_x)
        }
        if(get_y_transf=="linear") {
          y_data <- in_data_y
        } else if(get_y_transf=="asinh") {
          y_data <- base::asinh(in_data_y/as.numeric(transf_settings_y[1]))
        } else if(get_y_transf=="biexponential") {
          biexp_fun_y <- flowWorkspace::flowjo_biexp(pos = as.numeric(transf_settings_y[2]),
                                                     neg = as.numeric(transf_settings_y[3]),
                                                     widthBasis = (10^as.numeric(transf_settings_y[4]))*-1)
          y_data <- biexp_fun_y(in_data_y)
        } else if(get_y_transf=="hyperlog") {
          # in the future, make apply transform button unclickable if hyperlogtGml2 function is invalid, see error messages when certain values are selected
          # consider https://stackoverflow.com/questions/40621393/disabling-buttons-in-shiny
          hyperlog_fun_y <- flowCore::hyperlogtGml2(parameters = "tmp_y", T = as.numeric(transf_settings_y[5]), M = as.numeric(transf_settings_y[6]),
                                                    W = as.numeric(transf_settings_y[7]), A = as.numeric(transf_settings_y[8]))
          y_data <- eval(hyperlog_fun_y)(in_hyperlog_y)
        }
        plot(x = x_data, y = y_data, pch = 19, col = scales::alpha("black", 0.2),
             xlab = colnames(param_settings$reactive_data)[xcol],
             ylab = colnames(param_settings$reactive_data)[ycol], cex = 0.7)
        } else {
          plot.new()
        }
    }
  })
  output$render_2d_plot <- renderPlot({
    twoD_plot()
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
      return(map <- uwot::umap(X = Data_dynamic, init = input$init_input, min_dist = input$min_dist_input,
                               n_threads = parallel::detectCores(), verbose = TRUE,
                               n_neighbors = input$neighbors_input, spread = input$spread_input))
    }
  })
  output$sample_umap <- renderPlot({ # allow color by channel
    plot_data <- calc_umap()
    if(is.null(plot_data)) {
      plot.new()
    } else {
      plot(plot_data, pch = 19, col = scales::alpha("black",0.1), cex = 0.7, asp = 1, xlab = "UMAP1", ylab = "UMAP2")
    }
  })
  finalize <- eventReactive(input$finalize_transform, {
    # write transform choices to file, include a line to remove files in temp directory before exiting fcj_join function
    transform_algo_chosen <- param_settings$reactive_data[1,]
    if(sum(is.na(transform_algo_chosen))!=0) {
      showModal(modalDialog(
        title = "Error: no transform selected for one or more channels",
        "Check 'Selected Transforms' tab. Algo column should not have any NA values.", easyClose = TRUE
      ))
    } else {
      write.csv(param_settings$reactive_data, file = paste0(system.file(package = "FCSimple"),"/temp_files/tmp_transform_values.csv"), row.names = FALSE)
      showModal(modalDialog(
        title = "Transformation complete!",
        "Close app by clicking the X in the top right corner.", easyClose = FALSE)
      )
    }
    # placeholder value
    # this function defines what happens when the finalize button is clicked
    # transform choices should be passed along and used to apply the transforms to the user's data set
    # the app should close
  })
  observeEvent(input$finalize_transform, {
    finalize()
  })
  output$parameter_df <- renderTable(t(param_settings$reactive_data), rownames = TRUE)
}
