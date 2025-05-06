fcs_pca <- function(fcs_join_obj, pca_method = c("prcomp"), num_pc = NULL)
{
  require(stringr)
  require(ggplot2)

  if('batch_correction' %in% names(fcs_join_obj)) {
    cordat <- TRUE
    obj_data <- as.matrix(fcs_join_obj[['batch_correction']][['data']])
    print("batch_correction found in fcs_join_obj list. Using fcs_join_obj[['batch_correction']][['data']] for clustering.")
  } else {
    cordat <- FALSE
    obj_data <- as.matrix(fcs_join_obj[["data"]])
    print("batch_correction not found in fcs_join_obj. Using fcs_join_obj[['data']] for clustering.")
  }
  if(pca_method[1]=="prcomp") {
    set.seed(123) # not really needed
    PCA <- prcomp(x = obj_data)
    reset_seed <- sample(1:5, size = 1)
    rm(set.seed); rm(reset_seed)
    if(!is.null(num_pc)) {
      if(num_pc<1) {
        stop("'num_pc must be either NULL or a whole number greater than 0")
      }
      if(num_pc %% 1 != 0) {
        stop("'num_pc must be either NULL or a whole number greater than 0")
      }
      cvar <- cumsum(PCA$sdev^2 / sum(PCA$sdev^2))
      npc <- num_pc
    } else {
      cvar <- cumsum(PCA$sdev^2 / sum(PCA$sdev^2))
      cvar_inv <- 1 - cvar
      plot(cvar_inv, pch = 19, cex = 1.5)
      npc <- readline(prompt = "How many PCs should be used? Enter either a number or 'drawN' where N is the location on the x-axis where a vertical line will be drawn: ")
      while(grepl(pattern = "draw", x = npc)) {
        npc_draw <- as.numeric(stringr::str_extract(string = npc, pattern = "[0-9]+$"))
        if(npc_draw %% 1 != 0) {
          stop("N must be a whole number greater than 0")
        }
        print(paste0("user input read as: ",npc,". Drawing line at x = ",npc_draw))
        plot(cvar_inv, pch = 19, cex = 1.5); abline(v = npc_draw, lwd = 2, lty = 6, col = "red", main = paste0(npc_draw," PC"))
        npc <- readline(prompt = "How many PCs should be used? Enter either a number or 'drawN' where N is the location on the x-axis where a vertical line will be drawn: ")
        dev.off()
      }
      npc <- as.numeric(npc)
    }
    pc_data <- PCA$x[,c(1:npc)]
  }

  cvar_d <- data.frame(PC = 1:length(cvar), cumVar = cvar); cvar_d_npc <- data.frame(PCx = as.numeric(npc), PCy = cvar_d$cumVar[which(cvar_d$PC==npc)])
  cvar_plt <- ggplot() +
    geom_point(data = cvar_d, mapping = aes(x = PC, y = cumVar), pch = 21, size = 2.25, fill = "white") +
    geom_point(data = cvar_d_npc, mapping = aes(x = PCx, y = PCy), pch = 21, size = 3, fill = "red") +
    ggtitle(paste0(npc," PCs")) +
    xlab("PC") + ylab("cumulative variance") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 14))

  print("storing PCA information in $pca")
  fcs_join_obj[['pca']] <- list(pca_data = pca_data, list(PCs = npc, cumulative_variance = cvar, pca_method = pca_method, elbow_plot = cvar_plt))
  if(!'object_history' %in% names(fcs_join_obj)) {
    print("Consider running FCSimple::fcs_update() on the object.")
  } else {
    fcs_join_obj[['object_history']] <- append(fcs_join_obj[['object_history']], paste0("pca on ",ifelse(cordat,'corrected','uncorrected')," data: ",Sys.time()))
  }
  return(fcs_join_obj)
}
