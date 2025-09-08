#' @title Principal Component Analysis of Flow Cytometry Data
#'
#' @description
#'   Performs PCA on expression data (using batch‐corrected data if available)
#'   and stores the selected principal components, variance explained, and an
#'   elbow plot. If `num_pc` is `NULL`, shows an interactive plot of
#'   (1 – cumulative variance) vs. PC and prompts the user to enter the number
#'   of PCs to retain.
#'
#' @param fcs_join_obj
#'   A list returned by FCSimple::fcs_join(), optionally augmented by
#'   fcs_batch_correction, containing at minimum:
#'   - `data`: numeric matrix of events × channels
#'   - if batch correction was applied, `batch_correction$data`
#'
#' @param pca_method
#'   Character; PCA method to use. Currently only `"prcomp"` (default).
#'
#' @param num_pc
#'   Integer or `NULL` (default); number of principal components to extract.
#'   If `NULL`, the function displays a plot of 1 – cumulative variance vs. PC,
#'   prompts for input (e.g. a number or `"drawN"` to draw a vertical line at
#'   PC N), and repeats until a valid integer is provided.
#'
#' @details
#'   The function:
#'   1. Detects and uses `batch_correction$data` if present.
#'   2. Runs `stats::prcomp` on the selected data with a fixed seed.
#'   3. Computes cumulative variance explained (`cumVar`) and its complement.
#'   4. If `num_pc` is `NULL`, displays the inverse‐variance plot, captures
#'      user input, and highlights the chosen PC on the plot.
#'   5. Extracts the first `num_pc` principal components into `pca_data`.
#'   6. Builds an elbow plot (`ggplot2`) marking the selected PCs.
#'   7. Appends a message to `object_history` noting whether PCA was run on
#'      corrected or uncorrected data and the timestamp.
#'
#' @return
#'   The original `fcs_join_obj` augmented with a new element `pca`, a list with:
#'   - `pca_data`: numeric matrix (events × selected PCs)
#'   - `PCs`: integer number of PCs retained
#'   - `cumulative_variance`: numeric vector of cumulative variance explained
#'   - `pca_method`: character string of the method used
#'   - `elbow_plot`: a ggplot2 object illustrating cumulative variance and
#'                   highlighting the chosen PC count
#'
#' @examples
#' \dontrun{
#'   files   <- list(ff1, ff2)
#'   joined  <- FCSimple::fcs_join(files)
#'
#'   # Interactive PCA: choose number of PCs via prompt
#'   pca_obj <- FCSimple::fcs_pca(joined)
#'
#'   # Run PCA and keep first 5 PCs directly
#'   pca5 <- FCSimple::fcs_pca(joined, num_pc = 5)
#' }
#'
#' @seealso
#'   stats::prcomp, ggplot2::ggplot, FCSimple::fcs_join,
#'   FCSimple::fcs_batch_correction
#'
#' @importFrom stats prcomp
#' @importFrom stringr str_extract
#' @importFrom ggplot2 ggplot aes geom_point ggtitle xlab ylab theme_bw theme element_text
#' @export
fcs_pca <- function(fcs_join_obj, pca_method = c("prcomp"), num_pc = NULL, apply_scaling = TRUE)
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
    if(apply_scaling) {
      PCA <- prcomp(x = obj_data)
    } else {
      PCA <- prcomp(x = obj_data, center = FALSE, scale. = FALSE)
    }
    reset_seed <- sample(1:5, size = 1)
    rm(reset_seed)
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
  fcs_join_obj[['pca']] <- list(pca_data = pc_data, list(PCs = npc, cumulative_variance = cvar, pca_method = pca_method, elbow_plot = cvar_plt))
  if(!'object_history' %in% names(fcs_join_obj)) {
    print("Consider running FCSimple::fcs_update() on the object.")
  } else {
    fcs_join_obj[['object_history']] <- append(fcs_join_obj[['object_history']], paste0("pca on ",ifelse(cordat,'corrected','uncorrected')," data: ",Sys.time()))
  }
  return(fcs_join_obj)
}
