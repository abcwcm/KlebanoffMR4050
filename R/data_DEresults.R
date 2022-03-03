#' Load the list of DE test results (any direction) for all clonotypes present
#' in both conditions ('shared' clonotypes)
#'
#' @details The p- and q-values here represent a two-tailed test for any direction
#' of the logFC.
#' For details on how the DE analysis was done, see the vignette "DE_genes"
#' and the wrapper function \code{\link{run_DE}}.
#'
#' @format Nested list where names correspond to the abbreviated clonotype IDs.
#' For every clonotype, there's a list that contains:
#' \describe{
#' \item{findMarkers_results:}{A SimpleDataFrameList, i.e. the original result
#' of \code{scran::findMarkers()}, but only for the MUT comparison}
#' \item{marker_IDs:}{a data.table with the genes that passed the FDR threshold;
#' if that is NULL, this implies that there were no DEG for that particular
#' clonotype comparing MUT vs WT}
#' }
#' @usage load_DE_results("MR4050")
#' @examples \dontrun{
#' library(KlebanoffMR4050)
#'
#' sce.shared <- load_MR4050shared()
#' sce.shared$antigen <-  factor(gsub("\\..*","",sce.shared$Sample),
#'     levels = c("WT", "MUT"), ordered = TRUE)
#'
#' delist.both <- lapply( unique(sce.shared$id), function(x){
#' run_DE(
#'    sce.shared[, sce.shared$id == x],
#'    group_identifier = "antigen",
#'    direction = "any",
#'    FDR = 0.05, rank = Inf,
#'    comp_name = paste0(x, "_"))
#'    })
#' names(delist.both) <- unique(sce.shared$id)
#' }
#'
#'@export
#'
load_DE_results <- function(sample = "MR4050"){
    fn <- system.file("extdata", "delist.both",
        package = paste0("Klebanoff", sample), mustWork = TRUE)

    fin <- read.table(fn, stringsAsFactor = FALSE)[[1]]
    load_data_from_Box(fin, load_rds = FALSE)
}
