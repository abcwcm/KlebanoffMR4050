#' Identify genes to be kept
#'
#' @description This function is a wrapper around scater::perFeatureQCMetrics
#' which is run on each sample (defined in the sce.object's "Sample" entry)
#' individually.
#'
#' @param sce.object SingelCellExperiment object
#' @param min.cells minimum number of cells that must be expressing a given gene.
#' (Min. number of cells per Sample, that is). Default: 5
#' @param return_df whether the results of scater's QC should be returned in
#' a data.frame. Default: FALSE
#'
#' @return vector of TRUE/FALSE that indicates which genes passed the filters
#' and which ones did not
#'
#' @examples \dontrun{
#' keep_genes <-  do_gene_qc(sce.filt)
#' sce.filt[keep_genes,]
#' }
do_gene_qc <- function(sce.object, min.cells = 5, return_df = FALSE){

    if(!any(grepl("WT", sce.object$Sample))){stop("Seems like there are no WT samples in your data.")}
    if(!any(grepl("MUT", sce.object$Sample))){stop("Seems like there are no MUT samples in your data.")}
    if(!("Sample" %in% names(colData(sce.object)))){stop("There must be an entry named 'Sample' in your colData.")}

    ## define the subsets
    sbs.list <- lapply(unique(sce.object$Sample),
        function(x) colnames(sce.object[,sce.object$Sample == x]) )
    names(sbs.list) <- unique(sce.object$Sample)

    ## run scater's QC function
    feat.qc <- scater::perFeatureQCMetrics(sce.object,
        subsets = sbs.list, flatten = TRUE)
    feat.qc$Gene <- rownames(sce.object)

    ## define the threshold for gene exclusion based on percentage of
    ## cells with non-zero expression (as stored in "detected" by scater)
    thresh_wt <- min.cells/ncol(sce.object[,grepl("WT", sce.object$Sample)])
    thresh_mut <- min.cells/ncol(sce.object[,grepl("MUT", sce.object$Sample)])

    ## identify the columns of interest from the scater output
    col_mut <- grep("subsets.*MUT.*detected$", names(feat.qc), value=TRUE)
    col_wt <- grep("subsets.*WT.*detected$", names(feat.qc), value=TRUE)

    ## define the Booleans
    sel_mut <- feat.qc[[col_mut]] >= thresh_mut
    sel_wt <- feat.qc[[col_wt]] >= thresh_wt

    ## define the final gene set
    kp_gns  <- sel_mut & sel_wt

    if(return_df){
        hkg <- unique(c(grep("^mt-", rownames(sce.object),
            value=TRUE, ignore.case=TRUE),  # mitochondrial genes
            grep("^Rp[sl]", rownames(sce.object),
                value=TRUE, ignore.case=TRUE))) # ribosomal genes
        return(list(keep_genes = kp_gns, gene_qc = feat.qc, hk_genes = hkg))
        }else{
            return(kp_gns)
        }
}

#' Calculate cell cycle status per cell
#'
#' @description Wrapper function around \code{scran::cyclone()}. Assumes the
#' data is from a human data set.
#'
#' @details Users should not filter out low-abundance genes before applying
#' cyclone. Even if a gene is not expressed in any cell, it may still be useful
#' for classification if it is phase-specific. Its lack of expression relative
#' to other genes will still yield informative pairs, and filtering them out
#' would reduce power.
#'
#' @param sce.object SCE object
#' @param print_results if set to TRUE (default), the total number of cells
#' per cell cycle phase is printed to the screen
#' @return A list with \code{phases} (character vector of the predicted phase),
#' \code{scores} (data frame with the numeric phase scores for each phase and cell)
#' and \code{normalized.scores} (data frame with the row-normalized scores).
#' @seealso \code{\link{scran::cyclone()}}
get_cc_info <- function(sce.object, print_results = TRUE){
    set.seed(123)
    hg.prs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
    p = bpstart(MulticoreParam(12))
    cc.out <- scran::cyclone(sce.object, pairs=hg.prs, BPPARAM=p)
    names(cc.out$phases) <- colnames(sce.object)
    row.names(cc.out$scores) <- colnames(sce.object)
    row.names(cc.out$normalized.scores) <- colnames(sce.object)

    if(print_results){
        table(cc.out$phases)

        cols <- c("MUT", "WT")
        sapply(cols, function(x){
            table(cc.out$phases[grep(x, names(cc.out$phases), ignore.case = TRUE)])} )
    }

    return(cc.out)
}