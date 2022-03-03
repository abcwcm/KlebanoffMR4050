#' Add clonotype information to the colData of a SCE object
#' @param sce.object SCE object
#' @param clonotype_ino data.frame with the clonotype information for a given sample
#'
#' @return data.frame of the colData merged with the clonotype information
#'
prep_cell_info <- function(sce.object, clonotype_info = clono_info_wt){
    if(!("Barcode" %in% names(colData(sce.object)) & "barcode" %in% names(clonotype_info))){
        stop("Make sure that the sce.object has a colData entry named 'Barcode' and the clonotype info has a column named 'barcode'.")
    }

    ci <- colData(sce.object)
    ci <- merge(as.data.frame(ci),
        as.data.frame(clonotype_info), # add clonotype info
        by.x = "Barcode", by.y = "barcode", all.x = TRUE)
    row.names(ci) <- ci$cell
    return(ci)
}

#' Add chromosome information to the SCE object for each gene
get_mito <- function(sce.object){

  if("EnsDb.Hsapiens.v86" %in% rownames(installed.packages())){
    library(EnsDb.Hsapiens.v86)
    gn_location <- mapIds(EnsDb.Hsapiens.v86,
        keys = rownames(sce.object), column = "SEQNAME",
        keytype = "GENEID")
    rowData(sce.object)$chr <- paste0("chr", gn_location)
  }else{
    rowData(sce.object)$chr <- ifelse(grepl("MT-", rowData(sce)$gene_symbol), "chrMT", NA)
  }
    return(sce.object)
}

#' Add scater's cell QC to the SCE
#'
#' @param sce.object SCE object
#' @param is_mito vector of TRUE/FALSE indicating which rows correspond to mito-
#' chondrial genes
add_cell_qc <- function(sce.object, is_mito){
    if(length(is_mito) != nrow(sce.object)){stop("The vector indicating mitochondrial genes should have the same length as the number of rows of the sce.object.")}

    cell.qc <- scater::perCellQCMetrics(sce.object,
        subsets=list(mitochondrial=is_mito))
    colData(sce.object) <- cell.qc[, c("sum", "detected",
        "subsets_mitochondrial_sum", "subsets_mitochondrial_detected",
        "subsets_mitochondrial_percent")] %>% cbind(colData(sce.object), . )
    colData(sce.object)$log10_total_features <- log10(sce.object$detected)
    return(sce.object)
}
