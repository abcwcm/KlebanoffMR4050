#' Extract marker gene results from findMarkers into a data.table
#'
#' @param sce.object SCE object
#' @param marker_search_result result of \code{scran::findMarkers}, e.g. l.mks$HFDvsWT;
#' typically a list where the names indicate the comparison. This function here
#' assumes that findMarkers was run with \code{direction = "up"}.
#' @param FDR_thresh Default: 0.01
#' @param rank_thresh maximum rank that the extracted genes must have in addition
#' to \code{FDR_thresh}
#' @param logFC_thresh minimum logFC that the extracted genes must have (default: 0)
#' @param comp_name add a name for the comparison that underlies the marker search
#' result, e.g. "DB_clst"; if NULL (default), just the name from the \code{marker_search_result}
#' will be used.
#'
#' @author Friederike Dündar
#'
#' @import data.table
#' @import magrittr
#'
#' @return data.table with three columns:
#' \itemize{
#' \item \code{gene_symbol}
#' \item \code{up_in}: in which condition(s) the gene was upregulated
#' \item \code{classify}: whether the gene was found in just one comparison or
#' multiple ones
#' }
#'
#' @examples \dontrun{
#' gns_HFDvsWT <- extract_markers(sf, marker_search_result = l.mks$HFDvsWT,
#'                                FDR_thresh = 0.01, rank_thresh = 35)
#' }
#'
#' @seealso \code{\link{plot_marker_heatmap}}
#'
#'
#' @export
#'
extract_markers <- function(sce.object, marker_search_result, FDR_thresh = 0.01,
  rank_thresh = 35, logFC_thresh = 0,
  comp_name = NULL){

  if(!all(is.numeric(FDR_thresh) & is.numeric(logFC_thresh) & is.numeric(rank_thresh))){
    stop("One of your filters (FDR, rank, logFC) is not numeric.")
  }

  ## extract genes from the marker list-----------------------------------------
  gns <- lapply(seq_along(marker_search_result), function(i){
    ## use those genes that pass the defined thresholds
    tmp <- as.matrix(subset(marker_search_result[[i]],
      FDR <= FDR_thresh & Top <= rank_thresh))

    ## at least one logFC value (there may be multiple per gene) must pass the threshold
    tmp2 <- tmp[,grep("logFC",colnames(tmp), value = TRUE), drop=FALSE]
    keep <- apply(tmp2, 1, function(x) any(x > logFC_thresh))
    goi <- rownames(tmp[keep,])

    ## turn the results into a data.table
    if(!is.null(goi)){
      out <- data.table(gene_symbol = goi)
      out$cell_type <- paste0(comp_name, names(marker_search_result[i]))
      return(out)
    }else{
      message(paste("No genes passed the thresholds for", names(marker_search_result[i]) ))
    }
  }) %>% rbindlist

  if(dim(gns)[1] > 0){
    gns <- gns[,paste(cell_type, collapse = ","), gene_symbol]
    setnames(gns,"V1", "up_in")
    gns[, classify := ifelse(grepl(",", up_in), "shared", "unique")][] # [] is needed for the ata.table being printed immediately when the function is used

    # > head(gns)
    #gene_symbol      up_in classify
    #1:       Rps27 DB_clstHFD   unique
    #2:       Rpl41 DB_clstHFD   unique
    #3:       Rps29 DB_clstHFD   unique
    #4:      Malat1 DB_clstHFD   unique
    #5:       Rps28 DB_clstHFD   unique
    #6:       Pcsk2 DB_clstHFD   unique
    return(gns)
  }

}

#' Prepare a data.table in skinny format suitable for ggplotting
#'
#' @param object \code{SCE} object
#' @param exprs_values what type of expression values should be used from the
#' \code{SCE} object (default: "logcounts"); can be more than one name
#' (e.g., \code{names(assay(object))}).
#' @param features vector of feature names (e.g. gene names, HTOs, ADTs)
#' for which the \code{data.table} should be made
#' (optional; default: NULL)
#' @param include_metaData indicate whether the resulting \code{data.table}
#' should also include specified entries of colData(object). Cacn include specific
#' gene or other feature names. Default: NULL.
#'
#'
#' @return \code{data.table} in skinny format with additional columns
#' "cell" and whatever was specified via \code{include_metaData}
#'
#' @examples \dontrun{
#' # get a long data.table with values from all assays
#' exp.dt <- make_long_dt(sce, exprs_values = names(assays(sce)),
#'                       features = c("SHOX2", "MYH6"),
#'                       include_metaData= c("SampleName", "Condition","Replicate"))
#' }
#'
#' @seealso \code{\link{fx.extract_exprs}}
#'
#' @author Friederike Dündar
#'
#' @import data.table
#' @import scater
#'
#' @export
#'
make_long_dt <- function (object, exprs_values = "logcounts", features = NULL,
  include_metaData = NULL){

  ## define genes
  if (!is.null(features)) {
    features <- unique(features)
    genes <- features[features %in% rownames(object)]
    obj <- object[genes, ]
    gn <- genes
  }
  else {
    obj <- object
    gn <- rownames(object)
  }

  ## extract the values for all indicated exprs_value types -------------------
  dt_list <- lapply(exprs_values, function(x){
    fx.extract_exprs(obj, x, gn) } )
  out.dt <- rbindlist(dt_list)

  ## check altExp for remaining features
  if(length(features[!features %in% rownames(object)]) > 0){

    for(EXP in altExpNames(object)){

      obj.alt <- altExp(object, EXP)
      gn.alt <- features[features %in% rownames(obj.alt)]

      if(length(gn.alt) > 0){
        ev <- exprs_values[exprs_values %in% assayNames(obj.alt)]
        out.alt <- rbindlist(lapply(ev, function(x){
          fx.extract_exprs(obj.alt, x, gn.alt)}))
        out.dt <- rbindlist(list(out.dt, out.alt))
      }
    }
  }

  out.dt <- data.table::dcast(data = out.dt, feature_name + cell ~ value_type, value.var = "value")

  ## add colData -----------------------------------------------------
  if(!is.null(include_metaData)){
    found_mt <- FALSE
    varLabs <- names(colData(obj))
    if (any(include_metaData %in% varLabs )) {
      inc_mt <- include_metaData[include_metaData %in%  varLabs]
      meta_info <- as.data.table(colData(obj)[inc_mt])
      meta_info$cell <- colnames(obj)
      out.dt <- meta_info[out.dt, on = "cell"]
      found_mt <- TRUE
    }

    ## add gene values as meta data
    gnNames <- rownames(object)
    if (any(include_metaData %in% gnNames)) {
      inc_mt <- include_metaData[include_metaData %in% gnNames]
      if(!exists("out.dt")){out.dt <- NULL}
      out.dt <- fx.features_as_metaData(object, exprs_values = exprs_values,
        which_feat = inc_mt, meta.dt = out.dt)
      found_mt <- TRUE
    }

    ## add additional features' values as meta data
    for(EXP in altExpNames(object)){
      obj.alt <- altExp(object, EXP)
      gn.alt <- rownames(obj.alt)
      if(any(include_metaData %in% gn.alt)){
        inc_mt <- include_metaData[include_metaData %in% gn.alt]
        if(!exists("out.dt")){out.dt <- NULL}
        ev <- exprs_values[exprs_values %in% assayNames(obj.alt)]
        out.dt <- fx.features_as_metaData(obj.alt, exprs_values = ev,
          which_feat = inc_mt, meta.dt = out.dt)
        found_mt <- TRUE
      }
    }

    if (!found_mt) {
      if (!(all(include_metaData == FALSE))) {
        warning("Not all of the labels indicated by you via include_metadata are part of the sce.object. There will be no metadata returned.")
      }
    }
  }
  return(out.dt)
}

#' Features as meta data
#' @description Extracts expression values for individual features and adds
#' them as a column
#' @param object SCE object
#' @param exprs_values Types of expression values, i.e. which assay
#' to extract from (use \code{assayNames(sce)} if you're not sure
#' which ones to go with
#' @param which_feat feature name(s)
#' @param meta.dt data.table with already extracted metadata
fx.features_as_metaData <- function(object,exprs_values, which_feat, meta.dt = NULL){

  mt_gn <- lapply(exprs_values, function(x) {
    tmout <- as.data.frame(assay(object[unique(which_feat), ], x))
    tmout$feature_name <- row.names(tmout)
    tmout <- as.data.table(tmout)
    tmout$value_type <- x
    return(tmout)
  })

  mt_gn <- melt.data.table(
    rbindlist(mt_gn),
    variable.name = "cell",
    id.vars = c("value_type", "feature_name"))
  mt_gn <- dcast.data.table(mt_gn,
    cell ~ feature_name + value_type, sep = ".")

  ## merge with other metadata
  if(!is.null(meta.dt)){
    outt <- mt_gn[meta.dt, on = "cell"]
  }else{
    outt <- mt_gn
  }

  return(outt)
}

#' Extract expression values from SCE object
#'
#' @description Wrapper function for extracting expression values into a data.table
#' from an SCE object including checking the output.
#'
#' @param object SCE object
#' @param exprs_values which type of \code{assay} values to extract
#' @param genes vector of gene names (subset of rownames)
#'
#' @return molten data.table with one row per gene and cell and the corresponding
#' \code{exprs_values} and their type in a third and fourth column.
#' @seealso \code{\link{make_long_dt}}
#'
#' @import magrittr
#' @import data.table
#'
fx.extract_exprs <- function (object, exprs_values, genes){

  long.dt <- assay(object[unique(genes), ], exprs_values) %>% as.matrix %>%
    data.table
  if (is.null(long.dt)) {
    stop("There were no values to extract. Check the name for exprs_values.")
  }

  long.dt$feature_name <- unique(genes)
  long.dt <- melt(long.dt, id.vars = "feature_name", variable.name = "cell")
  long.dt$value_type <- exprs_values
  return(long.dt)
}

#' Download and locally cache an RDA or RDS file that's stored on a website
#'
#' @description Using the library \code{BiocFileCache}, this function takes as
#' input the direct link to a RDA (or RDS) file, downloads it and stashes it in a
#' local directory from where it can be retrieved (using this exact function) in
#' a later session. I.e. if the data has already been downloaded, it won't be down-
#' loaded again.
#'
#' @details Based on instructions from the BiocFileCache vignette: \url{https://bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html#local-cache-of-an-internet-resource}. This is a slightly more versatile form of
#' the function \code{load_RDSdata_from_Box()} as it also allows to load RDA files
#' directly into the global environment by setting the option \code{load_rds} to FALSE.
#'
#' @param shared_link The URL for the file as a string. NOTE: This should be the
#' \emph{direct} link, i.e. it should end with ".rda" or ".rdata", "rds" or the like
#' @param data_name A customized name that may be used to more easily refer to
#' the cached data (\code{rname} parameter of \code{BiocFileCache} functions),
#' e.g. "SCE2019" -- any kind of string will work, but it's important
#' to make a note of this name for reloading the data next time, otherwise the
#' data set will be downloaded again if there's a mismatch with the name.
#' If set to NULL (default), the \code{rname} in the cache will correspond to the
#' url indicated via \code{shared_link}.
#' @param cache_path If wanted, this allows to specify a path where the data
#' should be cached. if NULL (default), the default setting of \code{BiocFileCache()}
#' will be used.
#' @param check_for_update Indicate whether the data should be downloaded again.
#'
#' @return The final command here is \code{readRDS()}, so you should be assigning
#' the result of this function to an object name in your environment.
#'
#' @examples
#'
#' ## load RDS requires the user to define an object name (here: testingsce)
#' direct_link <- "https://wcm.box.com/shared/static/4887keeu77eskmdotve1sm7n6u7udckg.rds"
#' testingsce <- load_data_from_Box(shared_link = direct_link, data_name = "testRDS")
#'
#' ## loading RDA
#' sl <- "https://wcm.box.com/shared/static/mmf98o464n94fdh78x14sbwo7ts56ukq.rda"
#' load_data_from_Box(sl, load_rds = FALSE) # should add tiny_sce to global env.
#' ls()
#'
#' @export
#'
load_data_from_Box <- function(shared_link = "https://wcm.box.com/shared/static/mmf98o464n94fdh78x14sbwo7ts56ukq.rda",
    data_name = NULL, cache_path = NULL, check_for_update = FALSE, load_rds = TRUE){

    if(grepl("\\.rds$", shared_link, ignore.case = TRUE)){
        load_rds <- TRUE
    }

    suppressPackageStartupMessages({ library(BiocFileCache) })
    load_new <- FALSE

    if(is.null(data_name)){
        data_name <- shared_link
    }

    ## load the cache
    if(is.null(cache_path)){
        bfc <- BiocFileCache(ask=F)
    }else{
        bfc <- BiocFileCache(cache_path, ask=FALSE)
    }

    ## check if url is being tracked
    bfc_res <- bfcquery(bfc, shared_link)

    if(bfccount(bfc_res) == 0L){
        load_new <- TRUE
    }else{
        ## if it is in cache, get path to load
        rid <- bfc_res[bfc_res$fpath == shared_link & bfc_res$rname == data_name,]$rid

        ## double-check we're looking at what we want
        if(length(rid) > 0){
            if(length(rid) > 1){
                message(paste0("Multiple entries found for", data_name, shared_link,". We'll use only the first one."))
                rid <- rid[1]
            }

            ans <- bfcrpath(bfc, rids = rid)

            if(check_for_update){
                ## check to see if the resource needs to be updated
                check <- bfcneedsupdate(bfc, rids = rid)
                if( is.na(check)) upd <- TRUE ## 'check' can be NA if it cannot be determined, choose how to handle
                if(upd) ans <- bfcdownload(bfc, rid = rid)
            }
        }else{
            load_new <- TRUE
        }
    }

    ## if it wasn't in the cache, add
    if(load_new){
        ans <- bfcadd(bfc, rname = data_name, fpath = shared_link)
    }

    message(paste("Cached data here:", ans))

    if(load_rds){
        readRDS(ans)
    }else{
        load(ans, .GlobalEnv)
    }

}


#' Download and locally cache an RDS file that's stored on a website
#'
#' @description Using the library \code{BiocFileCache}, this function takes as
#' input the direct link to a RDS (!) file, downloads it and stashes it in a
#' local directory from where it can be retrieved (using this exact function) in
#' a later session. I.e. if the data has already been downloaded, it won't be down-
#' loaded again.
#'
#' @details Based on instructions from the BiocFileCache vignette: \url{https://bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html#local-cache-of-an-internet-resource}
#'
#' @param shared_link The URL for the file as a string. NOTE: This should be the
#' \emph{direct} link, i.e. it should end with ".rds" (!)
#' @param data_name A customized name that may be used to more easily refer to
#' the cached data (\code{rname} parameter of \code{BiocFileCache} functions),
#' e.g. "SCE2019" -- any kind of string will work, but it's important
#' to make a note of this name for reloading the data next time, otherwise the
#' data set will be downloaded again if there's a mismatch with the name.
#' If set to NULL (default), the \code{rname} in the cache will correspond to the
#' url indicated via \code{shared_link}.
#' @param cache_path If wanted, this allows to specify a path where the data
#' should be cached. if NULL (default), the default setting of \code{BiocFileCache()}
#' will be used.
#' @param check_for_update Indicate whether the data should be downloaded again.
#'
#' @return The final command here is \code{readRDS()}, so you should be assigning
#' the result of this function to an object name in your environment.
#'
#' @author Friederike Dündar
#'
#' @examples
#' direct_link <- "https://wcm.box.com/shared/static/4887keeu77eskmdotve1sm7n6u7udckg.rds"
#' testingsce <- load_RDSdata_from_Box(shared_link = direct_link, data_name = "testRDS")
#'
#' @export
#'
load_RDSdata_from_Box <- function(shared_link = "https://wcm.box.com/shared/static/mmf98o464n94fdh78x14sbwo7ts56ukq.rda",
  data_name = NULL, cache_path = NULL, check_for_update = FALSE){

  suppressPackageStartupMessages({ library(BiocFileCache) })
  load_new <- FALSE

  if(is.null(data_name)){
    data_name <- shared_link
  }

  ## load the cache
  if(is.null(cache_path)){
    bfc <- BiocFileCache(ask=F)
  }else{
    bfc <- BiocFileCache(cache_path, ask=FALSE)
  }

  ## check if url is being tracked
  bfc_res <- bfcquery(bfc, shared_link)

  if(bfccount(bfc_res) == 0L){
    load_new <- TRUE
  }else{
    ## if it is in cache, get path to load
    rid <- bfc_res[bfc_res$fpath == shared_link & bfc_res$rname == data_name,]$rid

    ## double-check we're looking at what we want
    if(length(rid) > 0){
      if(length(rid) > 1){
        message(paste0("Multiple entries found for", data_name, shared_link,". We'll use only the first one."))
        rid <- rid[1]
      }

      ans <- bfcrpath(bfc, rids = rid)

      if(check_for_update){
        ## check to see if the resource needs to be updated
        check <- bfcneedsupdate(bfc, rids = rid)
        if( is.na(check)) upd <- TRUE ## 'check' can be NA if it cannot be determined, choose how to handle
        if(upd) ans <- bfcdownload(bfc, rid = rid)
      }
    }else{
      load_new <- TRUE
    }
  }

  ## if it wasn't in the cache, add
  if(load_new){
    ans <- bfcadd(bfc, rname = data_name, fpath = shared_link)
  }

  message(paste("Cached data here:", ans))
  readRDS(ans)
}

#' Check for the presence of named columns
#'
#' @param which_names Column names to check, e.g. c("peaks", "genes")
#' @param input The names() of this will be checked.
#' @param input_name Name of the input, e.g. "start_params" just to make it
#' identifiable if an error message is returned.
#' @param function_name Name of function for which this test is carried out.
#' Again, just to make the returning error message a bit more readable. This
#' helps with debugging if you're wrapping many functions within each other.
#'
#' @return Returns an error and stop signal if entries of \code{which_names}
#'  are missing in the \code{input}.
#' @examples
#' \dontrun{
#' check_columns( c("cells", "sample", "condition"),
#'                long_df, "long_df", "plot_profile")
#' }
#'
#' @export
#'
check_columns <- function(which_names, input, input_name, function_name){

  check <- which_names %in% names(input)

  if( ! all( check )){
    stop(paste0("The input (", input_name , ") supplied to ", function_name,
      " is missing the following named columns: ",
      paste(which_names[!check], collapse = ", ") ))
  }
}
