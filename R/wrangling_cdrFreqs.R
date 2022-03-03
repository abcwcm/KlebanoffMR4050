
#' Add the occurrences of specific TRA/TRB chains to the colData of the SCE
#'
#' @param sce.object SCE object
add_chain_counts <- function(sce.object){

    if("cdr3s_aa" %in% names(colData(sce.object))){
        ## counting TRA, TRB
        sce.object$numTRA <-  stringr::str_count(sce.object$cdr3s_aa,"TRA")
        sce.object$numTRB <-  stringr::str_count(sce.object$cdr3s_aa,"TRB")

        ## multi-TRB
        sce.object$multiTRB <-  sce.object$numTRB  >= 2
        # sce.object$multiTRB  <- as.character(sce.object$multiTRB)
        #  colData(sce.object)[is.na( sce.object$multiTRB ),]$multiTRB <- as.character("NA")

        ## multi-TRA
        sce.object$multiTRA <-  sce.object$numTRA  >= 2
        #sce.object$multiTRA  <- as.character(sce.object$multiTRA)
        #colData(sce.object)[is.na( sce.object$multiTRA ),]$multiTRA <- as.character("NA")

        return(sce.object)
    }else{
        stop("The sce.object must have an entry in colData named 'cdr3s_aa'.")
    }

}

#' Count the number of times individual cdr3s_aa entires occur
#'
#' @param sce.object SCE object
#' @param count_what Whether IDs or amino acid sequences should be counted.
#' Default: "cdr3s_aa"
#' @param mode either one of "all" or "per.sample" -- determines how
#' the number of individual cdr3s_aa entries is counted
#' @return vector of numbers
#' @examples \dontrun{
#' sce.filt$freq_per_Sample <- add_frequencies(sce.filt, mode = "per.sample")
#' sce.filt$freq_across_all <- add_frequencies(sce.filt, mode = "all")
#' }
add_frequencies <- function(sce.object, count_what = "cdr3s_aa", mode = c("all","per.sample")){
    if(all( c(count_what, "Sample") %in% names(colData(sce.object)) )){
        if(mode == "per.sample"){
            cd <- colData(sce.object)[, c("Sample",count_what)] %>%
                as.data.frame %>% as.data.table(., keep.rownames=TRUE)
            cd <- cd[, .N,by=c(count_what,"Sample")] %>%
                merge(., cd, by = c(count_what,"Sample"), all = TRUE) %>% as.data.frame
            rownames(cd) <- cd$rn
        }else if(mode == "all"){
            cd <- colData(sce.object)[, count_what, drop=FALSE] %>%
                as.data.frame %>% as.data.table(., keep.rownames=TRUE)
            cd <- cd[, .N,by=count_what] %>%
                merge(., cd, by = count_what, all = TRUE) %>% as.data.frame
            rownames(cd) <- cd$rn
        }else{
            stop("Select a mode: 'all' or 'per.sample'.")
        }
        return(cd[colnames(sce.object),]$N)
    }else{
        stop(paste("The sce.object must have colData entries named 'Sample' and", count_what,"."))
    }
}