#' Function for printing pretty table
#'
#' @param data.frame data.frame to be printed
#' @param show_rownames TRUE or FALSE (default)
#'
#' @import kableExtra
#' @export
#'
my_table <- function(df, show_rownames = FALSE){

    if( !all(class(df) != "data.frame")){
        df <- as.data.frame(df)
    }

    kable(df, row.names=FALSE, padding = 0, longtable=TRUE) %>%
        kable_styling(bootstrap_options = c("striped", "hover", "condensed"))
}