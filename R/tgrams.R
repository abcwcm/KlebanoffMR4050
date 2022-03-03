#' findMarkers wrapper
#'
#' @description Wrapper function for running \code{scran::findMarkers}.
#' Will only return the results for the MUT condition, i.e. logFC values
#' represent MUT/WT.
#'
#' @param sce.object SingleCellExperiment object
#' @param group_identifier define the variable that will be used to determine which
#' samples belong to the mutant or the wildtype condition, default: "Sample"
#' @param direction parameter for \code{scran::findMarkers}; default: "up"
#' @param logfc parameter for code{scran::findMarkers}
#' @param FDR parameter for \code{extract_markers}, default: 0.01
#' @param rank parameter for \code{extract_markers}, default: Inf
#' @param comp_name parameter for \code{extract_markers}; this should be
#' the clonotype for which the comparison is being done
#'
#' #@return list with both the findMarkers results (which will consist of lists of 2)
#' and the extract_markers results, which will be a data.table with gene names,
#' but without the statistical metrics
#'
#' @import data.table
#' @import magrittr
#'
#' @export
run_DE <- function(sce.object, group_identifier = "Sample", direction = "up",
  logfc = 0.2, FDR = 0.01, rank = Inf, comp_name = NULL){

  if(!group_identifier %in% names(colData(sce.object))){
    stop(paste(group_identifier, "is not part of the colData of the sce.object."))
  }


  if( length(unique(colData(sce.object)[, group_identifier, drop = TRUE])) < 2 ){
    message(paste("run_DE: There is only one condition for", group_identifier, ", which makes the comparison impossible"))
  }else{
      if(all(table(colData(sce.object)[,group_identifier])>1)){
          fM <- scran::findMarkers(logcounts(sce.object),
              groups=colData(sce.object)[,group_identifier],
              direction=direction,
              lfc= logfc)

          deres <- extract_markers(sce.object,
              fM,
              FDR_thresh = FDR,
              rank_thresh = rank,
              comp_name = comp_name)
          fM$WT <- NULL ## removing the duplicated information from the WT results
          return(list(findMarkers_results = fM, marker_IDs = deres))
      }
      else{message(paste("We need more than 1 member per", group_identifier))}
  }
}



#' Prepare data.table for plotting
#'
#' @description This function needs the clonotypes that are going to be plotted
#' and for which findMarkers results have been obtained via \code{run_DE}. It
#' will return a data.table with expression values (logcounts) for the gene of
#' interest (\code{goi}) as well as for CD4 and CD8A.
#'
#' @param sce.object SingleCellExperiment object where the colData should contain "Sample",
#' "cdr3s_aa", and "id" (= unique clonotype ID)
#' @param which_clonotypes vector of clonotype IDs that are going to be shown, e.g.
#' "C19"
#' @param DEresult_list list of findMarkers and extract_markers results for all
#' clonotypes present in \code{sce.object}
#' @param ignore_zeros set to TRUE to remove zeros from the expression data (Default: FALSE)
#' @param goi gene of interest for which the asterisks will be defined, e.g. "IFNG"
#' @param additional_colData additional information from the colData of the SCE to be added
#' to the returned data.table, e.g. "antigen".
#'
#' @return data.table
#'
#' @seealso \code{run_DE}, \code{plot_tgram}
#'
#' @examples \dontrun{
#' > genes_of_interest <- c("IFNG","CD4","CD8A","CD8B","CD3E","TNF","IL2","IL23A")
#' > tt <- prep_data4tgram(sce.shared,
#'     which_clonotypes = unique(sce.shared$id),
#'     DEresult_list = delist,
#'     goi = genes_of_interest)
#' > tt$Sample <- factor(tt$Sample, levels = c("wt","mut"), ordered = TRUE)
#' > plot_tgram(tt[gene_name == "IFNG"], dot_shape = 20, col_wt = "limegreen",
#'     col_mut = "dodgerblue1", label_font_size = 10, show_legend = FALSE)
#' }
#'
#' @export
#'
prep_data4tgram <- function(sce.object, which_clonotypes = NULL, DEresult_list, ignore_zeros = FALSE, goi,
  additional_colData = NULL){

  check_columns(c("id", "cdr3s_aa", "Sample"),
    colData(sce.object), "sce.object","prep_data4tgram")

  ## expects that sample names are named MUT.patient | WT.patient
  if(!("antigen" %in% names(colData(sce.object)))){
    sce.object$antigen <- factor(gsub("\\..*","", sce.object$Sample),
      levels = c("WT","MUT"), ordered = TRUE)
    if(any(is.na(sce.object$antigen))){stop("Make sure that your Sample entries are labelled MUT.patient and WT.patient.")}
  }

  ## GET EXPRESSION VALUES ===============================================
  if(!is.null(which_clonotypes)){
    which_cells <- rownames(subset(colData(sce.object), id %in% unique(which_clonotypes)))
  }else{
    which_cells <- colnames(sce.object)
  }

  metdat <- c("Sample", "cdr3s_aa","id", "antigen")
  if(!is.null(additional_colData)){
    if(all(additional_colData %in% names(colData(sce.object)) )){
      metdat <- unique(c(metdat, additional_colData))
    }else{
      found <- additional_colData %in% names(colData(sce.object))
      warning(paste( paste(additional_colData[!found], collapse = ","), "not found in names(colData(sce.object))." ))
      additional_colData <- additional_colData[found]
      if(length(additional_colData) > 0){ metdat <- unique(c(metdat, additional_colData))}
    }
  }

  ldt <- make_long_dt(sce.object[, which_cells],
    exprs_values = "logcounts",
    features = goi, include_metaData = metdat)
  setnames(ldt, "feature_name", "gene_name")

  ## GET FDR VALUES ========================================================
  ### extract the findMarker results
  de_gns <- lapply(DEresult_list, function(x) x[[2]]) %>% rbindlist(., idcol = "id")

  ### extract the FDR values for the gene of interest
  if(any(goi %in% de_gns$gene_symbol)){ ## only if at least one of the genes is actually part of de_gns
    fdr_vals <- lapply(goi, function(gn){
        if( dim(de_gns[gene_symbol == gn])[1] >0 ){
          lapply(de_gns[gene_symbol == gn]$id,
          function(x){ # for every clonotype, for which the gene showed up as DE, extract the corresponding row
            ## since I tend to only test for upregulation, I need to find the FDR that's smallest
            gns.MUT <- DEresult_list[[x]][[1]]$MUT
            goival.MUT <- gns.MUT[grepl(paste0("^",gn,"$"), rownames(gns.MUT) ), ]
            gns.WT <- DEresult_list[[x]][[1]]$WT
            goival.WT <- gns.WT[grepl(paste0("^",gn,"$"), rownames(gns.WT) ), ]

            out <- data.table(
              id = x,
              FDR = min(goival.WT$FDR, goival.MUT$FDR, na.rm = TRUE),
              gene_name = gn)
          }) %>% rbindlist
      }
    }) %>% rbindlist
  }else{
    fdr_vals <- NULL
  }

  ## ADD INFO ABOUT STATS FOR GOI=========================================
  ### FDR
  if(ignore_zeros ==  TRUE){
    ldt <- ldt[ logcounts > 0 ]
  }

  if(!is.null(fdr_vals)){
    ldt <- merge(ldt, fdr_vals, on = c("gene_name","id"), all.x = TRUE)
    ldt[, star := ifelse(is.na(FDR), "",
      ifelse(FDR < 0.0001, "***",
        ifelse(FDR < 0.001, "**",
          ifelse(FDR <= 0.05, "*", ""))))]
  }else{
    ldt[, FDR := NA]
    ldt[, star := ""]
  }

  ldt[]
  #cdr3s_aa gene_name Sample    cell logcounts FDR star
  #1:  TRA:CAAASGGYQKVTF;TRB:CAWRTSGTYEQYF      CD8A    mut mut2239 3.7377803  NA
  #2:  TRA:CAAASGGYQKVTF;TRB:CAWRTSGTYEQYF      CD8A     wt  wt3040 4.5200104  NA
  #3:  TRA:CAAASGGYQKVTF;TRB:CAWRTSGTYEQYF      CD8A     wt  wt4364 4.8598710  NA
  #4: TRA:CAAFSGSARQLTF;TRB:CASSLNSGGYGYTF       CD4    mut mut1576 1.6009369  NA
  #5: TRA:CAAFSGSARQLTF;TRB:CASSLNSGGYGYTF       CD4    mut mut1699 1.9583623  NA
  #---
  #  1628:                   TRB:CSADNSTTNEKLFF      CD8A     wt  wt1606 2.1827861  NA
  #1629:                   TRB:CSADNSTTNEKLFF      IFNG    mut  mut955 0.5473646  NA
  #1630:                   TRB:CSADNSTTNEKLFF      IFNG    mut mut1762 0.8299638  NA
  #1631:                   TRB:CSADNSTTNEKLFF      IFNG     wt   wt430 2.1774460  NA
  #1632:                   TRB:CSADNSTTNEKLFF      IFNG     wt  wt1606 3.8054042  NA

  return(ldt)

}


#' Plot THE PLOT
#'
#' @description This function generates a plot that Chris Klebanoff calls
#' "Manhattan Plot". Basically, it will show the expression levels of the
#' gene(s) across different sample types and clonotypes. In addition, it
#' will add asterisks to those genes that were determined to be statistically
#' significantly different as provided in the results of \code{run_DE}.
#'
#' @param in_dt result of \code{prep_data4tgram}
#' @param boxplot_stroke thickness of the boxplot outline, default: 0.2
#'
#' @return ggplot2 object
#'
#' @examples \dontrun{
#' delist <- lapply(ct$cdr3s_aa, function(x) run_DE(sce.filt2[, sce.filt2$cdr3s_aa == x],
#'                                                  group_identifier = "Sample",
#'                                                  direction = "up",
#'                                                  logfc = 0.2,
#'                                                  FDR = 0.05, rank = Inf,
#'                                                  comp_name = paste0(x, "_")))
#' names(delist) <- ct$cdr3s_aa
#' tt <- prep_data4tgram(sce.filt2, which_clonotypes = ct$cdr3s_aa,
#'                       DEresult_list = delist, goi = "IFNG")
#' plot_tgram(tt, dot_shape = 20, col_wt = "limegreen",col_mut = "dodgerblue1", label_font_size = 10)
#' }
#'
#' @seealso \code{prep_data4tgram}, \code{run_DE}
#' @export
#'
plot_tgram <- function(in_dt = demks.dt,
  label_font_size = 12,
  gene_label_size = 16,
  dot_size = 1,
  dot_shape = 21,
  boxplot_stroke = .2,
  col_mut = "gray60",
  col_wt = "white",
  show_legend = TRUE){

  ## ERROR CHECKS ==========================================================
  check_columns(c("gene_name","id","logcounts","Sample","star", "antigen"),
    input = in_dt, input_name = "in_dt","plot_tgram")

  #goi <- unique(in_dt$gene_name)[!unique(in_dt$gene_name) %in% c("CD4","CD8A")]

  #if(length(goi) > 1){
  #  warning(paste("There were more than one gene different from CD4 and CD8A in the in_dt:", goi, "We will only use the first one." ))
  #  goi <- goi[1]
  #}
  #if(length(goi) == 0){
  #  warning("There was no other gene than CD4 or CD8A in the in_dt.")
  #}


  ## defining plot limits
  ymax <- max(in_dt$logcounts) + 1
  ystar <- max(in_dt$logcounts) + 0.3
  #med_goi <- median(in_dt[gene_name == goi]$logcounts)
  med_goi <- median(in_dt$logcounts)
  #in_dt$gene_name <- factor(in_dt$gene_name, levels = c(goi, "CD8A","CD4"), ordered = TRUE)
  in_dt$id <- factor(in_dt$id, levels = gtools::mixedsort(unique(in_dt$id)), ordered = TRUE)
  in_dt$antigen <- factor(in_dt$antigen, levels = c("WT", "MUT"), ordered = TRUE)

  ### GENERATE the PLOT ===============================================
  if(dot_shape %in% c(1:20)){
    filling <- FALSE
    pl_out <- ggplot(in_dt, aes(x = id,
      y = logcounts, color = antigen))
  }else{
    filling <- TRUE
    pl_out <- ggplot(in_dt, aes(x = id,
      y = logcounts, fill = antigen))
  }

  pl_out <- pl_out +
    ## defining size of the plot
    coord_cartesian(ylim = c(0, ymax)) +
    ## add points -------------
  ggbeeswarm::geom_quasirandom(groupOnX = TRUE, dodge.width=0.8,
    alpha = .6, size = dot_size, shape = dot_shape) +
    ## add boxplot -----------------
  geom_boxplot(alpha = 0, inherit.aes = FALSE,
    show.legend = FALSE, lwd = boxplot_stroke,
    aes(x = id, y = logcounts, fill = antigen))

  ## theme_bw -------------
  pl_out <- pl_out + theme_bw() +
    theme(strip.text.x = element_text(size = gene_label_size))

  ## modify legend and axis display ---------
  pl_out <- pl_out +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,
      vjust = 0.5, size = label_font_size))

  ## legend ----------------------
  if(show_legend == TRUE){
    pl_out <- pl_out +
      theme(legend.justification=c(0,1), legend.position=c(0,1),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"))
  }else{
    if(filling == TRUE){
      pl_out <- pl_out + guides(fill = FALSE)
    }else{pl_out <- pl_out + guides(color = FALSE)}
  }

  ## axis labels --------------
  pl_out <- pl_out +
    xlab("") + ylab("logcounts")
  #scale_y_continuous(breaks = seq(0, ceiling(max(in_dt$logcounts)), 1.5))
  ## defining shapes ----------
  #scale_shape_manual(values = c(1,15)) +

  ## defining colors -----------
  if(filling == TRUE){
    pl_out <- pl_out + scale_fill_manual(values = c(col_wt,col_mut))
    if(show_legend == TRUE){
      pl_out <- pl_out + guides(fill=guide_legend(ncol=2,
        override.aes = list(alpha = 1, size = 4)))
    }
  }else{
    pl_out <- pl_out +  scale_color_manual(values = c(col_wt,col_mut))
    if(show_legend == TRUE){
      pl_out <- pl_out + guides(color=guide_legend(ncol=2,
        override.aes = list(alpha = 1, size = 4)))
    }
  }

  ## grid---------------------------
  pl_out <- pl_out +
    #facet_wrap(~gene_name, scales = "free_y", ncol = 1) +
    facet_grid(~gene_name) +
    ## "threshold" line
    #geom_hline(data = in_dt[gene_name == goi], aes(yintercept = med_goi), linetype = "dashed", color = "gray") +
    geom_hline(data = in_dt, aes(yintercept = med_goi), linetype = "dashed", color = "gray") +
    ## add asterisks
    geom_text(data = unique(in_dt[, c("id", "star","Sample","gene_name","antigen"), with =FALSE]),
      inherit.aes = FALSE,
      aes(x = id, y = ystar, label=star), colour="black", size=8)

  return(pl_out)
}

#
# * xaxis labels must be enabled to be turned off
# * legend must be enabled to be turned off, too
