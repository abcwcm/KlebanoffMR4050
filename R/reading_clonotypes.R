#' Reading in clonotype info per barcode
#'
#' @details This function expects the files "filtered_contig_annotations.csv" and
#' "clonotypes.csv" to be found in the "outs" folder after following \code{cellRanger_result_dir}.
#' These files are used to match the cell barcodes with unique clonotype IDs.
#'
#' @param cellRanger_result_dir e.g. "/scratchLocal/frd2007/2018-11_Smita_Tcells/data/wildtype/".
#'  This should be the path to the un-tar'ed directory from cellRanger.
#'
#'  @export
#'
reading_in_clonotypes <- function(cellRanger_result_dir){

  ## reading barcodes (including empty ones)
  bcs <- fread(paste0(cellRanger_result_dir,"outs/filtered_contig_annotations.csv"), sep = ",")

  #> head(barcode_csv)
  #barcode is_cell                   contig_id high_confidence length
  #1 AAACCTGAGTGACATA-1    True AAACCTGAGTGACATA-1_contig_1            True    568
  #2 AAACCTGAGTGACATA-1    True AAACCTGAGTGACATA-1_contig_3            True    567
  #3 AAACCTGAGTGACATA-1    True AAACCTGAGTGACATA-1_contig_4            True    712
  #4 AAACCTGCACCGCTAG-1    True AAACCTGCACCGCTAG-1_contig_1            True    564
  #5 AAACCTGCACCGCTAG-1    True AAACCTGCACCGCTAG-1_contig_2            True    584
  #6 AAACCTGCACTTCTGC-1    True AAACCTGCACTTCTGC-1_contig_1            True    569
  #chain   v_gene d_gene  j_gene c_gene full_length productive
  #1   TRA   TRAV19   None  TRAJ50   TRAC        True       True
  #2   TRA TRAV13-1   None  TRAJ45   TRAC        True       None
  #3   TRB   TRBV19  TRBD2 TRBJ2-7  TRBC2        True       True
  #4   TRB    TRBV9  TRBD2 TRBJ2-6  TRBC2        True       True
  #5   TRA TRAV12-2   None  TRAJ44   TRAC        True       True
  #6   TRA   TRAV25   None  TRAJ54   TRAC        True       True
  #cdr3                                                   cdr3_nt
  #1     CALSSMKTSYDKVIF             TGTGCTCTGAGTAGCATGAAAACCTCCTACGACAAGGTGATATTT
  #2                None                                                      None
  #3      CATSLALAAYEQYF                TGTGCCACGAGTCTGGCGCTAGCGGCGTACGAGCAGTACTTC
  #4 CASSGLAGGPVSGANVLTF TGTGCCAGCAGCGGACTAGCGGGAGGGCCCGTTTCTGGGGCCAACGTCCTGACTTTC
  #5       CAGNTGTASKLTF                   TGTGCCGGAAATACCGGCACTGCCAGTAAACTCACCTTT
  #6        CAGLQGAQKLVF                      TGTGCAGGCCTTCAGGGAGCCCAGAAGCTGGTATTT
  #reads umis raw_clonotype_id         raw_consensus_id
  #1  2511    7     clonotype146 clonotype146_consensus_1
  #2   679    2     clonotype146                     None
  #3   850    4     clonotype146 clonotype146_consensus_2
  #4 10237   36       clonotype1   clonotype1_consensus_1
  #5  3470   10       clonotype1   clonotype1_consensus_2
  #6  1780    7     clonotype147 clonotype147_consensus_1

  ## extract those with a clonotype
  bcs.sub <- bcs[ raw_clonotype_id != "None", c("raw_clonotype_id","barcode"), with=FALSE] %>% unique

  ## reading clonotypes
  clotp <- fread(paste0(cellRanger_result_dir, "outs/clonotypes.csv"), sep = ",")
  setnames(clotp, "clonotype_id", "raw_clonotype_id")

  ## join the informationbiocon
  jnd <- merge(clotp, bcs.sub, by = "raw_clonotype_id", all = TRUE)

  ## sanity check
  tst <- jnd[raw_clonotype_id == "clonotype1", .N] ==  unique(jnd[raw_clonotype_id == "clonotype1"]$frequency)
  if(isFALSE(tst)){
    warning("The number of rows for clonotype1 does not match the frequency taken from the clonotypes.csv file.")
  }

  return(jnd)
}
