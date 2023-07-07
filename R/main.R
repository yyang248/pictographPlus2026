#' main steps to run pictograph2
#' @export
runPictograph2 <- function(){
  # read in files
  input_data <- importFiles('./inst/extdata/sim_v2_snv.csv', './inst/extdata/sim_v2_cn.csv')

  # separate by sample presence
  sep_list <- separateBySamplePresence(input_data)

  # initial integer CNA estimation by setting CNA to closest integer depending on the mean across all samples
  cna_est <- estimateCNA(input_data$tcn)

  # estimation of copy number cellular fraction
  cncf_est <- estimateCNCF(input_data$tcn, cna_est)

  # estimation of mutation cellular fraction for each mutation


  # return(input_data)
}
