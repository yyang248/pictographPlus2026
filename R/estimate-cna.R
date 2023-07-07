#' cluster CNAs and resolve CNA state
#' @export
estimateCNA <- function(tcn) {

  # initial integer CNA estimation by setting CNA to closest integer depending on the mean across all samples
  cna_est <- tcn
  for (i in seq_len(nrow(cna_est))) {
    # if mean CNA across all samples >= 2
    if (mean(cna_est[i,]) >= 2) {
      # set CNA to ceiling of max of CNAs
      cna_est[i,] = ceiling(max(cna_est[i,]))
    } else {
      # otherwise, set CNA to floor of min of CNA
      cna_est[i,] = floor(min(cna_est[i,]))
    }
  }
  warning("need to consider if tcn above/below 2 coexist; MCMC to sample cna")
  return(cna_est)
}

#' estimation of copy number cellular fraction
#' @export
estimateCNCF <- function(tcn, cna_est) {
  cncf_est <- (tcn - 2) / (cna_est - 2)
  warning("need to consider if tcn above/below 2 coexist; consider MCMC for sampling")
  return(cncf_est)
}
