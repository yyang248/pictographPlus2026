simulation <- function(mcf=0, icn=2, minor_cn=1, depth=30, num_SNP=30, seed=NULL, normal_prop=0) {
  # mcf=0.8
  # icn=1
  # minor_cn=0
  # depth=100
  # num_SNP=1500
  # normal_prop=0 # proportion of SNP is actually from CN-neutral region
  # seed=123

  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  tcn = 2 * (1-mcf) + icn * mcf
  
  # if a proportion of SNPs on CNA is actually on CN-neutral region
  num_neutral = round(num_SNP * normal_prop)
  SNP_depth_neutral = rpois(n = num_neutral, lambda = depth)
  SNP_alt_neutral = numeric(num_neutral)
  for (i in seq_len(num_neutral)) {
    SNP_alt_neutral[i] <- rbinom(n=1, size=SNP_depth_neutral[i], prob=0.5)
  }
  
  vaf_neutral = SNP_alt_neutral / SNP_depth_neutral
  
  # Actual tumor proportion
  num_SNP = num_SNP - num_neutral
  
  # generate depth for each SNP
  depth_total = rpois(n = num_SNP, lambda = depth * tcn / 2)
  
  # assign each SNP to one copy
  SNP_assignment = sample(c(1, 2), size = num_SNP, replace = TRUE) # which segment
  
  # generate depth for germline
  SNP_depth_germline = numeric(num_SNP)
  for (i in seq_len(num_SNP)) {
    SNP_depth_germline[i] <- rbinom(n=1, size=depth_total[i], prob=2 * (1-mcf)/tcn)
  }
  
  SNP_alt_germline = numeric(num_SNP)
  for (i in seq_len(num_SNP)) {
    SNP_alt_germline[i] <- rbinom(n=1, size=SNP_depth_germline[i], prob=0.5)
  }
  
  vaf_germline <- SNP_alt_germline/SNP_depth_germline
  
  # generate depth for tumor
  SNP_depth_tumor = depth_total - SNP_depth_germline 
  SNP_alt_tumor = numeric(num_SNP)
  for (i in seq_len(num_SNP)) {
    if (icn == 0) {
      SNP_alt_tumor[i] <- 0
    } else {
      tumor_vaf = minor_cn / icn
      if (SNP_assignment[i] == 1) {
        tumor_vaf = 1 - (minor_cn / icn)
      } 
      SNP_alt_tumor[i] <- rbinom(n=1, size=SNP_depth_tumor[i], prob=tumor_vaf)
    }
  }
  
  # vaf_tumor <- ifelse(SNP_depth_tumor==0, 0, SNP_alt_tumor/SNP_depth_tumor)
  vaf_tumor = c((SNP_alt_germline + SNP_alt_tumor) / depth_total, vaf_neutral)
  # return(vaf_neutral)
  title = paste("mcf: ", mcf, ", num_SNV: ", num_SNP, ", depth: ", depth, ", \nICN: ", icn, ", minor CN: ", minor_cn, 
                ", true_prop: ", 1-normal_prop, ", unimodal: ", is.unimodal(vaf_tumor), sep="")
  plot(density(vaf_tumor), xlim=c(0,1), main = title, cex.main=0.9)

  # plot(density(vaf_tumor), xlim=c(0,1), main = "VAF_tumor")
  # qqnorm(vaf_tumor, main = "tumor")
  # qqline(vaf_tumor, col="grey")
  
  is.unimodal(vaf_tumor)
  # return(um)
}
