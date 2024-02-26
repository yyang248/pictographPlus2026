simulation <- function(mcf=0, icn=2, minor_cn=1, depth=30, num_SNV=30, seed=NULL, normal_prop=0) {
  # mcf=0.5
  # icn=4
  # minor_cn=1
  # depth=100
  # num_SNV=100
  # normal_prop=0 # proportion of SNV is actually from CN-neutral region
  # seed=123

  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  tcn = 2 * (1-mcf) + icn * mcf
  
  # if a proportion of SNVs on CNA is actually on CN-neutral region
  num_neutral = round(num_SNV * normal_prop)
  SNV_depth_neutral = rpois(n = num_neutral, lambda = depth)
  SNV_alt_neutral = numeric(num_neutral)
  for (i in seq_len(num_neutral)) {
    SNV_alt_neutral[i] <- rbinom(n=1, size=SNV_depth_neutral[i], prob=0.5)
  }
  
  vaf_neutral = SNV_alt_neutral / SNV_depth_neutral
  
  # Actual tumor proportion
  num_SNV = num_SNV - num_neutral
  
  # generate depth for each SNV
  depth_total = rpois(n = num_SNV, lambda = depth * tcn / 2)
  
  # assign each SNV to one copy
  SNV_assignment = sample(c(1, 2), size = num_SNV, replace = TRUE) # which segment
  
  # generate depth for germline
  SNV_depth_germline = numeric(num_SNV)
  for (i in seq_len(num_SNV)) {
    SNV_depth_germline[i] <- rbinom(n=1, size=depth_total[i], prob=2 * (1-mcf)/tcn)
  }
  
  SNV_alt_germline = numeric(num_SNV)
  for (i in seq_len(num_SNV)) {
    SNV_alt_germline[i] <- rbinom(n=1, size=SNV_depth_germline[i], prob=0.5)
  }
  
  vaf_germline <- SNV_alt_germline/SNV_depth_germline
  
  # generate depth for tumor
  SNV_depth_tumor = depth_total - SNV_depth_germline 
  SNV_alt_tumor = numeric(num_SNV)
  for (i in seq_len(num_SNV)) {
    if (icn == 0) {
      SNV_alt_tumor[i] <- 0
    } else {
      tumor_vaf = minor_cn / icn
      if (SNV_assignment[i] == 1) {
        tumor_vaf = 1 - (minor_cn / icn)
      } 
      SNV_alt_tumor[i] <- rbinom(n=1, size=SNV_depth_tumor[i], prob=tumor_vaf)
    }
  }
  
  # vaf_tumor <- ifelse(SNV_depth_tumor==0, 0, SNV_alt_tumor/SNV_depth_tumor)
  vaf_tumor = c((SNV_alt_germline + SNV_alt_tumor) / depth_total, vaf_neutral)
  # return(vaf_neutral)
  title = paste("mcf: ", mcf, ", num_SNV: ", num_SNV, ", depth: ", depth, ", \nICN: ", icn, ", minor CN: ", minor_cn, 
                ", true_prop: ", 1-normal_prop, ", unimodal: ", is.unimodal(vaf_tumor), sep="")
  plot(density(vaf_tumor), xlim=c(0,1), main = title, cex.main=0.9)

  # plot(density(vaf_tumor), xlim=c(0,1), main = "VAF_tumor")
  # qqnorm(vaf_tumor, main = "tumor")
  # qqline(vaf_tumor, col="grey")
  
  is.unimodal(vaf_tumor)
  # return(um)
}

simulation_normal <- function(depth=30, num_SNV=30) {
  
  # generate depth for each SNV
  SNV_depth = rpois(n = num_SNV, lambda = depth)

  # generate alt read counts for each SNV
  SNV_alt = numeric(num_SNV)
  for (i in seq_len(num_SNV)) {
    SNV_alt[i] <- rbinom(n=1, size=SNV_depth[i], prob=0.5)
  }
  
  # calculate VAF for all SNVs
  vaf <- SNV_alt/SNV_depth
  
  # plot the SNVs
  title = paste("num_SNV: ", num_SNV, ", depth: ", depth, ", unimodal: ", is.unimodal(vaf), sep="")
  plot(density(vaf), xlim=c(0,1), main = title)
}

