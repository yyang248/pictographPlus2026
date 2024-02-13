#' read input data file and store in PICTOGRAPH input format
#' @export
#' @param mutation_file mutation data file that contains columns "sample", "mutation", "chrom", "start", "end", total_reads", and "alt_reads";
#' @param copy_number_file copy number file that contains columns "sample", "chrom", "start", "end", "tcn"
importFiles <- function(mutation_file, copy_number_file, outputDir, SNP_file=NULL, cytoband_file=NULL, alt_reads_thresh = 0, vaf_thresh = 0, cnv_max_dist=2000, cnv_max_percent=0.30, tcn_normal_range=c(1.8, 2.2), smooth_cnv=F, autosome=T) {
  
  # keep mutations if alt_reads >= alt_reads_thresh and vaf >= vaf_thresh
  mutation_data = importMutationFile(mutation_file, alt_reads_thresh, vaf_thresh)
  
  copy_number_data = importCopyNumberFile(copy_number_file, germline_SNP_file, outputDir, SNP_file, cnv_max_dist, cnv_max_percent, tcn_normal_range, smooth_cnv, autosome)

  mutation_data$tcn = copy_number_data$tcn
  
  # copy_number_info = importCopyNumberInfo(copy_number_file)

  # mutation_data$cytoband = copy_number_info$cytoband
  # mutation_data$drivers = copy_number_info$drivers
  # mutation_data$genes = copy_number_info$genes
  # 
  mutation_data$overlap = resolveOverlap(mutation_data)
  warning("resolveOverlap: need to check whether a mutation overlaps with two CNA segs")

  return(mutation_data)
}

#' check whether a mutation overlaps a CNA region
resolveOverlap <- function(mutation_data) {
  mut_count = nrow(mutation_data$position)
  cna_count = nrow(mutation_data$tcn)
  output = matrix(0, nrow=mut_count, ncol=cna_count)
  rownames(output) = mutation_data$position$mutation
  colnames(output) = rownames(mutation_data$tcn)
  # print(output)
  cnas = rownames(mutation_data$tcn)
  mutations = mutation_data$position
  for (i in seq_len(mut_count)) {
    for (j in seq_len(cna_count)) {
      chrSplit = strsplit(cnas[j], split='-')[[1]]
      chr = chrSplit[1]
      start_pos = strtoi(chrSplit[2])
      end_pos = strtoi(chrSplit[3])

      if (chr==mutations[i,2]) {
        if (strtoi(mutations[i,3]) >= start_pos && strtoi(mutations[i,4]) <= end_pos) {
          output[i,j] = 1
        }
      }
    }
  }
  return(output)
}

#' merge segments with similar coordinates
#' the maximum allowed distance between the start or end position of two segments is the max(cnv_max_dist, min(length of two segments)*cnv_max_percent)
smoothCNV <- function(data, cnv_max_dist=2000, cnv_max_percent=0.30) {
  # smooth_data <- as_tibble(data.frame(matrix(nrow=0, ncol=ncol(data))))
  # colnames(smooth_data) <- colnames(data)
  
  # samples <- unique(data$sample)
  
  message("improvements needed for smoothCNV under preprocessing.R; currently detecting two segments as the same if overlap and star/end position within certain distance; merge the two segments as the outcome")
  data <- data %>% 
    mutate(numeric_chrom = as.numeric(sub("chr", "", chrom))) %>% 
    arrange(numeric_chrom, start) %>%
    select(-numeric_chrom)
  
  indexList = seq_len(nrow(data))
  
  for (i in seq_len(nrow(data))) {
    if (i %in% indexList) {
      # print(data[i,])
      max_dis = max(cnv_max_dist, cnv_max_percent*(data[i,]$end-data[i,]$start))
      index_selected = which((data$chrom==data[i,]$chrom)&(data$start<=data[i,]$end)&(data[i,]$start<=data$end)&(abs(data$start-data[i,]$start)<max_dis)&(abs(data$end-data[i,]$end)<max_dis))
      # index_selected = which((data$chrom==data[i,]$chrom)&(abs(data$start-data[i,]$start)<max_dis)&(abs(data$end-data[i,]$end)<max_dis))
      if (length(index_selected) > 1) {
        # list of cnv segments able to merge
        cnv_selected = data[index_selected,]
        # start = min(cnv_selected$start)
        # end = max(cnv_selected$end)
        data[index_selected,]$start = min(cnv_selected$start)
        data[index_selected,]$end = max(cnv_selected$end)
        data[index_selected,]$drivers = paste(unique(unlist(strsplit(cnv_selected$drivers, ";"))), collapse=";")
        data[index_selected,]$genes = paste(unique(unlist(strsplit(cnv_selected$genes, ";"))), collapse=";")
      }
      
      indexList <- indexList[!(indexList %in% index_selected)]
    }
  }
  
  return(data)
}

#' import copy number file
#' @param cnv_max_dist: maximum of distance allowed between two segments to assign as the same one
#' @param cnv_max_percent: maximum percentage of distance allowed between two segments to assign as the same one
#' @param smooth_cnv: process input CNV to merge  segments with similar distance
importCopyNumberFile <- function(copy_number_file, outputDir, SNP_file=NULL, cnv_max_dist=2000, cnv_max_percent=0.30, tcn_normal_range=c(1.8, 2.2), smooth_cnv=F, autosome=T) {
  # read copy number csv file
  data <- read_csv(copy_number_file, show_col_types = FALSE)

  # Smooth segments so segments with different start/end position are treated the same, 
  # conditioned on overlapping and distance between positions
  if (smooth_cnv) {
    data <- smoothCNV(data, cnv_max_dist, cnv_max_percent)
  }
  
  data$CNA = paste(data$chrom, data$start, data$end, sep = '-')
  
  # keep only autosome
  if (autosome) {
    data <- data[!(data$chrom %in% c('chrX', 'chrY', 'chrM')),]
  }
  
  # convert to cna by sample matrix 
  output_data <- list()
  output_data = as.matrix(data[c("sample", "CNA", "tcn")] %>% pivot_wider(names_from = sample, values_from = tcn, values_fill = 2))
  
  # fill missing CNA with 2
  message("FILL MISSING TCN WITH 2; NEED TO FIX FOR CHRX/Y FOR VALUES_FILL in importCopyNumberFile in preprocessing.R")
  rownames(output_data) <- output_data[,'CNA']
  output_data <- output_data[,-1, drop=FALSE]
  rowname = rownames(output_data)
  colname = colnames(output_data)
  output_data <- matrix(as.numeric(output_data), ncol = ncol(output_data))
  rownames(output_data) = rowname
  colnames(output_data) = colname
  
  # check lOH using HET SNPs; 
  # if both germline and tumor SNP counts provided, will check LOH
  # otherwise, keep CNA outside normal range only
  if (!is.null(SNP_file)) {
    to_keep_index = check_LOH(output_data, outputDir, SNP_file, tcn_normal_range)
  } else {
    to_keep_index = which(rowSums(output_data<tcn_normal_range[1]|output_data>tcn_normal_range[2])>0)
  }
  

  
  output_data = output_data[to_keep_index,,drop=FALSE]
  
  
  return(output_data)
}

#' @import LaplacesDemon
check_LOH <- function(output_data, outputDir, SNP_file, tcn_normal_range=c(1.8, 2.2), pval=0.05) {
  print(output_data)
  
  tcn_normal_range=c(1.8, 2.2)
  pval=0.05
  SNP_data <- read_csv(SNP_file)
  samples <- colnames(output_data)
  output <- rownames_to_column(as.data.frame(output_data), var="chrom") 
  output <- output %>% separate(col="chrom", into = c("chrom", "start", "end"), sep="-")
  output <- output %>% mutate(start = as.integer(start), end = as.integer(end))
  
  pdf(paste(outputDir, "test.pdf", sep=""), width = (length(samples)+1)*6, height = (length(samples)+1)*9)
  par(mfrow=c(length(samples)*3,(length(samples)+1)*2), cex.main = 1, cex.axis = 1, cex.lab = 1, cex = 1)
  
  rowLen =10
  rowLen = nrow(output_data)
  
  tcn_ref = output_data
  tcn_ref[,] <- 0
  tcn_alt = output_data
  tcn_alt[,] <- 0
  
  to_keep_index = c()
  
  set.seed(123)
  for (i in 1:rowLen) {
    # 326; 11
    # set.seed(123)

    i = 326
    # print(output[i,])
    SNP_temp <- SNP_data %>% filter(chroms==output[i,]$chrom & position>=output[i,]$start & position<=output[i,]$end)
    SNP_temp
    if (nrow(SNP_temp) > 8) {
      # for a given segment, if tcn in all samples within tcn_normal_range, assuming copy neutral and test of LOH
      # using Kolmogorov-Smirnov Test
      vaf_germline <- SNP_temp$germline_alt / (SNP_temp$germline_alt + SNP_temp$germline_ref)
      if (is.unimodal(vaf_germline)) {
        # vaf_germline <- SNP_temp$germline_alt / (SNP_temp$germline_alt + SNP_temp$germline_ref)
        title = paste(output[i,]$chrom, output[i,]$start, output[i,]$end, sep="_")
        plot(density(vaf_germline), xlim=c(0,1), main = title)
        qqnorm(vaf_germline, main = title)
        qqline(vaf_germline, col="grey")
        
        ks_list = c()
        
        for (j in seq_len(length(samples))) {
          j = 1
          alt = paste(samples[j], "alt", sep="_")
          ref = paste(samples[j], "ref", sep="_")
          
          vaf <- SNP_temp[[alt]] / (SNP_temp[[alt]] + SNP_temp[[ref]])
          # title = paste(output[i,]$chrom, output[i,]$start, output[i,]$end, samples[j], sep="_")
          title = samples[j]
          
          # ks = ks.test(vaf_germline, vaf)
          ks = is.unimodal(vaf)
          
          if (is.unimodal(vaf)) { # if similar to germline distribution, alt and ref will be the sum of all
            tcn_alt[i,j] = sum(SNP_temp[[alt]])
            tcn_ref[i,j] = sum(SNP_temp[[ref]])
            plot(density(vaf), xlim=c(0,1), main = title)
            qqnorm(vaf, main = title)
            qqline(vaf, col="green")
          } else { # else, cluster the vaf into two and take one cluster
            kmeans_result <- kmeans(vaf, centers = 2)
            tcn_alt[i,j] = sum(SNP_temp[[alt]][kmeans_result$cluster == 1])
            tcn_ref[i,j] = sum(SNP_temp[[ref]][kmeans_result$cluster == 1])
            plot(density(vaf), xlim=c(0,1), main = title)
            qqnorm(vaf, main = title)
            qqline(vaf, col="red")
          }
          
          ks_list <- c(ks_list, ks)
        }
        
        # print(ks_list)
        
        if (all(output[i,4:ncol(output)] > tcn_normal_range[1] & output[i,4:ncol(output)] < tcn_normal_range[2])) {
          if (all(!ks_list)) {
            to_keep_index <- c(to_keep_index, i)
          } 
        } else {
          to_keep_index <- c(to_keep_index, i)
        }
        
      }
    } else {
      # message("what if no SNP in the region; DROP")
    }
  }
  dev.off()
  
  print(length(to_keep_index))
}

#' #' get copy number cytoband and genes information
#' importCopyNumberInfo <- function(copy_number_file) {
#'   drivers = c()
#'   genes = c()
#'   cytoband = c()
#'   
#'   data <- read_csv(copy_number_file, show_col_types = FALSE)
#'   data$CNA = paste(data$chrom, data$start, data$end, sep = '-')
#'   
#'   for (i in 1:nrow(data)) {
#'     cytoband[data[i,]$CNA] <- data[i,]$cytoband
#'     genes[data[i,]$CNA] <- data[i,]$genes
#'     drivers[data[i,]$CNA] <- data[i,]$drivers
#'   }
#'   return(list(cytoband=cytoband, drivers=drivers, genes=genes))
#' }

#' import mutation file
importMutationFile <- function(mutation_file, alt_reads_thresh = 0, vaf_thresh = 0) {
  data <- read_csv(mutation_file, show_col_types = FALSE)
  data <- data %>% filter(alt_reads/total_reads>=vaf_thresh, alt_reads>=alt_reads_thresh)
  output_data <- list()

  output_data$y <- as.matrix(data[c("mutation", "sample", "alt_reads")] %>% pivot_wider(names_from = sample, values_from = alt_reads, values_fill = 0))
  rownames(output_data$y) <- output_data$y[,'mutation']
  output_data$y <- output_data$y[,-1, drop=FALSE]
  rowname = rownames(output_data$y)
  colname = colnames(output_data$y)
  output_data$y <- matrix(as.numeric(output_data$y), ncol = ncol(output_data$y))
  rownames(output_data$y) = rowname
  colnames(output_data$y) = colname

  output_data$MutID <- rowname

  output_data$n <- as.matrix(data[c("mutation", "sample", "total_reads")] %>% pivot_wider(names_from = sample, values_from = total_reads, values_fill = 0))
  rownames(output_data$n) <- output_data$n[,'mutation']
  output_data$n <- output_data$n[,-1, drop=FALSE]
  rowname = rownames(output_data$n)
  colname = colnames(output_data$n)
  output_data$n <- matrix(as.numeric(output_data$n), ncol = ncol(output_data$n))
  rownames(output_data$n) = rowname
  colnames(output_data$n) = colname

  if (any((output_data$y - output_data$n) > 0)) {
    warning("Total read count must be equal or bigger than alt read count. Please check input data before proceeding!")
    stop()
  }

  if (any(output_data$n==0)) {
    print("Total read counts of 0 encoutered. Replaced 0 with mean total read count.")
    output_data$n[output_data$n==0] <- round(mean(output_data$n))
  }

  output_data$S = ncol(output_data$y)
  output_data$I = nrow(output_data$y)

  mutation_position = unique(data[c("mutation", "chrom", "start", "end")])
  if (!all(sort(unique(data[c("mutation", "chrom", "start", "end")]) %>% pull(mutation)) == sort(output_data$MutID))) {
    warning("Some mutations may have duplicated chromosome information; keeping the first occurence.")
    mutation_position = mutation_position[match(unique(mutation_position$mutation), mutation_position$mutation), ]
  }

  output_data$position = mutation_position
  return(output_data)
}

simulation <- function(mcf=0.9, icn=2, minor_cn=1, depth=30, num_SNP=30, seed=123, normal_prop=1) {
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
  depth_neutral = rpois(n = num_neutral, lambda = depth)
  SNP_depth_neutral = numeric(num_neutral)
  for (i in seq_len(num_neutral)) {
    SNP_depth_neutral[i] <- rbinom(n=1, size=depth_total[i], prob=0.5)
  }
  
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
  plot(density(vaf_tumor), xlim=c(0,1), main = "VAF_tumor")
  qqnorm(vaf_tumor, main = "tumor")
  qqline(vaf_tumor, col="grey")
  
  um = is.unimodal(vaf_tumor)
  return(um)
}

