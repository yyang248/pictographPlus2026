#' read input data file and store in PICTOGRAPH input format
#' @export
#' @param mutation_file mutation data file that contains columns "sample", "mutation", "chrom", "start", "end", total_reads", and "alt_reads";
#' @param copy_number_file copy number file that contains columns "sample", "chrom", "start", "end", "major", "minor"
importFiles <- function(mutation_file, copy_number_file=NULL, alt_reads_thresh = 0, vaf_thresh = 0) {
  mutation_data = importMutationFile(mutation_file, alt_reads_thresh, vaf_thresh)
  copy_number_data = importCopyNumberFile(copy_number_file)
  warning("CNA not checked for overlap yet; user need to make sure CNA seperate")
  mutation_data$major = copy_number_data$major
  mutation_data$minor = copy_number_data$minor
  mutation_data$tcn = mutation_data$major + mutation_data$minor
  mutation_data$overlap = resolveOverlap(mutation_data)

  return(mutation_data)
}

#' check whether a mutation overlaps a CNA region
resolveOverlap <- function(mutation_data) {
  mut_count = nrow(mutation_data$position)
  cna_count = nrow(mutation_data$major)
  output = matrix(0, nrow=mut_count, ncol=cna_count)
  rownames(output) = mutation_data$position$mutation
  colnames(output) = rownames(mutation_data$major)
  # print(output)
  cnas = rownames(mutation_data$major)
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

#' import copy number file
importCopyNumberFile <- function(copy_number_file) {
  data <- read_csv(copy_number_file, show_col_types = FALSE)
  data$CNA = paste(data$chrom, data$start, data$end, sep = '-')
  output_data <- list()
  output_data$major = as.matrix(data[c("sample", "CNA", "major")] %>% pivot_wider(names_from = sample, values_from = major, values_fill = 1))
  rownames(output_data$major) <- output_data$major[,'CNA']
  output_data$major <- output_data$major[,-1, drop=FALSE]
  rowname = rownames(output_data$major)
  colname = colnames(output_data$major)
  output_data$major <- matrix(as.numeric(output_data$major), ncol = ncol(output_data$major))
  rownames(output_data$major) = rowname
  colnames(output_data$major) = colname

  output_data$minor = as.matrix(data[c("sample", "CNA", "minor")] %>% pivot_wider(names_from = sample, values_from = minor, values_fill = 1))
  rownames(output_data$minor) <- output_data$minor[,'CNA']
  output_data$minor <- output_data$minor[,-1, drop=FALSE]
  rowname = rownames(output_data$minor)
  colname = colnames(output_data$minor)
  output_data$minor <- matrix(as.numeric(output_data$minor), ncol = ncol(output_data$minor))
  rownames(output_data$minor) = rowname
  colnames(output_data$minor) = colname

  return(output_data)
}

#' import mutation file
importMutationFile <- function(mutation_file, alt_reads_thresh = 0, vaf_thresh = 0) {
  data <- read_csv(mutation_file, show_col_types = FALSE)
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
    warning("Total read counts of 0 encoutered. Replaced 0 with mean total read count.")
    output_data$n[output_data$n==0] <- round(mean(output_data$n))
  }

  output_data$S = ncol(output_data$y)
  output_data$I = nrow(output_data$y)

  output_data$y[output_data$y / output_data$n < vaf_thresh] = 0
  output_data$y[output_data$y < alt_reads_thresh] = 0

  mutation_position = unique(data[c("mutation", "chrom", "start", "end")])
  if (!all(sort(unique(data[c("mutation", "chrom", "start", "end")]) %>% pull(mutation)) == sort(output_data$MutID))) {
    warning("Some mutations may have duplicated chromosome information; keeping the first occurence.")
    mutation_position = mutation_position[match(unique(mutation_position$mutation), mutation_position$mutation), ]
  }

  output_data$position = mutation_position
  return(output_data)
}

