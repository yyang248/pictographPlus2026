#' separate mutation and copy number data based on sample presence profile
#' @export
separateBySamplePresence <- function(input_data) {

  mutationTypes = separateMutationsBySamplePresence(input_data)
  # print(mutationTypes)

  cnaTypes = separateCNAsBySamplePresence(input_data)
  # print(cnaTypes)

  sep_list <- list()
  sep_list$mutation <- mutationTypes
  sep_list$cna <- cnaTypes
  print(sep_list)
  return(sep_list)
}

#' separate mutations by sample presence
separateMutationsBySamplePresence <- function(input_data) {
  # returns list of lists --
  # each item of list contains input data for a mutation sample presence set
  # original mutation indices from input_data are recorded in $mutation_indices
  pres <- ifelse(input_data$y > 0, 1, 0)
  pat <- apply(pres, 1, function(x) paste0(x, collapse=""))
  types <- sort(names(table(pat)), decreasing=TRUE)
  if (length(types) == 1) {
    type_indices <- list()
    type_indices[[types]] <- seq_len(input_data$I)
  } else {
    type_indices <- lapply(types, function(x) which(pat == x))
    names(type_indices) <- types
  }
  warning('consider merge one-sample cluster with other clusters')
  return(type_indices)
}

#' separate copy number alterations by sample presence
separateCNAsBySamplePresence <- function(input_data) {
  # returns list of lists --
  # each item of list contains input data for a mutation sample presence set
  # original mutation indices from input_data are recorded in $mutation_indices
  pres <- ifelse(input_data$tcn == 2, 0, 1)
  pat <- apply(pres, 1, function(x) paste0(x, collapse=""))
  types <- sort(names(table(pat)), decreasing=TRUE)
  if (length(types) == 1) {
    type_indices <- list()
    type_indices[[types]] <- seq_len(input_data$I)
  } else {
    type_indices <- lapply(types, function(x) which(pat == x))
    names(type_indices) <- types
  }
  return(type_indices)
}
