####
# Codon Pair Bias Analysis ####
####

#' getCodonPairBias
#'
#' Generate Codon pair bias measures for each codon pair and sequence
#'
#' @param pair_counts counts for each codon pair
#' @param single_counts counts for each codon
#' @param method method of codon pair bias calculation, e.g. RSCPU
#'
#' @export
#'
getCodonPairBias <- function(pair_counts, single_counts, method = "RSCPU") {

  # library
  require(Biostrings)

  # check if all codons are here
  GENETIC_CODE_PAIR <- getPairGENETIC_CODE()
  if(!all(colnames(pair_counts) %in% names(GENETIC_CODE_PAIR)))
    stop("Please make sure that your count data includes all codons")

  # choose a method for codon bias normalization
  if(method == "RSCPU"){
    normcounts <- getCodonBiasRSCPU(pair_counts, single_counts)
  }

  # return normalized counts
  return(normcounts)
}

# RSCU normalization
getCodonBiasRSCPU <- function(pair_counts, single_counts) {

  # get expected counts
  pair_counts_exp <- t(apply(single_counts, 1, function(x){
    as.vector(crossprod(as.matrix(x)))
  }))

  # get results
  norm_pair_counts <- pair_counts/pair_counts_exp

  # return
  return(norm_pair_counts)
}

####
# Codon Bias Analysis ####
####

#' getCodonBias
#'
#' Generate Codon bias measures for each codon and sequence
#'
#' @param counts counts for each codon
#' @param method method of codon pair bias calculation, e.g. RSCU
#' @param sizeFactor the size factor when size normalization is used
#'
#' @export
#'
getCodonBias <- function(counts, method = "RSCU", sizeFactor = 100) {

  # library
  require(Biostrings)

  # check if all codons are here
  if(!all(colnames(counts) %in% names(GENETIC_CODE)))
    stop("Please make sure that your count data includes all codons")

  # choose a method for codon bias normalization
  if(method == "RSCU"){
    normcounts <- getCodonBiasRSCU(counts)
  } else if(method == "CLR"){
    normcounts <- getCodonBiasCLR(counts)
  } else if(method == "Size"){
    normcounts <- getCodonBiasSize(counts, sizeFactor = sizeFactor)
  }

  # return normalized counts
  return(normcounts)
}

# RSCU normalization
getCodonBiasRSCU <- function(counts) {
  normcounts <- t(apply(counts, 1, function(x){
    agg_code <- aggregate(x, list(GENETIC_CODE), mean)
    agg_code_names <- agg_code$Group.1
    agg_code <- agg_code$x
    names(agg_code) <- agg_code_names
    return(x/agg_code[GENETIC_CODE])
  }))
  normcounts[is.nan(normcounts)] <- 1
  return(normcounts)
}

# CLR normalization
getCodonBiasCLR <- function(counts) {
  normcounts <- t(apply(counts, 1, function(x){
    return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
  }))
  return(normcounts)
}

# Simple size normalize
getCodonBiasSize <- function(counts, sizeFactor) {
  normcounts <- t(apply(counts, 1, function(x) return((x/sum(x))*sizeFactor)))
}
