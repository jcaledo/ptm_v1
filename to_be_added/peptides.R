

## ---------- peptides.R ------------ ##
#                                      #
#     dipept                           #
#     oligo.pept                       #
#     seq2mass                         #
#                                      #
## ---------------------------------- ##


## ---------------------------------------------------------------- ##
#      dipept <- function(target, aa = 'M', uniprot = TRUE)          #
## ---------------------------------------------------------------- ##
#' Can random explain the ocurrence of dipeptides?
#' @description Contrast the null hypothesis that the occurrence of dipeptides is due to chance
#' @usage dipept <- function(target, aa = 'M', uniprot = TRUE)
#' @param target a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param aa one letter code of the amino acid to be analyzed.
#' @param uniprot logical, if TRUE the argument 'target' shoud be an ID.
#' @details The code computes the random variable X (number of occurrences of the dipeptide) and contrats that value with the expected value under the null hypothesis: X ~ Bin(p = (m/N)^2, n), where 'p' is the squared frequency of the amino acid 'aa' in that protein, 'm' is the frequency of 'aa', and n = N - 1 (being N the length of the protein). The function computes the expected value (E[X]) and the p-Value (the probability of obtaining a number of occurrences equal or higher to the observed value).
#' @return Returns a named vector with the protein id, N, m, Xo, Xe, p-character vector or a as a character string.
#' @author Juan Carlos Aledo
#' @examples dipep('P01009')
#' dipept(target = 'P00004', aa = 'G')
#' dipept(target = 'P0DP23')
#' dipept(target = "ASSSPASTSSV", aa = 'S', uniprot = FALSE)
#' @export

dipept <- function(target, aa = 'M', uniprot = TRUE){

  if (uniprot){
    seq <- ptm::get.seq(id = target) # Protein sequence
    id <- target
  } else {
    seq <- target
    id <- 'user-provided sequence'
  }

  N <- nchar(seq) # Protein length
  m <- length(gregexpr(aa, seq)[[1]]) # Number of residues aa
  p <- (m/N)^2 # Probability of aa-aa dipeptide under the null hypothesis
  n <- N - 1 # Number of different ways of observing aa-aa

  ## ---- Computing the ocurrence of dipeptides ------ ##
  sseq <- strsplit(seq, split = "")[[1]] # Sequence as vector
  Xo <- 0
  for (i in 1:n){
    if (sseq[i] == aa & sseq[i+1] == aa){
      Xo <- Xo + 1
    }
  }

  ## ---- Expectation according to null hypothesis ------ ##
  Xe <- round(n*p, 3) # Expectation of X ~ Bi(p, n)
  p_value <- 1 - pbinom(Xo - 1, n, p) # Probability of finding Xo or higher

  ## -------------------- Output ------------------------ ##
  output <- c(id, N, m, Xo, Xe, round(p_value, 3))
  names(output) <- c("id", "N", 'm',"Xo", "Xe", "p-value")
  return(output)
}

## ---------------------------------------------------------------- ##
#      oligo.pept <- function(target, aa = 'M', uniprot = TRUE)      #
## ---------------------------------------------------------------- ##
#' Maximum length of oligopeptides
#' @description Computes the maximum length of the indicated oligopeptide found in a sequence
#' @usage oligo.pept <- function(target, aa = 'M', uniprot = TRUE)
#' @param target a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param aa one letter code of the amino acid to be analyzed.
#' @param uniprot logical, if TRUE the argument 'target' shoud be an ID.
#' @details The code computes the maximum length of oligopeptides of the indicated amino acid.
#' @return Returns a positive integer.
#' @author Juan Carlos Aledo
#' @examples dipep('P01009')
#' dipept(target = 'P00004', aa = 'G')
#' dipept(target = 'P0DP23')
#' dipept(target = "ASSSPASTSSV", aa = 'S', uniprot = FALSE)
#' @export

oligo.pept <- function(target, aa = 'M', uniprot = TRUE){

  if (uniprot){
    seq <- ptm::get.seq(id = target) # Protein sequence
    id <- target
  } else {
    seq <- target
    id <- 'user-provided sequence'
  }
  motif <- paste(aa, "*", sep = "")
  t <- gregexpr(motif, seq)[[1]]
  return(max(attributes(t)$match.length))
}

## ---------------------------------------------------------------- ##
#         seq2mass <- function(seq)                                  #
## ---------------------------------------------------------------- ##
#' Find the mass of a peptide
#' @description Computes the mass of a given peptide
#' @usage seq2mass <- function(seq)
#' @param seq a character string specifying the sequence of the peptide.
#' @details
#' @return Returns the peptide mass in Da.
#' @author Juan Carlos Aledo
#' @examples seq2mass('MPSSVSWGILLLAGLCCLVPVSLAEDPQGDAAQK')
#' @export

seq2mass <- function(seq){

  a <- aa.comp(seq, uniprot = FALSE)

  w <- c(71.079018, 157.196106, 114.104059, 114.080689,
         103.143407, 128.131048, 128.107678, 57.05203, 137.141527,
         113.159985, 113.159985, 129.18266, 131.197384, 147.177144,
         97.117044, 87.078323, 101.105312, 186.213917, 163.176449,
         99.132996)
  a$contribution <- a$frequency * w
  output <- sum(a$contribution)
  attr(output, 'unit') <- 'Da'
  return(output)
}
