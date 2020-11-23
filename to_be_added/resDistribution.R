## ---- resDistribution.R --------- ##
#                                    #
#       cdf.dbir                     #
#       res.barcode                  #
#       make.bins
#                                    #
## -------------------------------- ##

## --------------------------------------------------------------- ##
#     cdf.dbir(target, aa = 'M', uniprot = TRUE, p_value = 0.05)    #
## --------------------------------------------------------------- ##
#' CDF of The Distance Between Identical Residues
#' @description Returns the theoretical, empirical and simulated cumulative distribution function for the distances between identical residues
#' @usage cdf.dbir(target, aa = 'M', uniprot = TRUE, p_value = 0.05)
#' @param target a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param aa the amino acid of interest.
#' @param uniprot logical, if TRUE the argument 'target' shoud be an ID.
#' @param p_value level of significance.
#' @details The null hypothesis is that the residues of the amino acid 'aa' are randomly distributed through the sequence. Thus, if we have 'n' residues of the amino acid 'aa' randomly distributed along a sequence of length 'N', to analyse the observed spatial pattern as a function of the distance, d, between two identical residues, we define a random variable, X, as the number of residue different to 'aa' found between two residues of type 'aa'. Under the conditions of the null hypothesis, the theoretical CDF is that corresponding to the geometric distribution: G(d) = 1 - (1-p)^{|d|}, where p = n/N.
#' @return Returns a plot with the theoretical (dashed line), the empirical (red line) and the simulated (yellow lines) CDFs. In addition, the function returns a dataframe with the X values that disagree with the null hypothesis, as well as their corresponding empirical CDF values.
#' @author Juan Carlos Aledo
#' @examples cdf.dbir('P01009')
#' @seealso res.barcode()
#' @export

cdf.dbir <- function(target, aa = 'M', uniprot = TRUE, p_value = 0.05){

  if (uniprot){
    seq <- ptm::get.seq(id = target) # Protein sequence
    id <- target
  } else {
    seq <- target
    id <- 'user-provided sequence'
  }

  at <- gregexpr(aa, seq)[[1]] # positions at which aa is found
  N <- nchar(seq) # sequence's length
  n <- length(at) # number of aa residues
  p <- n/N # relative frequency of aa in the sequence


  ## --------- Theoretical distribution
  plot(0:N, pgeom(0:N, p), ty = 'l', lty = 2, xlab = 'd', ylab = 'G(d) = P[X <= d]')

  ## --------- Empirical distribution
  G <- as.data.frame(matrix(rep(NA, n*4), ncol = 4))
  names(G) <- c('i', 'pos', 'd_CT', 'd')
  G$i <- 1:n
  G$pos <- at
  for (i in 1:nrow(G)){
    G$d_CT[i] <- N - G$pos[i]
    t <- at[i+1] - at[i]
    G$d[i] <- t - 1 # number of fails instead of number trials
  }
  G <- G[-nrow(G),] # last aa is removed

  numerador <- c()
  eCDF <- c() # empirical CDF
  for (d in 0:N){
    counter <- 0
    for (i in 1:nrow(G)){
      if (G$d[i] <= d){
        counter <- counter + 1
      }
    }
    numerador <- c(numerador, counter)
  }
  eCDF <- numerador/nrow(G)

  ## ------ Simulated distributions
  set.seed(123)
  Simul <- as.data.frame(matrix(rep(NA, 100 * (N+1)), ncol = N+1))
  names(Simul) <- 0:N

  for (s in 1:100){
    Gs <- as.data.frame(matrix(rep(NA, n*3), ncol = 3))
    names(Gs) <- c('i', 'rpos', 'd')
    Gs$i <- 1:n
    rpos <- sample(1:N, size = n, replace = FALSE)
    rpos <- rpos[order(rpos)]

    for (i in 1:nrow(Gs)){
      t <- rpos[i+1] - rpos[i]
      Gs$d[i] <- t - 1 # number of fails instead of number trials
    }
    Gs$rpos <- rpos
    Gs <- Gs[-nrow(Gs),] # last aa is removed

    numerador <- c()
    sCDF <- c() # simulated CDF
    for (d in 0:N){
      counter <- 0
      for (i in 1:nrow(Gs)){
        if (Gs$d[i] <= d){
          counter <- counter + 1
        }
      }
      numerador <- c(numerador, counter)
    }
    sCDF <- numerador/nrow(Gs)
    points(0:(length(sCDF)-1), sCDF, ty = 'l', lty = 2, col = 'yellow')
    Simul[s,] <- sCDF
  }

  output <- as.data.frame(matrix(rep(NA, (N+1) * 4), ncol = 4))
  colnames(output) <- c('X','eCDF', 'higher_or_equal', 'lower_or_equal')
  output[,1] <- 0:N
  output[,2] <- eCDF
  for (i in 1:length(eCDF)){
    output[i,3] <- length(Simul[which(Simul[,i] >= eCDF[i])])
    output[i,4] <- length(Simul[which(Simul[,i] <= eCDF[i])])
  }

  points(0:(length(eCDF)-1), eCDF, ty = 'l', lwd = 2, col = 'red')

  ## --- Significant points
  output <- output[which(output$higher_or_equal <= 100 * p_value | output$lower_or_equal <= 100 * p_value),]

  he <- output[which(output$higher_or_equal <= 100 * p_value),  c(1,2)]
  points(he$X, he$eCDF, ty = 'p')

  le <- output[which(output$lower_or_equal <= 100 * p_value),  c(1,2)]
  points(le$X, le$eCDF, ty = 'p')
  points(0:N, pgeom(0:N, p), ty = 'l', lty = 2, xlab = 'd', ylab = 'G(d) = P[X <= d]')

  attr(output, 'seq') <- id
  return(output)

}


## --------------------------------------------------------------- ##
#     res.barcode(target, aa = 'M', uniprot = TRUE, title = FALSE   #
## --------------------------------------------------------------- ##
#' Residue Barcode
#' @description Barcode representation of the distribution along the sequence of a given amino acid
#' @usage res.barcode(target, aa = 'M', uniprot = TRUE, title = FALSE)
#' @param target a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param aa the amino acid of interest.
#' @param uniprot logical, if TRUE the argument 'target' shoud be an ID.
#' @param title logical, if TRUE additional data are shown in the plot title.
#' @details The null hypothesis is that the residues of the amino acid 'aa' are randomly distributed through the sequence. Thus, if we have 'n' residues of the amino acid 'aa' randomly distributed along a sequence of length 'N', to analyse the observed spatial pattern as a function of the distance, d, between two identical residues, we define a random variable, X, as the number of residue different to 'aa' found between two residues of type 'aa'. Under the conditions of the null hypothesis, the theoretical CDF is that corresponding to the geometric distribution: G(d) = 1 - (1-p)^{|d|}, where p = n/N.
#' @return Returns a barcode plot as well as the corresponding vector of zeros (a residue different to 'aa') and ones (residue equal to 'aa').
#' @author Juan Carlos Aledo
#' @examples res.barcode('P01009')
#' @seealso cdf.dbir()
#' @importFrom bio3d aa123
#' @export

res.barcode <- function(target, aa = 'M', uniprot = TRUE, title = FALSE){

  if (uniprot){
    seq <- ptm::get.seq(id = target) # Protein sequence
    id <- target
  } else {
    seq <- target
    id <- 'user-provided sequence'
  }

  at <- gregexpr(aa, seq)[[1]] # positions at which aa is found
  barcode <- rep(0, nchar(seq))
  barcode[at] <- 1

  f <- round(length(at)/length(barcode), 3)

  if (title){
    res <- bio3d::aa123(aa)
    plot(1:length(barcode), barcode, ty = 'h', xlab = 'Residue Position', ylab = '',
         yaxt = 'n', main = paste(res, '-', 100 * f, '%', sep = ''))
  } else {
    plot(1:length(barcode), barcode, ty = 'h', xlab = 'Residue Position', ylab = '', yaxt = 'n')
  }


  attr(barcode, 'seq') <- id
  attr(barcode, 'residue' ) <- aa
  attr(barcode, 'frequency') <- f

  return(barcode)
}

## --------------------------------------------------------------- ##
#     make.bins(target, aa = 'M', uniprot = TRUE, w)   #
## --------------------------------------------------------------- ##
#' Segment a Protein Sequence in Bins
#' @description Segments a protein sequence in bins of the wished size and computes the occurrence of a given amino acid in each bins.
#' @usage make.bins(target, aa = 'M', uniprot = TRUE)
#' @param target a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param aa the amino acid of interest.
#' @param uniprot logical, if TRUE the argument 'target' shoud be an ID.
#' @param w widow size.
#' @details
#' @return Returns a dataframe with the numbers of occurrences of the given amino acid in each segment.
#' @author Juan Carlos Aledo
#' @examples make.bins('P01009')
#' @seealso aa.aggr()
#' @export

make.bins <- function(target, aa = 'M', uniprot = TRUE, w){

  if (uniprot){
    seq <- ptm::get.seq(id = target) # Protein sequence
    id <- target
  } else {
    seq <- target
    id <- 'user-provided sequence'
  }
  n <- length(gregexpr(aa, seq)[[1]]) # number of aa residues
  N <- nchar(seq) # protein length
  p <- n/N # frequency of aa
  seq_ <- substr(seq, start = (N %% w) + 1, stop = N) # remove from NT the leftover residues

  output<- as.data.frame(matrix(rep(NA, (N %/% w)*2), ncol = 2))
  names(output) <- c('segment', 'counts')
  output$segment <- 1:nrow(output)
  counter <- 0

  for (i in seq(1, nchar(seq_), by = w)){
    counter <- counter + 1
    t <- substr(seq_, start = i, stop = i + w - 1)
    occurrences <- gregexpr(aa, t)[[1]]
    if (occurrences[1] == -1){
      output$counts[counter] <- 0
    } else {
      output$counts[counter] <- length(occurrences)
    }
  }

  attr(output, 'mean') <- mean(output$counts)
  attr(output, 'variance') <- sd(output$counts)^2

  return(output)
}

## --------------------------------------------------------------- ##
#     aa.aggr(target, aa = 'M', uniprot = TRUE, up = 100)           #
## --------------------------------------------------------------- ##
#' Search for The Trend of a Given Amino Acid to Aggregate in The Primary Structure
#' @description Plots the Variance/Mean for the occurrence of an amino acid as function of the segment lengths.
#' @usage make.bins(target, aa = 'M', uniprot = TRUE)
#' @param target a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param aa the amino acid of interest.
#' @param uniprot logical, if TRUE the argument 'target' shoud be an ID.
#' @param up segments of sizes ranging from 1 to up.
#' @details
#' @return Returns a plot as well as the minimum and maximum of the Variance/Mean ratio.
#' @author Juan Carlos Aledo
#' @examples aa.aggr('P01009')
#' @seealso make.bins(), sliding.window()
#' @export

aa.aggr <- function(target, aa = 'M', uniprot = TRUE, up = 100){

  if (uniprot){
    seq <- ptm::get.seq(id = target) # Protein sequence
    id <- target
  } else {
    seq <- target
    id <- 'user-provided sequence'
  }



}


## --------------------------------------------------------------- ##
#     sliding.window(target, aa = 'M', uniprot = TRUE, sw = 10)     #
## --------------------------------------------------------------- ##
#' Search for the region in the primary structure enriched with a given amino acid
#' @description Plots the occurrences of a given amino acid withind a sliding window.
#' @usage sliding.window(target, aa = 'M', uniprot = TRUE, sw = 10)
#' @param target a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param aa the amino acid of interest.
#' @param uniprot logical, if TRUE the argument 'target' shoud be an ID.
#' @param sw size of the window
#' @details
#' @return Returns a plot as well as start and stop positions for the window with maximal occurrences of the given amino acid.
#' @author Juan Carlos Aledo
#' @examples sliding.window('P01009')
#' @seealso make.bins(), sliding.window()
#' @export

sliding.window <- function(target, aa = 'M', uniprot = TRUE, sw = 10){

  if (uniprot){
    seq <- ptm::get.seq(id = target) # Protein sequence
    id <- target
  } else {
    seq <- target
    id <- 'user-provided sequence'
  }

  N <- nchar(seq) # protein length
  up <- N - sw + 1 # start position for the last segment

  output <- as.data.frame(matrix(rep(NA, up * 3), ncol = 3))
  names(output) <- c("start", "stop", "occurrences")

  for (i in 1:up){
    output$start[i] <- i
    output$stop[i] <- i + sw - 1
    t <- gregexpr(aa, substr(seq, start = i, stop = i + sw -1))[[1]]
    if (t[1] == -1){
      output$occurrences[i] <- 0
    } else {
      output$occurrences[i] <- length(t)
    }
  }

 t <- gregexpr(aa, seq)[[1]]
 if (t[1] != -1){
   first <- t[1]
   last <- t[length(t)]
   range <- last - first
 } else {
   first <- "aa not found"
   last <- "aa not found"
   range <- "aa not found"
 }

 attr(output, "seq") <- id
 attr(output, "aa") <- aa
 attr(output, "first") <- first
 attr(output, "last") <- last
 attr(output, "range") <- range

 return(output)
}
