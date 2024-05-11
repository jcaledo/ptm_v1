## ---------- environment.R ------------ ##
#                                         #
#     env.extract                         #
#     env.matrices                        #
#     env.Ztest                           #
#     env.plot                            #
#                                         #
## ------------------------------------- ##


## ------------------------------------------------------------------------------- ##
#   env.extract <- function(prot, db = 'none', c, r, ctr = 'none', exclude = c())   #
## ------------------------------------------------------------------------------- ##
#' Sequence Environment Around a Given Position
#' @description Extracts the sequence environment around a given position.
#' @usage env.extract(prot, db = 'none', c, r, ctr = 'none', exclude = c())
#' @param prot either a uniprot id or a string sequence.
#' @param db a character string specifying the desired database; it must be one of 'uniprot', 'metosite', 'none'.
#' @param c center of the environment.
#' @param r radius of the environment.
#' @param ctr the type of control environment; it must be one of 'random', 'closest', or 'none'.
#' @param exclude a vector containing the positions to be excluded as control.
#' @details The random control returns an environment center at a random position containing the same type or amino acid than the positive environment. The closest control searches for the closest position where such a type of amino acid is found and returns its environment.
#' @return Returns a  list of two strings (environments).
#' @author Juan Carlos Aledo
#' @examples env.extract('P01009', db = 'uniprot', 271, 10, ctr = 'random')
#' @references Aledo et al. Sci Rep. 2015; 5: 16955. (PMID: 26597773)
#' @seealso env.matrices(), env.Ztest() and env.plot()
#' @export

env.extract <- function(prot, db = 'none', c, r, ctr = 'none', exclude = c()){

  ## ------------ Preparing the sequence ------------ ##
  if (db == 'uniprot'){
    seq <- tryCatch(
      {
        get.seq(prot)
      },
      error = function(cond){
        return(NULL)
      }
    )
    if (is.null(seq)){
      message("Sorry, get.seq failed")
      return(NULL)
    }

  } else if (db == 'metosite'){
    seq <- tryCatch(
      {
        get.seq(prot, db = 'metosite')
      },
      error = function(cond){
        return(NULL)
      }
    )
    if (is.null(seq)){
        message("Sorry, get.seq failed")
        return(NULL)
    }
  } else {
    seq <- prot
  }

  ## ---- Checking for coherence of the request ---- ##
  if (c > nchar(seq)){
    stop("The requested center is not found in the provided sequence")
  }

  ## -- Ancillary function to isolate environment -- ##
  extract <- function(seq, c, r){
    envL <- substring(seq, c-r, c-1)
    if (nchar(envL) < r){
      X <- paste(rep('X', r-nchar(envL)), collapse = "")
      envL <- paste(X, envL, sep = "")
    }

    envR <- substring(seq, c+1, c+r)
    if (nchar(envR) < r){
      X <- paste(rep('X', r-nchar(envR)), collapse = "")
      envR <- paste(envR, X, sep = "")
    }

    envC <- tolower(substring(seq, c, c))
    env <- paste(envL, envC, envR, sep = "")
    return(env)
  }

  ## ------- Building the output file ------- ##

  ## -- The positive:
  env_pos <- extract(seq, c, r)

  ## -- The control:
  res <- aa.at(c, seq, FALSE)
  all_res <- gregexpr(res, seq)[[1]] # Positions containing the same amino acid finds at 'c'
  other_res <- all_res[-which(all_res == c)] # Positions containing the same amino acid except the positive
  other_res <- setdiff(other_res, exclude) # Exclude some positions if requested
  n <- length(other_res)

  if (n == 0){ # There are no other residues in the protein to be used as control
    env_ctr <- ""
  } else {
    if (ctr == 'random'){
      c_random <- other_res[sample(1:n, 1)] # Random position where there is a residue equal to that find at 'c'
      env_ctr <- extract(seq, c_random, r)
    } else if (ctr == 'closest'){
      closest <- which(abs(other_res - c) == min(abs(other_res - c)))[1] # there may be two equidistant residues
      c_closest <- other_res[closest]
      env_ctr <- extract(seq, c_closest, r)
    } else {
      env_ctr <- ""
    }
  }

  output <- list(env_pos, env_ctr)
  names(output) <- c("Positive", "Control")

  return(output)
}


## ---------------------------------------------------------------- ##
#               env.matrices <- function(env)                        #
## ---------------------------------------------------------------- ##
#' Environment Matrices
#' @description Provides the frequencies of each amino acid within the environment.
#' @usage env.matrices(env)
#' @param env a character string vector containing the environments.
#' @return Returns a list of two dataframes. The first, shown the environment in matrix form. The second provides the frequencies of each amino acid within the environments.
#' @author Juan Carlos Aledo
#' @references Aledo et al. Sci Rep. 2015; 5: 16955. (PMID: 26597773)
#' @examples env.matrices(c('ANQRmCTPQ', 'LYPPmQTPC', 'XXGSmSGXX'))
#' @seealso env.extract(), env.Ztest() and env.plot()
#' @export

env.matrices <- function(env){

  ## --------- Checking input data --------- ##
  env <- env[!is.na(env)] # remove NA if needed
  L <- unique(nchar(env)) # environment's length
  if (length(L) != 1){
    stop("The length of sequence environments is not homogeneous in your input data")
  }
  # non-canonical amino acid to X:
  env <- gsub("[^A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y,
               a, c, d, e, f, g, h, i, k, l, m, n, p, q, r, s, t, v, w, y]",
              "X", env)

  aa <- c("A","C","D","E","F","G","H","I","K","L","M",
          "N","P","Q","R","S","T","V","W","Y", "X")
  n_env <- length(env) # number of environments being analyzed

  r <- L %/% 2 # environment's radius

  ## ---- Dealing with the amino acids DF ---- ##
  g <- mapply(strsplit,  env, "")
  aaDF <- data.frame(do.call(rbind, g))
  names(aaDF) <- -r:r
  aaDF_ <- as.data.frame(apply(aaDF, 2, toupper)) # working copy
  aaDF_ <- as.data.frame(lapply(aaDF_, factor, levels = aa))

  ## ---- Dealing with the frequencies DF ---- ##
  fDF <- as.data.frame(matrix(rep(NA, L*21), ncol = L))
  names(fDF) <- -r:r
  rownames(fDF) <- aa

  for (i in 1:L){
    fDF[,i] <- table(aaDF_[,i])
  }
  t <- toupper(aaDF[,r+1][1]) # central residue
  if (fDF[which(aa == t), r+1] != n_env){
    warn <- paste("The amino acid found at the central position of the environment \n",
                  "may not be always the same in your input data", sep = "")
    warning(warn)
  }
  attr(fDF, 'number_sequences_analyzed') <- n_env
  output <- list(aaDF, fDF)

  return(output)
}


## ---------------------------------------------------------------- ##
#            env.Ztest <- function(pos, ctr, alpha = 0.05)           #
## ---------------------------------------------------------------- ##
#' Preferred/Avoided Amino Acids Within an Environment
#' @description Searches for amino acids either overrepresented or avoided at given positions within a sequence environment.
#' @usage env.Ztest(pos, ctr, alpha = 0.05)
#' @param pos a 21 x m matrix containing the absolute frequencies of 21 amino acids at the m positions, in the positive environments.
#' @param ctr a 21 x m matrix containing the absolute frequencies of 21 amino acids at the m positions, in the control environments.
#' @param alpha significance level.
#' @details Please, note that in addition to the 20 proteinogenetic amino acid we are using the symbol X when the target (central) residue is closer to the N-terminal or C-terminal of the protein than the radius used.
#' @return Returns a list with three elements: (1) a matrix with the values of the Z statistical. (2) A dataframe with information regarding amino acid overrepresented in the positive environments, and (3) a dataframe similar to the previous one, but for amino acids avoided from the positive environments.
#' @author Juan Carlos Aledo
#' @examples pos = env.matrices(hmeto$positive)[[2]][,-11]; ctr = env.matrices(hmeto$control)[[2]][,-11]
#' @examples env.Ztest(pos, ctr, alpha = 0.0001)
#' @references Aledo et al. Sci Rep. 2015; 5: 16955. (PMID: 26597773)
#' @seealso env.extract(), env.matrices() and env.plot()
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @importFrom stats complete.cases
#' @export

env.Ztest <- function(pos, ctr, alpha = 0.05){

  ## --- Absolute to relative frequencies --- ##
  if (!is.null(attributes(pos)$number_sequences_analyzed)){
    n_pos <- attributes(pos)$number_sequences_analyzed
  } else {
    n_pos <- sum(pos[,1])
  }
  if (!is.null(attributes(ctr)$number_sequences_analyzed)){
    n_ctr <- attributes(ctr)$number_sequences_analyzed
  } else {
    n_ctr <- sum(ctr[,1])
  }

  pos <- as.matrix(pos)/n_pos
  ctr <- as.matrix(ctr)/n_ctr
  ## ------------ Data quality control ------ ##
  poor <- pos[which(pos > 1)]
  if (length(poor) > 0){
    pos[which(pos > 1)] <- 1
    warning(paste("The value ", poor, " has been approached to 1.0"))
  }
  ## ------------ Normalization ------------- ##
  difference <- pos - ctr # frequency pos - frequency ctr
  var_difference <- (pos*(1-pos)/n_pos) + (ctr*(1-ctr)/n_ctr)
  var_difference[which(var_difference == 0)] <- .Machine$double.xmin # When there is not variance
  Z <- difference/sqrt(var_difference)

  critical_value <- stats::qnorm((1-alpha))

  ## -----------  Overrepresented ----------- ##
  over <- as.data.frame(matrix(rep(NA, 400*4), ncol = 4))
  names(over) <- c('aa', 'at', 'Z', 'pValue')
  counter <- 1
  n_positions <- dim(Z)[2] # number of positions within the environment

  for (j in 1:n_positions){
    for (i in 1:21){ # number of amino acids
      if (Z[i,j] > critical_value){
        over$aa[counter] <- rownames(Z)[i]
        over$at[counter] <- colnames(Z)[j]
        over$Z[counter] <- Z[i,j]
        over$pValue[counter] <- 1 - stats::pnorm(Z[i,j])
        counter <- counter + 1
      }
    }
  }
  over <- over[complete.cases(over),]
  over <- over[order(over$Z, decreasing = TRUE),]

  ## ----------  Underrepresented ----------- ##
  under <- as.data.frame(matrix(rep(NA, 400*4), ncol = 4))
  names(under) <- c('aa', 'at', 'Z', 'pValue')
  counter <- 1
  for (j in 1:n_positions){
    for (i in 1:21){ # number of amino acids
      if (Z[i,j] < -critical_value){
        under$aa[counter] <- rownames(Z)[i]
        under$at[counter] <- colnames(Z)[j]
        under$Z[counter] <- Z[i,j]
        under$pValue[counter] <- stats::pnorm(Z[i,j])
        counter <- counter + 1
      }
    }
  }
  under <- under[complete.cases(under),]
  under <- under[order(under$Z, decreasing = FALSE),]

  ## ----------  Output ----------- ##
  results <- list(Z, over, under)
  return(results)
}


## ----------------------------------------------------------------- ##
#    env.plot <- function(Z, aa, pValue = 0.001, ylim = c(-8,6),      #
#                                      ty = 'p', title = "")          #
## ----------------------------------------------------------------- ##
#' Differential Sequence Environment Plot
#' @description Plots the Z statistics at each position within the environment for the requested amino acid.
#' @usage env.plot(Z, aa, pValue = 0.001, ylim = c(-8,6), ty = 'p',title = "")
#' @param Z a matrix containing the standardized difference in frequencies (positive - control).
#' @param aa the amino acid of interest.
#' @param pValue the p-Value chosen to confer statistical significance.
#' @param ylim range of the dependent variable. If we pass the argument 'automatic', the function will choose a suitable range for you.
#' @param ty what type of plot should be drawn ("p": points, "l": lines, "b": both).
#' @param title character string giving a title for the plot.
#' @details  The p-Value is used to draw two horizontal lines delimiting the region supporting the null hypothesis: no significant differences. Points laying above or below of these lines cannot be explained by randomness.
#' @return This function returns a plot for the requested amino acid.
#' @author Juan Carlos Aledo
#' @examples \dontrun{## Get the matrices
#' pos = env.matrices(hmeto$positive)[[2]][,-11]
#' ctr = env.matrices(hmeto$control)[[2]][,-11]
#' ## Run the test
#' Z = env.Ztest(pos, ctr, alpha = 0.0001)[[1]]
#' ## Plot the results
#' env.plot(Z, aa = 'E', pValue = 0.05)}
#' @references Aledo et al. Sci Rep. 2015; 5: 16955. (PMID: 26597773)
#' @seealso env.extract(), env.matrices() and env.Ztest()
#' @importFrom stats qnorm
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @export

env.plot <- function(Z, aa, pValue = 0.001, ylim = c(-8,6), ty = 'p', title = ""){

  amino_acids <- c("A","C","D","E","F","G","H","I","K","L","M",
          "N","P","Q","R","S","T","V","W","Y", "X")

  cv <- stats::qnorm(1 - pValue) # critical value
  AA <- which(amino_acids == aa) # alphabetical order of the amino acid
  positions <- as.numeric(colnames(Z))

  if (ylim[1] == 'automatic'){
    u <- max(Z[AA, ], cv)
    l <- min(Z[AA, ], -cv)
    ylim <- c((l + 0.1*l), (u + 0.1*u))
  }
  plot(positions, Z[AA,], ylab="Z value", ylim = ylim, ty = ty )
  title(main = aa, sub = title)
  abline(-cv, 0, lty=2)
  abline(cv, 0, lty=2)
}
