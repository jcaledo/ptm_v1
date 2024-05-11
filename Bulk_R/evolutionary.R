## ----------- envolutionary.R ------------ ##
#                                            #
#     msa                                    #
#                                            #
## ---------------------------------------- ##

## ----------------------------------------------------------------- ##
#   msa <- function(sequences, ids = names(sequences),                #
#                     sfile = FALSE, inhouse = FALSE)                 #
## ----------------------------------------------------------------- ##
#' Multiple Sequence Alignment
#' @description Aligns multiple protein sequences.
#' @usage msa(sequences, ids = names(sequences), sfile = FALSE, inhouse = FALSE)
#' @param sequences vector containing the sequences.
#' @param ids vector containing the sequences' ids.
#' @param sfile if different to FALSE, then it should be a string indicating the path where to save a fasta alignment file.
#' @param inhouse logical, if TRUE the in-house MUSCLE software is used. It must be installed on your system and in the search path for executable.
#' @return Returns a list of four elements. The first one ($seq) provides the sequences analyzed, the second element ($ids) returns the identifiers, the third element ($aln) provides the alignment in fasta format and the fourth element ($ali) gives the alignment in matrix format.
#' @examples \dontrun{msa(sequences = sapply(c("P19446", "P40925", "P40926"), ptm::get.seq),
#'  ids = c("wmelon", "cyt", "mit"))}
#' @references Edgar RC. Nucl. Ac. Res. 2004 32:1792-1797.
#' @references Edgar RC. BMC Bioinformatics 5(1):113.
#' @seealso custom.aln(), list.hom(), parse.hssp(), get.hssp(), shannon()
#' @importFrom bio3d seqbind
#' @importFrom bio3d seqaln
#' @importFrom bio3d write.fasta
#' @export

msa <- function(sequences, ids = names(sequences), sfile = FALSE, inhouse = FALSE){

  ## --- Checking the inputs
  if (length(sequences) < 2){
    stop("At least two sequences are required!")
  } else if (length(sequences) != length(ids)){
    stop("The number of sequences and sequences' ids doesn't match!")
  }

  ## --- Using the program MUSCLE, which must be in the user system path
  if (inhouse){
    seqs <- lapply(sequences, function(x) strsplit(x, split = "")[[1]])
    sqs <- bio3d::seqbind(seqs[[1]], seqs[[2]], blank = "-")
    c <- 2
    while (c < length(sequences)){
      c <- c + 1
      sqs <- bio3d::seqbind(sqs, seqs[[c]], blank = "-")
    }
    aln <- bio3d::seqaln(sqs, id = ids, exefile = "muscle")
    aln$seq <- sequences

    if (sfile != FALSE){
      bio3d::write.fasta(aln, file = sfile)
    }
    if (file.exists("aln.fa")){
      system("rm aln.fa")
    }

    return(aln)

  } else {
    ## --- Using the Biostring and muscle R packages
    # seq <- Biostrings::AAStringSet(sequences)

    if (requireNamespace('Biostrings', quietly = TRUE)){
      seq <- Biostrings::AAStringSet(sequences)
    } else {
      stop("You must install the package Biostrings in order to use this function")
    }

    # aln1 <- muscle::muscle(seq)

    if (requireNamespace('muscle', quietly = TRUE)){
      aln1 <- muscle::muscle(seq)
    } else {
      stop("You must install the package muscle in order to use this function")
    }

    aln <- list()
    aln$seq <- sequences
    aln$ids <- ids
    aln$aln <- as.character(aln1)
    l <- sapply(aln$aln, function (x) strsplit(x, split = ""))
    aln$ali <- matrix(unlist(l), nrow = length(sequences), byrow = TRUE)

    if (sfile != FALSE){
      for (i in 1:length(aln$aln)){
        t <- paste(">", aln$ids[i], sep = "")
        cat(t, file = sfile, append = TRUE)
        tt <- paste("\n", aln$aln[i], "\n", sep = "" )
        cat(tt, file = sfile, append = TRUE)
      }
    }
    return(aln)
  }
}

