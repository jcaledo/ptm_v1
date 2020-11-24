## ---------- series.R ------------ ##
#                                    #
#     series.01                      #
#     series.property                #
#                                    #
#                                    #
## -------------------------------- ##


## ---------------------------------------------------------------------- ##
#     series.01(seq, residues = 'M', uniprot = TRUE, plot = TRUE, ...)          #
## ---------------------------------------------------------------------- ##
#' Convert Amino Acid Sequences into a 01 Absence/Presence Serie
#' @description The amino acid sequence is translated into a numerical (0 and 1) sequence indicating the absence (0) or presence (1) of the indicated amino acid(s).
#' @usage series.01(seq, residues = 'M', uniprot = TRUE, plot = TRUE)
#' @param seq a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param residues a character vector specifying the amino acids to be analyzed.
#' @param uniprot logical, if TRUE the argument 'seq' shoud be an ID.
#' @param plot logical. If TRUE (the default) the serie is plotted.
#' @details
#' @return Returns a numeric vector with the values of the serie.
#' @author Juan Carlos Aledo
#' @examples series.01(seq = 'P0DP23', residues = c('M', 'Q'))
#' @seealso series.property(), plot.seq()
#' @export

series.01 <- function(seq, residues = 'M', uniprot = TRUE, plot = TRUE, ...){

  ## ------- The sequence to be converted ------ ##
  if (uniprot){
    seq <- ptm::get.seq(seq, as.string = FALSE)[1,]
  }
  if (length(seq) == 1){
    seq <- strsplit(seq, split = "")[[1]]
  }

  ## -------------- The translation ------------ ##
  series <- rep(0, length(seq))
  series[which(seq %in% residues)] <- 1
  attr(series, 'residues') <- residues

  ## ---------------- Plot Series -------------- ##
  z <- list(...)
  if (plot){
    plot(1:length(series), series, ylab = "Presence/Absence",
         xlab = "Residue Position", typ = 'h', main = z$title)
  }

  return(series)

}

## ---------------------------------------------------------------------- ##
#   series.property(seq, property = 'eiip', uniprot = TRUE, plot = TRUE)   #
## ---------------------------------------------------------------------- ##
#' Convert Amino Acid Sequences into Numerical Serie
#' @description Translate an amino acid sequence into a numerical serie
#' @usage series.property(seq, property = 'eiip', uniprot = TRUE, plot = TRUE)
#' @param seq a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param property a character string specifying the amino acid physicochemical property of interest (see details).
#' @param uniprot logical, if TRUE the argument 'seq' shoud be an ID.
#' @param plot ogical. If TRUE (the default) the serie is plotted.
#' @details The property chosen should be one among: "uniprot_2019" (relative frequencies in database), "human" (relative frequencies in human proteins), "eiip" (Electro-Ion Interaction Potentials),  "volume" (side chain volume), "polarizability", "av_hyd" (averaged hydrophobicity), "fasman_hyd" (Fasman's hydrophobicity), "pi_hel" (pi-helix propensity),  "hel" (helix propensity), "a_hel" (alpha-helix propensity), "b_sheet1" (beta-sheet propensity),  "b_sheet2" (beta-sheet propensity), "B_factor". Nevertheless, there are more than half a thousand of amino acid indices that have been tabulated and published, so the user can provide his/her own index as a named vector of dimension 20 with the following order: "A", "R","N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V".
#' @return Returns a numeric vector with the values of the chosed property along the primary structure.
#' @author Juan Carlos Aledo
#' @examples series.property(seq = 'P01009', property = 'eiip')
#' series.property(seq = 'P0DP23', property = myownindex)
#' @references https://www.genome.jp/aaindex
#' @seealso series.01(), plot.seq()
#' @export

series.property <- function(seq, property = 'eiip', uniprot = TRUE, plot = TRUE){


  ## ------- The sequence to be converted ------ ##
  if (uniprot){
    seq <- ptm::get.seq(seq, as.string = FALSE)[1,]
  }
  if (length(seq) == 1){
      seq <- strsplit(seq, split = "")[[1]]
  }


  ## ---------- The property to be used -------- ##
  if (length(property == 20)){
    p <- property
  } else {
    aai <- ptm:::aai
    p <- aai[,names(aai) == property]
    names(p) <- aai$aa
  }

  ## -------------- The translation ------------ ##
  series <- sapply(seq, function(x) p[which(names(p) == x)])
  names(series) <- seq
  attr(series, 'property') <- property

  ## ---------------- Plot Series -------------- ##
  if (plot){
    plot(1:length(seq), series, ylab = "EIPP", xlab = "Residue Position", typ = 'b')
    points(which(seq == "M"), series[which(names(series) == "M")], pch = 19, col = 'red')
  }

  return(series)
}


## ---------------------------------------------------------------------- ##
#   series.property(seq, property = 'eiip', uniprot = TRUE, plot = TRUE)   #
## ---------------------------------------------------------------------- ##
#' Convert Amino Acid Sequences into Numerical Serie
#' @description Translate an amino acid sequence into a numerical serie
#' @usage series.property(seq, property = 'eiip', uniprot = TRUE, plot = TRUE)
#' @param seq a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param property a character string specifying the amino acid physicochemical property of interest (see details).
#' @param uniprot logical, if TRUE the argument 'seq' shoud be an ID.
#' @param plot ogical. If TRUE (the default) the serie is plotted.
#' @details The property chosen should be one among .... Nevertheless, there are more than half a thousand of amino acid indices that have been tabulated and published (https://www.genome.jp/aaindex), so the user can provide his/her own index a named vector of dimension 20 with the following order: "A", "R","N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V".
#' @return Returns a nume vector with the values of the chosed property along the primary structure.
#' @author Juan Carlos Aledo
#' @examples series.property(seq = 'P01009', property = 'eiip')
#' series.property(seq = 'P0DP23', property = myownindex)
#' @seealso
#' @importFrom bio3d read.fasta
#' @export



