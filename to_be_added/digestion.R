# TRYPSIN:
# trypsin predominantly cleaves proteins at the C-terminal side of the amino acids lysine and arginine
# except when either is bound to a C-terminal proline (cut after arginine and lysine, but not before proline), although large-scale mass spectrometry data
# suggest cleavage occurs even with proline (10.1021/pr0705035).
# Missed cleavage sites are frequently observed when acidic amino acids, aspartic and glutamic acids, are present near the cleavage site. 
# Also, the sequence motifs with successive lysine and/or arginine residues represent a source of missed cleaved sites.
# Anal Chem. 2015 Aug 4;87(15):7636-43. doi: 10.1021/acs.analchem.5b00866. Epub 2015 Jul 21.

# Lys-C:
# Endoproteinase Lys-C is a protease that cleaves proteins on the C-terminal side of lysine residues.
# This enzyme is naturally found in the bacterium Lysobacter enzymogenes and is commonly used in
# protein sequencing.Lys-C activity is optimal in the pH range 7.0 - 9.0.

# https://web.expasy.org/peptide_cutter/

## ------- distances.R ------------ ##
#                                    #
#   digest                           #
#                                    #
## -------------------------------- ##

## ------------------------------------------------------------------------------ ##
# digest <- function(target, protease = 'trypsin', kr = TRUE, uniprot = TRUE, ...) #
## ------------------------------------------------------------------------------ ##
#' Find Cleveage Patterns
#' @description Predicts potential cleavage sites cleaved by proteases or chemicals in a given protein sequence
#' @usage function(target, protease = 'trypsin', kr = TRUE, uniprot = TRUE, partial = FALSE)
#' @param target a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param protease a character string among ('trypsin', 'lys-c') indicating the protease to use.
#' @param kr stands for Keil rules. When set to FALSE trypsin is allowed to cut before proline.
#' @details Bla, bla, bla.
#' @return Bla, bla
#' @examples
#' @seealso
#' @export

digest <- function(target, protease = 'trypsin', kr = TRUE, uniprot = TRUE, partial = FALSE){

  if (uniprot){
    seq <- ptm::get.seq(id = target) # Protein sequence
    id <- target
  } else {
    seq <- target
    id <- 'user-provided sequence'
  }

  if (protease == 'trypsin' & kr == FALSE){
    at <- gregexpr("[K|R]", seq)[[1]] # cut regardless a Pro is or not found towards CT
    if (partial == FALSE){
      if (substr(seq, nchar(seq), nchar(seq)) %in% c('K', 'R')){
        at <- at[-length(at)]
      }
      output <- as.data.frame(matrix(rep(NA, (1+length(at))*4), ncol = 4))
      names(output) <- c('start', 'end', 'peptide', 'mass')
      output$start <- c(1, (at + 1))
      output$end <- c(at, nchar(seq))
      for (i in 1:nrow(output)){
        output$peptide[i] <- substr(seq, output$start[i], output$end[i])
      }
    } else {
      # Digestiones parciales
    }
  } else if (protease == 'trypsin' & kr == TRUE){
    at <- ggregexpr("[K|R](?!P)", seq, perl = TRUE)[[1]]
  } else if (protease == 'lys-c'){

  }
  return(output)
}

