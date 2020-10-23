## ------- pdb_ancillary.R ------ ##
#                                  #
#   pdb.quaternary                 #
#   pdb.chain                      #
#   pdb2uniprot                    #
#   uniprot2pdb                    #
#   pdb.res                        #
#   pdb.pep                        #
#   pdb.select                     #
#                                  #
## ------------------------------ ##


## ---------------------------------------------------------------- ##
#      pdb.quaternary <- function(pdb, keepfiles = FALSE)                  #
## ---------------------------------------------------------------- ##
#' Protein Subunit Composition
#' @description Determines the subunit composition of a given protein.
#' @usage pdb.quaternary(pdb, keepfiles = FALSE)
#' @param pdb the path to the PDB of interest or a 4-letter identifier.
#' @param keepfiles logical, if TRUE the fasta file containing the alignment of the subunits is saved in the current directory, as well as the splitted pdb files.
#' @details A fasta file containing the alignment among the subunit sequences can be saved in the current directory if required.
#' @return This function returns a list with four elements: (i) a distances matrix, (ii) the sequences, (iii) chains id, (iv) the PDB ID used.
#' @author Juan Carlos Aledo
#' @examples pdb.quaternary('1bpl')
#' @export

pdb.quaternary <- function(pdb, keepfiles = FALSE){

  t <- strsplit(pdb, split = "\\/")[[1]]
  pdb_id <- substring(t[length(t)], 1,4)
  chains <- suppressWarnings(pdb.chain(pdb))

  sequences <-  character(length(chains))
  oldw <- getOption("warn")
  options(warn = -1) # avoids unnecessary warnings: 'pdb exists. Skipping download'
  for (i in seq_len(length(chains))){
    t <- paste(pdb_id, ":", chains[i], sep = "")
    sequences[i] <-  get.seq(t, db = 'pdb')
  }
  options(warn = oldw) # restores the warnings

  outfile <- paste(pdb_id, ".fa", sep = "")

  if (length(chains) == 1){ # The protein is a single monomer
    output <- list()
    output[[1]] <- NA
    output[[2]] <- sequences
    output[[3]] <- chains
    output[[4]] <- pdb_id
  } else {
    oldw <- getOption("warn")
    options(warn = -1) # avoids unnecessary warnings: 'pdb exists. Skipping download'
    myaln <- msa(sequences, ids = chains, sfile = outfile)
    if (requireNamespace('seqinr', quietly = TRUE)){
      aln <- seqinr::read.alignment(file = outfile, format = 'fasta' )
      d <- seqinr::dist.alignment(aln, "identity")
    } else {
      stop("The package seqinr must be installed to use this function")
    }
    options(warn = oldw) # restores the warnings
    output <- list()
    output[[1]] <- as.matrix(d)
    output[[2]] <- sequences
    output[[3]] <- chains
    output[[4]] <- pdb_id
  }

  if (!keepfiles){
    file.remove(outfile)
    if (file.exists(paste(pdb, ".pdb", sep = ""))){
      file.remove(paste(pdb, ".pdb", sep = ""))
    }
    unlink('split_chain', recursive = TRUE)
  }

  return(output)
}


## ---------------------------------------------------------------- ##
#           pdb.chain <- function(pdb, keepfiles = FALSE)            #
## ---------------------------------------------------------------- ##
#' Download and/or Split PDB Files.
#' @description Downloads a PDB file (if requiered) and splits it to provide a file by chain.
#' @usage pdb.chain(pdb, keepfiles = FALSE)
#' @param pdb the path to the PDB of interest or a 4-letter identifier.
#' @param keepfiles logical, if TRUE the function makes a 'temp' directory within the current directory and save in it a pdb file for each chain present in the given structure.
#' @return The function returns a chr vector where each coordinate is a chain from the structure.
#' @author Juan Carlos Aledo
#' @examples pdb.chain('1bpl')
#' @importFrom stats na.omit
#' @importFrom bio3d read.pdb
#' @importFrom bio3d get.pdb
#' @importFrom bio3d pdbsplit
#' @export

pdb.chain <- function(pdb, keepfiles = FALSE){

  if (nchar(pdb) < 4){
    stop("A proper pdb argument must be given!")
  }
  t <- bio3d::read.pdb(pdb)
  chains <- unique(names(t$seqres))

  if (keepfiles == TRUE & nchar(pdb) == 4){ # the input is a PDB ID
    bio3d::get.pdb(pdb, split = TRUE)
  } else if (keepfiles == TRUE & nchar(pdb) > 4){
    bio3d::pdbsplit(pdb)
  }
  return(chains)
}


## ---------------------------------------------------------------- ##
#          pdb2uniprot <- function(pdb, chain)              #
## ---------------------------------------------------------------- ##
#' Returns the UniProt ID Given the PDB and Chain IDs
#' @description Return the uniprot id of a given chain within a PDB structure.
#' @usage pdb2uniprot(pdb, chain)
#' @param pdb the 4-letter PDB identifier.
#' @param chain letter identifying the chain.
#' @return The function returns the UniProt ID for the chain of interest.
#' @author Juan Carlos Aledo
#' @examples pdb2uniprot('1u8f', 'O')
#' @seealso uniprot2pdb(), id.mapping()
#' @export

pdb2uniprot <- function(pdb, chain){

  pdb <- tolower(pdb)

  d <- NULL
  baseUrl <- "https://github.com/jcaledo/pdb_chain_uniprot/blob/master/"
  call <- paste(baseUrl, pdb, ".Rda?raw=true",sep = "")
  resp <- try(load(url(call)), silent = TRUE)
  if (inherits(resp, "try-error")) {
    text <- "Sorry, no results were found for this entry"
    return(text)
  }

  t <- d[which(d$PDB == pdb & d$CHAIN == chain),]
  t <- data.frame(lapply(t, as.character), stringsAsFactors = FALSE)
  output <- as.character(t$SP_PRIMARY)
  attr(output, 'RES_BEG') <- t$RES_BEG
  attr(output, 'RES_END') <- t$RES_END
  attr(output, 'PDB_BEG') <- t$PDB_BEG
  attr(output, 'PDB_END') <- t$PDB_END
  return(output)
}
## ---------------------------------------------------------------- ##
#                  uniprot2pdb <- function(up_id)                    #
## ---------------------------------------------------------------- ##
#' Return the PDB and Chain IDs of the Provided UniProt Protein
#' @description Returns the PDB and chain IDs of the provided protein.
#' @usage uniprot2pdb(up_id)
#' @param up_id the UniProt ID.
#' @return The function returns a dataframe with info about the PDB related to the protein of interest.
#' @author Juan Carlos Aledo
#' @examples uniprot2pdb("P04406")
#' @seealso pdb2uniprot(), id.mapping()
#' @export

uniprot2pdb <- function(up_id){

  d <- NULL
  baseUrl <- "https://github.com/jcaledo/uniprot_pdb_chain/blob/master/"
  call <- paste(baseUrl, up_id, '.Rda?raw=true', sep = "")
  resp <- try(load(url(call)), silent = TRUE)
  if (inherits(resp, "try-error")) {
    text <- "Sorry, no results were found for this entry"
    return(text)
  }
  output <- d[which(d$SP_PRIMARY == up_id),]
  return(output)
}

## ---------------------------------------------------------------- ##
#             pdb.res <- function(at, up, pdb, chain)                #
## ---------------------------------------------------------------- ##
#' Check Whether a Given Residue is Found in the PDB
#' @description Checks whether or not a given residue is resolved in the PDB structure.
#' @usage pdb.res(at, up, pdb, chain)
#' @param at the position in the primary structure of the residue (according to the UniProt sequence).
#' @param up the UniProt ID.
#' @param pdb the 4-letter PDB identifier.
#' @param chain letter identifying the chain.
#' @details This function checks if a given residue in the Uniprot sequence is found in the PDB.
#' @return The functions returns TRUE if the residue is found in the PDB sequence (and gives the position in that sequence). If the residue of interest is not found in the PDB the function returns FALSE.
#' @author Juan Carlos Aledo
#' @examples pdb.res(at = 361, up = 'P48163', pdb = '2aw5', chain = 'A')
#' @seealso pdb.pep()
#' @export
pdb.res <- function(at, up, pdb, chain){
  a <- renum.pdb(pdb, chain, up)
  if(a$pdb[which(a$uni_pos == at)] == aa.at(at, up)){
    output <- TRUE
    attr(output, 'pdb_pos') <- a$pdb_pos[which(a$uni_pos == at)]
  } else {
    output <- FALSE
  }

  return(output)
}

## ---------------------------------------------------------------- ##
#                 pdb.pep <- function(pep, pdb)                      #
## ---------------------------------------------------------------- ##
#' Check Whether an Oligopeptide is Found in the PDB
#' @description Checks whether a given oligopeptide is resolved in the PDB.
#' @usage pdb.pep(pep, pdb)
#' @param pep character string corresponding to the sequence of the oligopeptide.
#' @param pdb the 4-letter PDB identifier.
#' @details The oligopeptide sequence must be given in one letter amino acid code.
#' @return The functions returns TRUE if peptide is found in the PDB sequence (and gives the starting position), and FALSE otherwise.
#' @author Juan Carlos Aledo
#' @examples pdb.pep(pep = 'IVKGRASLTQEQ' , pdb = '2aw5')
#' @seealso pdb.res()
#' @importFrom bio3d read.pdb
#' @importFrom bio3d aa321
#' @export
pdb.pep <- function(pep, pdb){
  pep <- toupper(pep)

  ## ----- Get de pdb sequence
  oldw <- getOption("warn")
  options(warn = -1) # avoids unnecessary warnings: 'pdb exists. Skipping download'
  t <- bio3d::read.pdb(pdb)$atom
  options(warn = oldw) # restores the warnings
  t <- t[which(t$elety == 'CA' & t$type == 'ATOM'), ]
  seq <- paste(bio3d::aa321(t$resid), collapse = "") # sequence resolved in PDB


  ## ----- Search the peptide in the sequence
  if (gregexpr(pep, seq)[[1]][1] != -1){
    output <- TRUE
    attr(output, 'at') <- gregexpr(pep, seq)[[1]]
    attr(output, 'chain') <- t$chain[gregexpr(pep, seq)[[1]]]
  } else {
    output <- FALSE
  }

  return(output)
}

## ---------------------------------------------------------------- ##
#           pdb.select <- function(up_id, threshold = 0.9)           #
## ---------------------------------------------------------------- ##
#' Select the PDB with the Optimal Coverage to the UniProt Sequence
#' @description Select the PDB and chain with the optimal coverage to a given UniProt sequence.
#' @usage pdb.select(up_id, threshold = 0.9)
#' @param up_id the UniProt ID.
#' @param threshold coverage value that when reached the search is halted.
#' @return A list of two elements: (i) the PDB ID and (ii) the chain. The coverage with the UniProt sequence is given as an attribute.
#' @author Juan Carlos Aledo
#' @examples pdb.select('P01009', threshold = 0.8)
#' @seealso pdb.quaternary(), pdb.chain(), pdb.res(), pdb.pep(), uniprot2pdb(), pdb2uniprot()
#' @export
pdb.select <- function(up_id, threshold = 0.9){

  pdbs <- uniprot2pdb(up_id)
  if (grepl('Sorry', pdbs[1])){
    return("NO PDB FOUND")
  } else { ## --------- Choose the PDB and chain when they exist
    oldw <- getOption("warn") # To avoid some irrelevant warnings
    options(warn = -1)
    identity <- rep(NA, nrow(pdbs))
    match <- 0 # fraction of the uniprot sequence contained in the pdb sequence
    c <- 0
    while (match < threshold & c < nrow(pdbs)){ # ---- Stop search if 90 % uniprot is found in pdb
      c <- c + 1
      t <- renum.pdb(pdb = as.character(pdbs$PDB[c]),
                     chain = as.character(pdbs$CHAIN[c]),
                     uniprot = up_id)
      match <- sum(t$uniprot == t$pdb)/nrow(t)
      identity[c] <- match
    }
    options(warn = oldw) # Restores warnings

    ind <- which(identity == max(identity, na.rm = TRUE))
    pdb <-  as.character(pdbs$PDB[ind])
    chain <- as.character(pdbs$CHAIN[ind])
    coverage <- identity[ind]

    output <- list(pdb, chain)
    attr(output, 'coverage') <- round(coverage, 3)
    return(output)
  }
}
