## ------- pdb_ancillary.R ------ ##
#                                  #
#   pdb.seq                        #
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
#             pdb.seq <- function(pdb)                               #
## ---------------------------------------------------------------- ##
#' Get Chain Sequences
#' @description Gets the sequences of the chain find in a given PDB.
#' @usage pdb.seq(pdb)
#' @param pdb the 4-letter PDB identifier.
#' @return  Returns a dataframe with as many rows as different chains are present in the PDB. For each row six variables are returned: (i) the entry id, (ii) the entity id, (iii) the chain, (iv) the protein name, (v) the species and (vi) the sequence.
#' @author Juan Carlos Aledo
#' @examples pdb.seq('1bpl')
#' @importFrom httr GET
#' @importFrom httr content
#' @export

pdb.seq <- function(pdb){

  ## ----------------- Check PDB argument ------------ #
  if (nchar(pdb) != 4){
    stop( "Please, provide a proper PDB ID")
  } else {
    call <- paste('https://www.rcsb.org/fasta/entry/',
                  pdb, '/download', sep = "")
  }

  ## -------- Client <-> Server Communication -------- #
  if (!is.null(call)){
    resp <- try(httr::GET(call), silent = FALSE)
    if (inherits(resp, "try-error")) {
      text <- "LOST CONNECTION"
    } else {
      text <- httr::content(resp, "text", encoding = "utf-8")
    }
    retry = 0
    while ((grepl("^LOST CONNECTION", text) || httr::http_error(resp) ||
            grepl("^ERROR", text) || grepl("Nothing has been found",
                                           text)) && retry < 3) {
      retry = retry + 1
      Sys.sleep(5)
      resp <- try(httr::GET(url), silent = TRUE)
      if (inherits(resp, "try-error")) {
        text <- "LOST CONNECTION"
      }
      else {
        text <- httr::content(resp, "text", encoding = "utf-8")
      }
    }
  }
  ## ----------- Parsing the response ------------- #
  t <- strsplit(text, split = ">")[[1]][-1]
  seq <- data.frame(entry = rep(NA, length(t)),
                    entity = rep(NA, length(t)),
                    chain = rep(NA, length(t)),
                    name = rep(NA, length(t)),
                    species = rep(NA, length(t)),
                    sequence = rep(NA, length(t)))

  for (i in 1:nrow(seq)){
    tt <- strsplit(t[i], split = "\n")[[1]]
    z <- strsplit(tt[1], "\\|")[[1]]
    seq$entry[i] <- strsplit(z[1], split = "_")[[1]][1]
    seq$entity[i] <- strsplit(z[1], split = "_")[[1]][2]
    chains <- strsplit(z[2], split = " ")[[1]][2]

    # chains <- trimws(gsub('[Chains|Chain]', "", z[2]))

    seq$chain[i] <- trimws(chains)
    seq$name[i] <- z[3]
    seq$species[i] <- z[4]
    seq$sequence[i] <- tt[2]
  }
  ## ------------------ Output ------------------- #
  return(seq)
}

## ---------------------------------------------------------------- ##
#      pdb.quaternary <- function(pdb, keepfiles = FALSE)            #
## ---------------------------------------------------------------- ##
#' Protein Subunit Composition
#' @description Determines the subunit composition of a given protein.
#' @usage pdb.quaternary(pdb, keepfiles = FALSE)
#' @param pdb the path to the PDB of interest or a 4-letter identifier.
#' @param keepfiles logical, if TRUE the fasta file containing the alignment of the subunits is saved in the current directory, as well as the split pdb files.
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

  for (i in seq_len(length(chains))){
    t <- paste(pdb_id, ":", chains[i], sep = "")
    sequences[i] <-  get.seq(t, db = 'pdb')
  }

  outfile <- paste(pdb_id, ".fa", sep = "")

  if (length(chains) == 1){ # The protein is a single monomer
    output <- list()
    output[[1]] <- NA
    output[[2]] <- sequences
    output[[3]] <- chains
    output[[4]] <- pdb_id
  } else {
    myaln <- msa(sequences, ids = chains, sfile = outfile)
    if (requireNamespace('seqinr', quietly = TRUE)){
      aln <- seqinr::read.alignment(file = outfile, format = 'fasta' )
      d <- seqinr::dist.alignment(aln, "identity")
    } else {
      stop("The package seqinr must be installed to use this function")
    }
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
#' @description Downloads a PDB file (if required) and splits it to provide a file by chain.
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
#          pdb2uniprot <- function(pdb, chain)                       #
## ---------------------------------------------------------------- ##
#' Return the UniProt ID Given the PDB and Chain IDs
#' @description Returns the uniprot id of a given chain within a PDB structure.
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
#' @examples \dontrun{pdb.res(at = 361, up = 'P48163', pdb = '2aw5', chain = 'A')}
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
#' @examples \dontrun{pdb.pep(pep = 'IVKGRASLTQEQ' , pdb = '2aw5')}
#' @seealso pdb.res()
#' @importFrom bio3d read.pdb
#' @importFrom bio3d aa321
#' @export
pdb.pep <- function(pep, pdb){
  pep <- toupper(pep)

  ## ----- Get de pdb sequence
  t <- suppressWarnings(bio3d::read.pdb(pdb)$atom) # avoids warning: 'pdb exists. Skipping download'
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
#' @examples \dontrun{pdb.select('P01009', threshold = 0.8)}
#' @seealso pdb.quaternary(), pdb.chain(), pdb.res(), pdb.pep(), uniprot2pdb(), pdb2uniprot()
#' @export
pdb.select <- function(up_id, threshold = 0.9){

  pdbs <- uniprot2pdb(up_id)
  if (grepl('Sorry', pdbs[1])){
    return("NO PDB FOUND")
  } else { ## --------- Choose the PDB and chain when they exist
    pdbs <- pdbs[which(nchar(pdbs$CHAIN)==1),]
    if (nrow(pdbs) == 0){
      return("NO PDB FOUND")
    }
    identity <- rep(NA, nrow(pdbs))
    match <- 0 # fraction of the uniprot sequence contained in the pdb sequence
    c <- 0
    while (match < threshold & c < nrow(pdbs)){ # ---- Stop search if 90 % uniprot is found in pdb
      c <- c + 1
      t <- suppressWarnings(renum.pdb(pdb = as.character(pdbs$PDB[c]),
                     chain = as.character(pdbs$CHAIN[c]),
                     uniprot = up_id))
      match <- sum(t$uniprot == t$pdb)/nrow(t)
      identity[c] <- match
    }

    ind <- which(identity == max(identity, na.rm = TRUE))
    pdb <-  as.character(pdbs$PDB[ind])
    chain <- as.character(pdbs$CHAIN[ind])
    coverage <- identity[ind]

    output <- list(pdb, chain)
    attr(output, 'coverage') <- round(coverage, 3)
    return(output)
  }
}
