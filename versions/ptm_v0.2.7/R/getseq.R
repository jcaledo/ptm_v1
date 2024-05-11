## ---------- getseq.R ------------ ##
#                                    #
#     get.seq                        #
#     pdb.seq                        #
#     pdb.quaternary                 #
#     pdb.chain                      #
#     pdb2uniprot                    #
#     uniprot2pdb                    #
#     prot2codon                     #
#     pdb.select                     #
#     id.mapping                     #
#     id.features                    #
#     species.mapping                #
#     species.kegg                   #
#     kegg.uniprot                   #
#     uniprot.kegg                   #
#     uniprot.pdb                    #
#     pdb.uniprot                    #
#                                    #
## -------------------------------- ##

## ---------------------------------------------------------------- ##
#      get.seq <- function(id, db = 'uniprot', as.string = TRUE)     #
## ---------------------------------------------------------------- ##
#' Import a Protein Sequence from a Database
#' @description Imports a protein sequence from a selected database.
#' @usage get.seq(id, db = 'uniprot', as.string = TRUE)
#' @param id the identifier of the protein of interest.
#' @param db a character string specifying the desired database; it must be one of 'uniprot', 'metosite', 'pdb', 'kegg-aa', 'kegg-nt'.
#' @param as.string logical, if TRUE the imported sequence will be returned as a character string.
#' @details MetOSite uses the same type of protein ID than UniProt. However, if the chosen database is PDB, the identifier should be the 4-character unique identifier characteristic of PDB, followed by colon and the chain of interest. For instance, '2OCC:B' means we are interested in the sequence of chain B from the structure 2OCC. KEGG used its own IDs (see examples).
#' @return Returns a protein (or nucleotide) sequence either as a character vector or a as a character string.
#' @author Juan Carlos Aledo
#' @examples get.seq('P01009')
#' @importFrom jsonlite fromJSON
#' @export

get.seq <- function(id, db = 'uniprot', as.string = TRUE){

  db <- tolower(db)

  if (db == 'uniprot'){
    baseUrl <- "http://uniprot.org/uniprot/"
    call <- paste(baseUrl, id, ".fasta", sep = "")

  } else if (db == 'metosite'){
    baseUrl <- 'https://metosite.uma.es/api/proteins/scan/'
    call <- paste(baseUrl, id, sep = "")

  } else if (db == 'pdb'){
    text <- 'no text'
    id <- strsplit(id, split=':')[[1]]
    df <- tryCatch(
      {
        pdb.seq(id[1])
      },
      error = function(cond){
        return(NULL)
      }
    )

    if (is.data.frame(df)){ # pdb.seq successed
      call <- NULL
      if (!is.na(id[2])){ # a given chain
        chains <- strsplit(gsub(",", "", paste(df$chain, collapse = "")), split = "")[[1]]
        if (id[2] %in% chains){
          seq <- df$sequence[grepl(id[2], df$chain)]
        } else if (grepl(id[2], chains)){
          seq <- df$sequence[grepl(id[2], df$chain)]
        } else {
          stop("The chosen chain is not found in the PDB")
        }
      } else { # no chain specified
        seq <- paste(df$sequence, collapse = "")
      }
    }

  } else if (db == 'kegg-aa' | db == 'kegg-nt'){
    t <- paste(strsplit(db, split = '-')[[1]][2], 'seq', sep = '')
    if (requireNamespace('KEGGREST', quietly = TRUE)){
      text <- tryCatch(
        {
          as.character(KEGGREST::keggGet(id, t))
        },
        error = function(cond){
          # message(cond)
          return(NULL)
        }
      )
    } else {
      message("Please, install the KEGGREST pakage to get this functionality")
      return(NULL)
    }
    call <- NULL

  } else{
    stop('You should indicate a proper DB')
  }

  ## -------- Client <-> Server Communication -------- ##
  if (!is.null(call)){
    text <- gracefully_fail(call)
  }

  ## -------------- Parsing the response --------------- ##
  if (!is.null(text)){
    if (db == 'uniprot'){
      seq <- strsplit(text, split = "\\n")[[1]][-1]
      seq <- paste(seq, collapse = "")
    } else if (db == 'metosite'){
      data <- jsonlite::fromJSON(text, flatten = TRUE)
      seq <- data$prot_seq
    } else if (db == 'kegg-aa' | db == 'kegg-nt'){
      seq <- text
    }

    if (is.null(seq) || is.na(seq)){
      message(paste("The entry", id, "is not found in the", toupper(db), "database"))
      output <- NULL
    } else {
      output <- seq
      if (!as.string){
        output <- strsplit(seq, split ="")
      }
    }

  } else {
    message("Sorry, no result could be retrieved")
    output <- NULL
  }

  if(!is.null(output)){
    attr(output, "ID") <- id
    attr(output, "DB") <- db
  }

  return(output)
}


## ---------------------------------------------------------------- ##
#             pdb.seq <- function(pdb)                               #
## ---------------------------------------------------------------- ##
#' Get Chain Sequences
#' @description Gets the sequences of the chain find in a given PDB.
#' @usage pdb.seq(pdb)
#' @param pdb the 4-letter PDB identifier.
#' @return  Returns a dataframe with as many rows as different chains are present in the PDB. For each row six variables are returned: (i) the entry id, (ii) the entity id, (iii) the chain, (iv) the protein name, (v) the species and (vi) the sequence.
#' @author Juan Carlos Aledo
#' @examples \dontrun{pdb.seq('1bpl')}
#' @export

pdb.seq <- function(pdb){

  ## ----------------- Check PDB argument ------------ #
  if (nchar(pdb) != 4){
    stop( "Please, provide a proper PDB ID")
  } else {
    call <- paste('https://www.rcsb.org/fasta/entry/',
                  pdb, '/download', sep = "")
  }

  text <- gracefully_fail(call)
  if (is.null(text)){
    message("pdb.seq failed because it could not communicate with the rcsb server")
    return(NULL)
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
    x <- gsub("Chain", "", z[2])
    x <- gsub("[a-z]", "", x)
    x <- gsub("],", ",", x)
    x <- gsub("\\[", ",", x)
    x <- gsub("]", "", x)
    chains <- trimws(x)

    # chains <- strsplit(z[2], split = " ")[[1]][2]
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
#' @examples \dontrun{pdb.quaternary('1bpl')}
#' @importFrom seqinr read.alignment
#' @export

pdb.quaternary <- function(pdb, keepfiles = FALSE){

  t <- strsplit(pdb, split = "\\/")[[1]]
  pdb_id <- substring(t[length(t)], 1,4)

  chains <- suppressWarnings(pdb.chain(pdb_id))
  if (is.null(chains)){
    message("pdb.chain failed")
    return(NULL)
  }
  sequences <-  character(length(chains))

  for (i in seq_len(length(chains))){
    t <- paste(pdb_id, ":", chains[i], sep = "")
    sequences[i] <-  tryCatch(
      {
        get.seq(t, db = 'pdb')
      },
      error = function(cond){
        return(NA)
      })
  }

  outfile <- paste(pdb_id, ".fa", sep = "")

  if (length(chains) == 1){ # The protein is a single monomer
    output <- list()
    output[[1]] <- NA
    output[[2]] <- sequences
    output[[3]] <- chains
    output[[4]] <- pdb_id
  } else {
    myaln <- msa(sequences, ids = chains, sfile = paste("./", outfile, sep = "")) # <------- CHANGE TO TEMP DIR
    aln <- seqinr::read.alignment(file = outfile, format = 'fasta') # <------- CHANGE TO TEMP DIR
    d <- seqinr::dist.alignment(aln, "identity")
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
#' @examples \dontrun{pdb.chain('1bpl')}
#' @importFrom stats na.omit
#' @importFrom bio3d read.pdb
#' @importFrom bio3d get.pdb
#' @importFrom bio3d pdbsplit
#' @export

pdb.chain <- function(pdb, keepfiles = FALSE){

  if (nchar(pdb) < 4){
    stop("A proper pdb argument must be given!")
  }

  t <- tryCatch(
    {
      suppressWarnings(bio3d::read.pdb(pdb, verbose = FALSE))
    },
    error = function(cond){
      return(-999)
    }
  )

  if (is.numeric(t)){
    message("Sorry, read.pdb failed")
    return(NULL)
  }

  chains <- unique(names(t$seqres))

  if (keepfiles == TRUE & nchar(pdb) == 4){ # the input is a PDB ID
    tryCatch(
      {bio3d::get.pdb(pdb, split = TRUE)},
      error = function(cond){
        return(NULL)
      }
    )
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
#' @examples \dontrun{pdb2uniprot('1u8f', 'O')}
#' @seealso uniprot2pdb(), id.mapping()
#' @export

pdb2uniprot <- function(pdb, chain){

  pdb <- tolower(pdb)

  d <- NULL
  baseUrl <- "https://github.com/jcaledo/pdb_chain_uniprot/blob/master/"
  call <- paste(baseUrl, pdb, ".Rda?raw=true",sep = "")

  tryCatch(
    {
      load(url(call))
    },
    error = function(cond){
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )

  if (is.null(d)){
    message("Sorry, pdb_chain_uniprot couldn't be accessed")
    return(NULL)
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
#' @examples \dontrun{uniprot2pdb("P04406")}
#' @seealso pdb2uniprot(), id.mapping()
#' @export

uniprot2pdb <- function(up_id){

  d <- NULL
  baseUrl <- "https://github.com/jcaledo/uniprot_pdb_chain/blob/master/"
  call <- paste(baseUrl, up_id, '.Rda?raw=true', sep = "")
  tryCatch(
    {
      load(url(call))
    },
    error = function(cond){
      message("Sorry, pdb_chain_uniprot couldn't be accessed")
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )

  if (is.null(d)){
    message("Sorry, no result could be retrieve")
    return(NULL)
  }
  output <- d[which(d$SP_PRIMARY == up_id),]
  return(output)
}


## ---------------------------------------------------------------- ##
#       prot2codon <- function(prot, chain = "", laxity = TRUE)      #
## ---------------------------------------------------------------- ##
#' Find the Coding Triplets for a Given Protein
#' @description Finds the codons corresponding to a given protein.
#' @usage prot2codon(prot, chain = "", laxity = TRUE)
#' @param prot is either a UniProt or PDB id, or the path to a pdb file.
#' @param chain when prot corresponds to a pdb, the chain of interest must be provided.
#' @param laxity logical, if FALSE the program stop when a mismatch between the protein and the gene sequences is detected. Otherwise the program doesn't stop and at the end points out the mismatches.
#' @return Returns a dataframe with as many rows as residues has the protein.
#' @author Juan Carlos Aledo
#' @examples \dontrun{prot2codon('P01009')}
#' @importFrom bio3d read.pdb
#' @importFrom bio3d aa321
#' @importFrom seqinr translate
#' @importFrom seqinr s2c
#' @export

prot2codon <- function(prot, chain = "", laxity = TRUE){

  ## ---------- Primary structure of the protein ----------- ##
  if (regexpr("/",prot)[[1]] != -1){ # input is path to a local pdb file
    mypdb <- bio3d::read.pdb(prot)$atom
    seq <- bio3d::aa321(mypdb$resid[which(mypdb$elety == 'CA' &
                                            mypdb$chain == chain)])
    t <- strsplit(prot, split = "/")[[1]]
    t <- t[length(t)]
    t <- substring(t, 1,4)
    source <- 'pdb'

  } else if (nchar(prot) == 4){ # input is PDB ID
    id <- paste(prot, chain, sep = ":")
    seq <- get.seq(id, db = 'pdb', as.string = FALSE)
    if (is.null(seq)){
      message("Sorry, get.seq failed")
      return(NULL)
    } else {
      seq <- seq[[1]]
    }
    t <- prot
    source <- 'pdb'

  } else { # input should be a uniprot ID
    seq <- get.seq(prot, as.string = FALSE)
    if (is.null(seq)){
      message("Sorry, get.seq failed")
      return(NULL)
    } else {
      seq <- seq[[1]]
    }
    t <- prot
    source <- 'uniprot'
  }

  ## ------------------- Output structure -------------------- ##
  output <- as.data.frame(matrix(rep(NA, length(seq)*5), ncol = 5))
  names(output) <- c('id', 'chain', 'pos', 'aa', 'codon')
  output$id <- prot
  if (chain == ""){
    output$chain <- NA
  } else {
    output$chain <- chain
  }
  output$pos <- 1:length(seq)
  output$aa <- seq

  ## ----------------- Finding the codons ------------------ ##
  organism <- NULL
  organism <- species.mapping(t, db = source)
  if (is.null(organism)){
    message("Sorry, species.mapping failed")
    return(NULL)
  }

  if (source == 'pdb'){
    t <- pdb2uniprot(t, chain = chain)
    if (is.null(t)){
      message("Sorry, pdb2uniprot failed")
      return(NULL)
    }
  }

  kegg_id <- id.mapping(t, from = 'uniprot', to = 'kegg')
  if (is.null(kegg_id)){
    message("Sorry, we couldn't get a DNA seq from KEGG")
    return(NULL)
  } else {
    kegg_id <- kegg_id[1]
  }
  dna <- get.seq(kegg_id, 'kegg-nt')
  if (is.null(dna)){
    message("Sorry, get.seq failed")
    return(dna)
  }

  codon <- gsub("(.{3})", "\\1 ", dna)
  codon <- strsplit(codon, split = " ")[[1]]
  # Remove stop if present:
  codon <- codon[-which(codon %in% c("UAA", "UGA", "UAG", "TAA", "TGA", "TAG"))]

  translated <- seqinr::translate(seqinr::s2c(dna))

  sequences <- c(paste(seq, collapse = ""),
                 paste(translated, collapse = ""))
  aln <- msa(sequences, c("query", "translated"))

  k <- 1
  j <- 1
  for (i in 1:dim(aln$ali)[2]){
    if (aln$ali[1,i] == aln$ali[2,i]){
      output$codon[k] <- codon[j]
      k <- k + 1
      j <- j + 1
    } else if (aln$ali[1,i] == '-'){
      j <- j +1
    } else if (aln$ali[2,i] == '-'){
      output$codon[k] <- NA
      k <- k + 1
    } else { # mismatch
      output$codon[k] <- codon[j]
      k <- k + 1
      j <- j + 1
    }
  }

  ## ------------ Checking the translation -------------------- ##
  triplet <- c("GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG", "AAU", "AAC",
               "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG", "GGU", "GGC", "GGA", "GGG",
               "CAU", "CAC", "AUU", "AUC", "AUA", "UUA", "UUG", "CUU", "CUC", "CUA", "CUG", "AAA",
               "AAG", "AUG", "UUU", "UUC", "CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG",
               "AGU", "AGC", "ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA",
               "GUG", "UAA", "UGA", "UAG")
  aa <- c("A","A","A","A","R","R","R","R","R","R","N","N","D","D","C","C","Q","Q","E","E",
          "G","G","G","G","H","H","I","I","I","L","L","L","L","L","L","K","K","M","F","F",
          "P","P","P","P","S","S","S","S","S","S","T","T","T","T","W","Y","Y","V","V","V",
          "V","Stop","Stop","Stop")
  names(triplet) <- aa
  triplet <- gsub("U", "T", triplet)

  output$check <- NA
  for (i in 1:nrow(output)){
    if (!is.na(output$codon[i])){
      t <- as.character(output$codon[i])
      if (!laxity){
        if (output$aa[i] != names(triplet)[which(triplet == t)]){
          stop(paste("Problem at ", i))
        }
      }
      output$check[i] <- (output$aa[i] == names(triplet)[which(triplet == t)])
    }

  }

  dif <- nrow(output) - sum(output$check)
  message(paste(dif, "  mismatches found"))

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

  pdbs <- tryCatch(
    {
      uniprot2pdb(up_id)
    },
    error = function(cond){
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )

  if (is.null(pdbs)){
    message("NO PDB FOUND")
    return(NULL)
  } else { ## --------- Choose the PDB and chain when they exist
    pdbs <- pdbs[which(nchar(pdbs$CHAIN)==1),]
    if (nrow(pdbs) == 0){
      message("NO PDB FOUND")
      return(NULL)
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


## ---------------------------------------------------------------- ##
#               id.mapping <- function(id, from, to)                 #
## ---------------------------------------------------------------- ##
#' Identifier Mapping
#' @description Mapping between protein identifiers.
#' @usage id.mapping(id, from, to)
#' @param id the identifier to be converted.
#' @param from the type for the identifier of origin; it must be one of 'uniprot', 'pdb', or 'kegg'.
#' @param to the type for the identifier of destination; it must be one of 'uniprot', 'pdb', or 'kegg'.
#' @return Returns a character string corresponding to the requested identifier.
#' @author Juan Carlos Aledo
#' @examples \dontrun{id.mapping('P01009', from = 'uniprot', to = 'pdb')}
#' @seealso pdb2uniprot(), uniprot2pdb()
#' @export

id.mapping <- function(id, from, to){
  output <- NULL
  if (from == 'uniprot' & to == 'pdb'){
    output <- uniprot.pdb(id)
    if (is.null(output)){
      message("Sorry, id.mapping failed. Please, try again")
      return(NULL)
    }

  } else if (from == 'uniprot' & to == 'kegg'){
    output <- uniprot.kegg(id)
    if (is.null(output)){
      message('Sorry, no KEGG ID could be found!')
      return(NULL)
    }

  } else if (from == 'kegg' & to == 'uniprot'){
    output <- kegg.uniprot(id)
    if (is.null(output)){
      message('Sorry, no result could be retrieved')
      return(NULL)
    }

  } else if (from == 'kegg' & to == 'pdb'){
    output <- kegg.uniprot(id)
    if (is.null(output)){
      message('Sorry, no result could be retrieved')
      return(NULL)
    } else {
      output <- uniprot.pdb(output[1])
    }

  } else if (from == 'pdb' & to == 'uniprot'){
    output <- pdb.uniprot(id)
    if (is.null(output)){
      message('Sorry, no result could be retrieved')
      return(NULL)
    }

  } else if (from == 'pdb' & to == 'kegg'){
    output <- pdb.uniprot(id)
    if (is.null(output)){
      message('Sorry, no result could be retrieved')
      return(NULL)
    } else {
      output <- lapply(output, uniprot.kegg)
      if (identical(output, character(0))){
        message('Sorry, no KEGG ID could be found!')
        return(NULL)
      }
    }
  }

  if (!is.null(output)){
    attr(output, 'from') <- from
    attr(output, 'to') <- to
  }
  return(output)
}


## ---------------------------------------------------------------- ##
#           id.features <- function(id, features = "")               #
## ---------------------------------------------------------------- ##
#' Features Related to the Protein Entry
#' @description Obtains features related to the provided id.
#' @usage id.features(id, features = "")
#' @param id the UniProt identifier of the protein of interest.
#' @param features a string identifying the features (comma separated) to be recovered.
#' @details By default the the function provides info regarding the following features: id, reviewed, entry name and organism. If wished, this list of features can be expanded using the argument 'features'. There is a large list of features that can be retrieved. You can look up your relevant feature's name in the full list of UniProtKB found at https://www.uniprot.org/help/return_fields.
#' @return Returns a named list with the requested features.
#' @author Juan Carlos Aledo
#' @examples \dontrun{id.features('P04406', features = 'ec,keyword,xref_pdb')}
#' @seealso species.mapping()
#' @importFrom httr GET
#' @importFrom httr POST
#' @importFrom httr accept_json
#' @importFrom httr content
#' @importFrom httr add_headers
#' @importFrom utils read.table
#' @export

id.features <- function(id, features = "" ){

  # https://www.uniprot.org/help/id_mapping
  # https://www.uniprot.org/help/return_fields

  ## -------- Building the query ---------- ##
  col <- '?fields=accession%2Creviewed%2Cprotein_name%2Corganism_name'
  if (features != ""){
    features <- gsub(', ', ',', features)
    features <- gsub(',', "%2C", features)
    col <- paste(col, "%2C", features, sep = '')
  }

  isJobReady <- function(jobId) {
    pollingInterval = 5
    nTries = 20
    for (i in 1:nTries) {
      url <- paste("https://rest.uniprot.org/idmapping/status/", jobId, sep = "")
      r <- GET(url = url, accept_json())
      status <- content(r, as = "parsed")
      if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
        return(TRUE)
      }
      if (!is.null(status[["messages"]])) {
        print(status[["messages"]])
        return (FALSE)
      }
      Sys.sleep(pollingInterval)
    }
    return(FALSE)
  }

  getResultsURL <- function(redirectURL) {
    if (grepl("/idmapping/results/", redirectURL, fixed = TRUE)) {
      url <- gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
    } else {
      url <- gsub("/results/", "/results/stream/", redirectURL)
    }
  }

  files = list(
    ids = id,
    from = "UniProtKB_AC-ID",
    to = "UniProtKB"
  )

  my_headers <- httr::add_headers('User-Agent' = paste('R', 'metosite@uma.es'))
  r <- POST(url = "https://rest.uniprot.org/idmapping/run", body = files, my_headers, encode = "multipart", accept_json())
  submission <- content(r, as = "parsed")

  if (isJobReady(submission[["jobId"]])) {
    url <- paste("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]], sep = "")
    r <- GET(url = url, accept_json())
    details <- content(r, as = "parsed")
    url <- getResultsURL(details[["redirectURL"]])
    # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available

    url <- paste(url, col, "&format=tsv&size=500", sep = "")
    r <- GET(url = url, accept_json())
    output <- read.table(text = content(r), sep = "\t", header=TRUE)
  }

  if (!is.null(output)){
    return(output)

  } else {
    message("Uniprot server could not respond")
    output <- NULL
  }

  if (length(output) == 0){
    output <- NULL
  }
  return(output)
}


## ---------------------------------------------------------------- ##
#           species.mapping <- function(id, db = 'uniprot')          #
## ---------------------------------------------------------------- ##
#' Map Protein ID to Species
#' @description Maps a protein ID to its corresponding organism.
#' @usage species.mapping(id, db = 'uniprot')
#' @param id the identifier of the protein of interest.
#' @param db a character string specifying the corresponding database. Currently, only 'uniprot' or 'pdb' are valid options.
#' @return Returns a character string identifying the organism to which the given protein belong.
#' @author Juan Carlos Aledo
#' @examples \dontrun{species.mapping('P01009')}
#' @seealso id.features()
#' @export

species.mapping <- function(id, db = 'uniprot'){

  db <- tolower(db)
  if (!(db %in% c("uniprot", "pdb"))){
    stop("Option db should be one of uniprot or pdb")
  }
  if (db == 'pdb'){
    id <- id.mapping(id, from = 'pdb', to = 'uniprot')[1]
    if (is.null(id)){
      message("Sorry, id.mapping failed")
      return(NULL)
    }
  }

  url <- paste("https://www.uniprot.org/uniprot/", id, ".fasta", sep = "")
  text <- gracefully_fail(url)

  if (is.null(text)){
    message("Sorry, no result could be retrieved")
    return(NULL)
  } else {
    start <- gregexpr("OS=", text)[[1]]
    stop <- gregexpr("OX=", text)[[1]]
    species <- substr(text, start+3, stop-2)
  }

  if (nchar(species) == 0){
    message("Sorry, no result could be retrieved")
    species <- NULL
  }
  return(species)
}


## ---------------------------------------------------------------- ##
#     species.kegg <- function(organism, from = 'scientific')         #
## ---------------------------------------------------------------- ##
#' Convert Between Species Name and KEGG 3-Letter Code Format
#' @description Converts between species name and KEGG 3-letter code format.
#' @usage species.kegg(organism, from = 'scientific')
#' @param organism character string defining the organisms.
#' @param from string indicating the character of the provided name. It should be one of 'vulgar', 'scientific', '3-letter'.
#' @return Returns a dataframe with the entries matching the request.
#' @author Juan Carlos Aledo
#' @examples \dontrun{species.kegg('chempanzee', from = 'vulgar')}
#' @examples \dontrun{species.kegg('Pan paniscus')}
#' @examples \dontrun{species.kegg('ppo', from = '3-letter')}
#' @seealso id.features(), species.mapping()
#' @importFrom httr GET
#' @importFrom httr content
#' @export

species.kegg <- function(organism, from = 'scientific'){

  organisms <- NULL
  tryCatch(
    {
      load(url("https://github.com/jcaledo/kegg_species/blob/master/organisms.Rda?raw=true"))
    },
    error = function(cond){
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )

  if (is.null(organisms)){
    message("Sorry, no result could be retrieved")
    return(NULL)
  }

  organisms$organism <- as.character(organisms$organism)
  organisms$species <- as.character(organisms$species)

  if (from == 'vulgar'){
    ## ----- Try full match
    t <- paste("\\W", organism, "\\W", sep = "")
    output <- organisms[which(grepl(t, organisms$species)), ]
  } else if (from == 'scientific'){
    output <- organisms[which(regexpr(organism, organisms$species) != -1), ]
  } else if (from == '3-letter'){
    output <- organisms[which(organisms$organism == organism), ]
  } else {
    stop("A proper value for the parameter 'from' must be provided")
  }

  if (nrow(output) == 0 & from == 'vulgar'){
    organism <- tolower(organism)
    t <- paste("\\W", organism, "\\W", sep = "")
    output <- organisms[which(grepl(t, organisms$species)), ]
  }  else if (nrow(output) == 0 & from == 'scientific'){  # Try partial match using binomial name
    binomial <- strsplit(organism, " ")[[1]]
    binomial <- paste(binomial[1], binomial[2])
    output <- organisms[which(regexpr(binomial, organisms$species) != -1), ]
  }

  if (nrow(output) == 0 & from == 'scientific'){ # Try an even shorter partial match
    genus <- strsplit(organism, " ")[[1]][1]
    output <- organisms[which(regexpr(genus, organisms$species) != -1), ]
  }

  if (nrow(output) == 0){
      message('Sorry, could find a code for this organism')
      return(NULL)
  } else {
    return(output)
  }
}


## ---------------------------------------------------------------- ##
#                   kegg.uniprot <- function(id)                     #
## ---------------------------------------------------------------- ##
#' Identifier Mapping From KEGG to UniProt
#' @description Mapping between KEGG and UniProt protein identifiers.
#' @usage kegg.uniprot(id)
#' @param id the identifier to be converted.
#' @return Returns a character string corresponding to the requested identifier.
#' @author Juan Carlos Aledo
#' @examples \dontrun{kegg.uniprot('hsa:5265')}
#' @seealso id.mapping()
#' @export

kegg.uniprot <- function(id){

  if (requireNamespace("KEGGREST", quietly = TRUE)){
    output <- tryCatch(
      {
        KEGGREST::keggConv('uniprot', id)
      },
      error = function(cond){
        return(NULL)
      }
    )
  } else {
    message("The package KEGGREST must be installed to get this mapping")
    output <- NULL
  }

  if (length(output) == 0){
    output <- NULL
  } else {
    output <- substr(output, 4, nchar(output))
  }

  return(output)
}

## ---------------------------------------------------------------- ##
#                   uniprot.kegg <- function(id)                     #
## ---------------------------------------------------------------- ##
#' Identifier Mapping From UniProt to KEGG
#' @description Mapping between UniProt and KEGG protein identifiers.
#' @usage uniprot.kegg(id)
#' @param id the identifier to be converted.
#' @return Returns a character string corresponding to the requested identifier.
#' @author Juan Carlos Aledo
#' @examples \dontrun{uniprot.kegg('P01009')}
#' @seealso id.mapping()
#' @export

uniprot.kegg <- function(id){

  if (requireNamespace("KEGGREST", quietly = TRUE)){
    sp <- tryCatch(
      {
        species.mapping(id)
      },
      error = function(cond){
        return(NULL)
      }
    )

    if (is.null(sp)){
      return(NULL)
    } else {
      organisms <- NULL
      tryCatch(
        {
          load(url("https://github.com/jcaledo/kegg_species/blob/master/organisms.Rda?raw=true"))
        },
        error = function(cond){
          return(NULL)
        },
        warning = function(w) conditionMessage(w)
      )
      if (is.null(organisms)){
        return(NULL)
      } else {
        org <- organisms$organism[which(regexpr(sp, organisms$species) != -1)[1]]
      }

      if (is.na(org)){ # if there is not a full match, try a partial match
        sp_ <- strsplit(sp, " ")[[1]]
        sp_ <- paste(sp_[1], sp_[2])
        org <- organisms$organism[which(regexpr(sp_, organisms$species) != -1)[1]]
      }
      output <- tryCatch(
        {
          KEGGREST::keggConv(org, paste('uniprot:', id, sep = ""))
        },
        error = function(cond){
          return(NULL)
        }
      )
    }

  } else {
    message("The package KEGGREST must be installed to get this mapping")
    return(NULL)
  }
  return(output)
}

## ---------------------------------------------------------------- ##
#                   uniprot.pdb <- function(id)                     #
## ---------------------------------------------------------------- ##
#' Identifier Mapping From UniProt to PDB
#' @description Mapping between UniProt and PDB protein identifiers.
#' @usage uniprot.pdb(id)
#' @param id the identifier to be converted.
#' @return Returns a character string corresponding to the requested identifier.
#' @author Juan Carlos Aledo
#' @examples \dontrun{uniprot.pdb('P01009')}
#' @seealso id.mapping()
#' @importFrom httr GET
#' @importFrom httr POST
#' @importFrom httr content
#' @importFrom httr add_headers
#' @importFrom httr accept_json
#' @importFrom utils read.table
#' @export

uniprot.pdb <- function(id){

  # https://www.uniprot.org/help/id_mapping
  isJobReady <- function(jobId) {
    pollingInterval = 5
    nTries = 20
    for (i in 1:nTries) {
      url <- paste("https://rest.uniprot.org/idmapping/status/", jobId, sep = "")
      r <- GET(url = url, accept_json())
      status <- content(r, as = "parsed")
      if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
        return(TRUE)
      }
      if (!is.null(status[["messages"]])) {
        print(status[["messages"]])
        return (FALSE)
      }
      Sys.sleep(pollingInterval)
    }
    return(FALSE)
  }

  getResultsURL <- function(redirectURL) {
    if (grepl("/idmapping/results/", redirectURL, fixed = TRUE)) {
      url <- gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
    } else {
      url <- gsub("/results/", "/results/stream/", redirectURL)
    }
  }


  files = list(
    ids = id,
    from = "UniProtKB_AC-ID",
    to = "PDB"
  )

  my_headers <- httr::add_headers('User-Agent' = paste('R', 'metosite@uma.es'))
  r <- POST(url = "https://rest.uniprot.org/idmapping/run", body = files, my_headers, encode = "multipart", accept_json())
  submission <- content(r, as = "parsed")

  if (isJobReady(submission[["jobId"]])) {
    url <- paste("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]], sep = "")
    r <- GET(url = url, accept_json())
    details <- content(r, as = "parsed")
    url <- getResultsURL(details[["redirectURL"]])
    # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available
    url <- paste(url, "?format=tsv", sep = "")
    r <- GET(url = url, accept_json())
    resultsTable = read.table(text = content(r), sep = "\t", header=TRUE)
    output <- resultsTable$To
  }

  if (!is.null(output)){
    return(output)

  } else {
    message("Uniprot server could not respond")
    output <- NULL
  }

  if (length(output) == 0){
    output <- NULL
  }
  return(output)
}

## ---------------------------------------------------------------- ##
#                   pdb.uniprot <- function(id)                     #
## ---------------------------------------------------------------- ##
#' Identifier Mapping From PDB to UniProt
#' @description Mapping between PDB and UniProt protein identifiers.
#' @usage pdb.uniprot(id)
#' @param id the identifier to be converted.
#' @return Returns a character string corresponding to the requested identifier.
#' @author Juan Carlos Aledo
#' @examples \dontrun{pdb.uniprot('3cwm')}
#' @seealso id.mapping()
#' @importFrom httr GET
#' @importFrom httr POST
#' @importFrom httr content
#' @importFrom httr add_headers
#' @importFrom httr accept_json
#' @importFrom utils read.table
#' @export

pdb.uniprot <- function(id){

  isJobReady <- function(jobId) {
    pollingInterval = 5
    nTries = 20
    for (i in 1:nTries) {
      url <- paste("https://rest.uniprot.org/idmapping/status/", jobId, sep = "")
      r <- GET(url = url, accept_json())
      status <- content(r, as = "parsed")
      if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
        return(TRUE)
      }
      if (!is.null(status[["messages"]])) {
        print(status[["messages"]])
        return (FALSE)
      }
      Sys.sleep(pollingInterval)
    }
    return(FALSE)
  }

  getResultsURL <- function(redirectURL) {
    if (grepl("/idmapping/results/", redirectURL, fixed = TRUE)) {
      url <- gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
    } else {
      url <- gsub("/results/", "/results/stream/", redirectURL)
    }
  }

  files <- list(
    ids = id,
    from = "PDB",
    to = "UniProtKB"
  )

  my_headers <- httr::add_headers('User-Agent' = paste('R', 'metosite@uma.es'))
  r <- POST(url = "https://rest.uniprot.org/idmapping/run", body = files, my_headers, encode = "multipart", accept_json())
  submission <- content(r, as = "parsed")

  if (isJobReady(submission[["jobId"]])) {
    url <- paste("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]], sep = "")
    r <- GET(url = url, accept_json())
    details <- content(r, as = "parsed")
    url <- getResultsURL(details[["redirectURL"]])
    # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available
    url <- paste(url, "?format=tsv", sep = "")
    r <- GET(url = url, accept_json())
    resultsTable = read.table(text = content(r), sep = "\t", header=TRUE)
    output <- resultsTable$Entry
  }

  if (!is.null(output)){
    return(output)

  } else {
    message("Uniprot server could not respond")
    output <- NULL
  }

  if (length(output) == 0){
    output <- NULL
  }
  return(output)
}
