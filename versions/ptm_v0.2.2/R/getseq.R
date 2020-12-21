## ---------- getseq.R ------------ ##
#                                    #
#     get.seq                        #
#     prot2codon                     #
#     id.mapping                     #
#     id.features                    #
#     species.mapping                #
#     species.kegg                   #
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
#' @examples \dontrun{get.seq("hsa:5265", db = "kegg-aa")}
#' @examples \dontrun{get.seq("1u8f:P", db = "pdb")}
#' @importFrom httr GET
#' @importFrom httr content
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
    df <- pdb.seq(id[1])
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

  } else if (db == 'kegg-aa' | db == 'kegg-nt'){
    t <- paste(strsplit(db, split = '-')[[1]][2], 'seq', sep = '')
    if (requireNamespace('KEGGREST', quietly = TRUE)){
      text <- tryCatch(
        {
          as.character(KEGGREST::keggGet(id, t))
        },
        error = function(cond){
          # message(cond)
          return(paste("Sorry, the entry ", id, " couldn't be retrieved!"))
        }
      )
    } else {
      stop("Please, install the KEGGREST pakage to get this functionality")
    }
    call <- NULL

  } else{
    stop('You should indicate a proper DB')
  }

  ## -------- Client <-> Server Communication -------- ##
  if (!is.null(call)){
    resp <- try(httr::GET(call), silent = TRUE)
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


  ## -------------- Parsing the response --------------- ##
  if (text != "LOST CONNECTION"){
    if (db == 'uniprot'){
      seq <- strsplit(text, split = "\\n")[[1]][-1]
      seq <- paste(seq, collapse = "")
    } else if (db == 'metosite'){
      if (requireNamespace("jsonlite", quietly = TRUE)){
        data <- jsonlite::fromJSON(text, flatten = TRUE)
        seq <- data$prot_seq
      } else {
        warning("If you install the package 'jsonlite' this ugly output would be nicely parsed")
        return(text)
      }
    } else if (db == 'kegg-aa' | db == 'kegg-nt'){
      seq <- text
    }

    if (is.null(seq) || is.na(seq)){
      output <- paste("The entry", id, "is not found in the", toupper(db), "database")
    } else {
      output <- seq
      if (!as.string){
        output <- strsplit(seq, split ="")
      }
    }
  } else {
    output <- text
  }

  attr(output, "ID") <- id
  attr(output, "DB") <- db
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
#' @examples prot2codon('P01009')
#' @examples \dontrun{prot2codon(prot = '1ATU', chain = 'A')}
#' @importFrom bio3d read.pdb
#' @importFrom bio3d aa321
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
    seq <- ptm::get.seq(id, db = 'pdb', as.string = FALSE)[[1]]

    t <- prot
    source <- 'pdb'

  } else { # input should be a uniprot ID
    seq <- ptm::get.seq(prot, as.string = FALSE)[[1]]
    if (length(seq) == 0) {return("Sorry, no Uniprot seq could be found!")}
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
  organism <- species.mapping(t, db = source)

  if (source == 'pdb'){
    t <- pdb2uniprot(t, chain = chain)
  }

  kegg_id <- id.mapping(t, from = 'uniprot', to = 'kegg')[1]
  if (kegg_id[1] == "Sorry, no KEGG ID could be found!"){
    output <- "Sorry, we couldn't get a DNA seq from KEGG"
    return(output)
  }
  dna <- ptm::get.seq(kegg_id, 'kegg-nt')
  codon <- gsub("(.{3})", "\\1 ", dna)
  codon <- strsplit(codon, split = " ")[[1]]
  if (requireNamespace('seqinr', quietly = TRUE)){
    translated <- seqinr::translate(seqinr::s2c(dna))
  } else {
    stop("The package seqinr must be installed to use this function")
  }

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

  ## ------------ Cheking the translation -------------------- ##
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
  if (!laxity){
    if (sum(output$check) != nrow(output)) { stop("Translation problem")}
  }
  return(output)
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
#' @examples id.mapping('P01009', from = 'uniprot', to = 'pdb')
#' @seealso pdb2uniprot(), uniprot2pdb()
#' @importFrom httr GET
#' @importFrom httr content
#' @importFrom httr add_headers
#' @export

id.mapping <- function(id, from, to){

  ## ----- From KEGG to UniProt ------------- ##
  kegg.uniprot <- function(id){
    if (requireNamespace("KEGGREST", quietly = TRUE)){
      output <- KEGGREST::keggConv('uniprot', id)
    } else {
      stop("The package KEGGREST must be installed to get this mapping")
    }
    output <- substr(output, 4, nchar(output))
    return(output)
  }
  ## ----- From UniProt to KEGG ------------- ##
  uniprot.kegg <- function(id){
    if (requireNamespace("KEGGREST", quietly = TRUE)){
      sp <- species.mapping(id)
      organisms <- NULL
      load(url("https://github.com/jcaledo/kegg_species/blob/master/organisms.Rda?raw=true"))
      org <- organisms$organism[which(regexpr(sp, organisms$species) != -1)[1]]
      if (is.na(org)){ # if there is not a full match, try a partial match
        sp_ <- strsplit(sp, " ")[[1]]
        sp_ <- paste(sp_[1], sp_[2])
        org <- organisms$organism[which(regexpr(sp_, organisms$species) != -1)[1]]
      }
      output <- KEGGREST::keggConv(org, paste('uniprot:', id, sep = ""))
    } else {
      stop("The package KEGGREST must be installed to get this mapping")
    }
    return(output)
  }

  uniprot_url <- "http://www.uniprot.org/uploadlists/"
  my_headers <- httr::add_headers('User-Agent' = paste('R', 'metosite@uma.es'))

  ## ----- From UniProt to PDB ------------- ##
  uniprot.pdb <- function(id){
    params <- list(from = 'ACC',
                   to = 'PDB_ID',
                   format = 'tab',
                   query = id)
    request <- try(httr::GET(uniprot_url, query = params, my_headers),
                   silent = TRUE)
    if (inherits(request, "try-error")) {
      ans <- "LOST CONNECTION"
    } else {
      ans <- httr::content(request, 'text', encoding = "ISO-8859-1")
    }
    retry = 0
    while ((grepl("^LOST CONNECTION", ans) || httr::http_error(request) ||
            grepl("^ERROR", ans) || grepl("Nothing has been found",
                                          ans)) && retry < 3) {
      retry = retry + 1
      Sys.sleep(5)
      request <- try(httr::GET(uniprot_url, query = params, my_headers),
                     silent = TRUE)
      if (inherits(request, "try-error")) {
        ans <- "LOST CONNECTION"
      }
      else {
        ans <- httr::content(request, 'text', encoding = "ISO-8859-1")
      }
    }

    pos <- gregexpr(id, ans)[[1]]
    output <- character(length(pos))
    for (i in seq_len(length(pos))){
      n_id <- (substr(ans, pos[i]+nchar(id)+1, pos[i]+nchar(id)+4))
      output[i] <- n_id
    }
    return(output)
  }

  ## ----- From PDB to UniProt ------------- ##
  pdb.uniprot <- function(id){
    params <- list(from = 'PDB_ID',
                   to = 'ACC',
                   format = 'tab',
                   query = id)
    request <- try(httr::GET(uniprot_url, query = params, my_headers),
                   silent = TRUE)
    if (inherits(request, "try-error")) {
      ans <- "LOST CONNECTION"
    } else {
      ans <- httr::content(request, 'text', encoding = "ISO-8859-1")
    }
    retry = 0
    while ((grepl("^LOST CONNECTION", ans) || httr::http_error(request) ||
            grepl("^ERROR", ans) || grepl("Nothing has been found",
                                          ans)) && retry < 3) {
      retry = retry + 1
      Sys.sleep(5)
      request <- try(httr::GET(uniprot_url, query = params, my_headers),
                     silent = TRUE)
      if (inherits(request, "try-error")) {
        ans <- "LOST CONNECTION"
      }
      else {
        ans <- httr::content(request, 'text', encoding = "ISO-8859-1")
      }
    }


    temp <- strsplit(ans, split = "\n")[[1]][-1]
    output <- character(length(temp))
    for (i in seq_len(length(temp))){
      t <- strsplit(temp[i], split = "\t")[[1]][2]
      output[i] <- t

    }
    return(output)
  }

  ## -------- Mapping using the above functions --------- ##
  if (from == 'uniprot' & to == 'pdb'){
    output <- uniprot.pdb(id)
    if (output[1] == "" | output[1] == "To\n"){
      return("No PDB found")
    } else {
      return(output)
    }
  } else if (from == 'uniprot' & to == 'kegg'){
    output <- uniprot.kegg(id)
    if (identical(output, character(0))){
      output <- 'Sorry, no KEGG ID could be found!'
    }
    return(output)
  } else if (from == 'kegg' & to == 'uniprot' ){
    output <- kegg.uniprot(id)
    return(output)
  } else if (from == 'kegg' & to == 'pdb'){
    output <- kegg.uniprot(id)
    output <- uniprot.pdb(output[1])
    return(output)
  } else if (from == 'pdb' & to == 'uniprot'){
    output <- pdb.uniprot(id)
    return(output)
  } else if (from == 'pdb' & to == 'kegg'){
    output <- pdb.uniprot(id)
    output <- lapply(output, uniprot.kegg)
    if (identical(output, character(0))){
      output <- 'Sorry, no KEGG ID could be found!'
    }
    attr(output, 'from') <- from
    attr(output, 'to') <- to
    return(output)
  }
}


## ---------------------------------------------------------------- ##
#           id.features <- function(id, features = "")               #
## ---------------------------------------------------------------- ##
#' Features Related to the Protein Entry
#' @description Obtains features related to the provided id.
#' @usage id.features(id, features = "")
#' @param id the UniProt identifier of the protein of interest.
#' @param features a string identifying the features (comma separated) to be recovered.
#' @details By default the the function provides info regarding the following features: id, reviewed, entry name and organism. If wished, this list of features can be expanded using the argument 'features'. There is a larga list of features that can be retrieved. You can look up your relevant feature's name in the full list of UniProtKB found at https://www.uniprot.org/help/uniprotkb_column_names.
#' @return Returns a named list with the requested features.
#' @author Juan Carlos Aledo
#' @examples \dontrun{id.features('P04406', features = 'ec,keywords,database(PDB)')}
#' @seealso species.mapping()
#' @importFrom httr GET
#' @importFrom httr content
#' @importFrom httr add_headers
#' @export

id.features <- function(id, features = "" ){

  ## -------- Building the query ---------- ##
  col <- 'id,reviewed,entry name,organism'
  if (features != ""){
    features <- gsub(', ', ',', features)
    col <- paste(col, features, sep = ',')
  }
  base_url <- "http://www.uniprot.org/uploadlists/"
  my_headers <- httr::add_headers('User-Agent' = paste('R', 'metosite@uma.es'))
  params <- list(from = 'ACC+ID',
                 to = 'ACC',
                 format = 'tab',
                 columns = col,
                 query = id)

  ## --------- Sending the query ----------- ##
  request <- try(httr::GET(base_url, query = params, my_headers),
                 silent = TRUE)
  if (inherits(request, "try-error")) {
    text <- "LOST CONNECTION"
  } else {
    cont <- httr::content(request, "text", encoding = "ISO-8859-1")
  }
  retry = 0
  while ((grepl("^LOST CONNECTION", cont) || httr::http_error(request) ||
          grepl("^ERROR", cont) || grepl("Nothing has been found",
                                         cont)) && retry < 3) {
    retry = retry + 1
    Sys.sleep(5)
    request <- try(httr::GET(base_url), silent = TRUE)
    if (inherits(request, "try-error")) {
      cont <- "LOST CONNECTION"
    }
    else {
      cont <- httr::content(request, "text", encoding = "ISO-8859-1")
    }
  }

  ## ------- Parsing the response ---------- ##
  l <- strsplit(cont, split = "\n")[[1]]
  headers <- strsplit(l, split = "\t")[[1]]
  headers <- gsub(" ", "_", headers)
  output <- as.list(strsplit(l[2], split = "\t")[[1]])
  names(output) <- gsub(" ", "_", headers)
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
#' @examples species.mapping('P01009')
#' @seealso id.features()
#' @importFrom httr GET
#' @importFrom httr content
#' @export

species.mapping <- function(id, db = 'uniprot'){

  db <- tolower(db)
  if (!(db %in% c("uniprot", "pdb"))){
    stop("Option db should be one of uniprot or pdb")
  }
  if (db == 'pdb'){
    id <- id.mapping(id, from = 'pdb', to = 'uniprot')[1]
  }

  url <- paste("https://www.uniprot.org/uniprot/", id, ".fasta", sep = "")
  resp <- try(httr::GET(url), silent = TRUE)
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

  start <- gregexpr("OS=", text)[[1]]
  stop <- gregexpr("OX=", text)[[1]]
  species <- substr(text, start+3, stop-2)

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
  load(url("https://github.com/jcaledo/kegg_species/blob/master/organisms.Rda?raw=true"))
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
      output <- paste('Sorry, could find a code for the organism:', organism)
  }

  return(output)

}
