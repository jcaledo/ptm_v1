## ---------------------------------------------------------------- ##
#      get.seq <- function(id, db = 'uniprot', as.string = TRUE)     #
## ---------------------------------------------------------------- ##
#' Import a Protein Sequence from a Database
#' @description Imports a protein sequence from a selected database.
#' @usage get.seq(id, db = 'uniprot', as.string = TRUE)
#' @param id the identifier of the protein of interest.
#' @param db a character string specifying the desired database; it must be one of 'uniprot', 'metosite', 'ncbi','pdb', 'kegg-aa', 'kegg-nt'.
#' @param as.string logical, if TRUE the imported sequence will be returned as a character string.
#' @details MetOSite and NCBI use the same type of protein ID than UniProt. However, if the chosen database is PDB, the identifier should be the 4-character unique identifier characteristic of PDB, followed by colon and the chain of interest. For instance, '2OCC:B' means we are interested in the sequence of chain B from the structure 2OCC. KEGG used its own IDs (see examples).
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

  } else if (db == 'ncbi'){
    baseUrl <- "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=protein"
    call <- paste(baseUrl, "&val=", id, "&report=fasta", sep = "")

  } else if (db == 'pdb'){
    id <- strsplit(id, split=':')[[1]]
    baseUrl <- "https://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList="
    call <- paste(baseUrl, toupper(id[1]), "&compressionType=uncompressed", sep = "")

  } else if (db == 'kegg-aa' | db == 'kegg-nt'){
    t <- paste(strsplit(db, split = '-')[[1]][2], 'seq', sep = '')
    if (requireNamespace('KEGGREST', quietly = TRUE)){
      text <- as.character(KEGGREST::keggGet(id, t))
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
    } else if (db == 'pdb'){
      ID <- toupper(id)
      ID <- paste(ID, collapse = ":")
      if (nchar(ID) == 4){ # No chain selector
        seq <- strsplit(text, split = "\n")[[1]]
        seq <- paste(seq[!grepl("^>", seq)], collapse = "")
      } else if (nchar(ID) == 6){ # chain needs to be selected
        seq <- strsplit(text, split = ">")[[1]][-1]
        seq <- seq[grepl(ID, seq)]
        seq <- strsplit(seq, "\n")[[1]][-1]
        seq <- paste(seq, collapse = "")
      } else {
        stop("The PDB ID must be formated properly")
      }
    } else if (db == 'ncbi'){
      seq <- strsplit(text, "\n")[[1]]
      seq <- paste(seq[-1], collapse = "")
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
