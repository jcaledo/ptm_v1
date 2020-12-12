## ---------------  dssp.R  ------------------ ##
#                                               #
#   parse.dssp                                  #
#   compute.dssp                                #
#   mkdssp                                      #
#                                               #
## ------------------------------------------- ##


## ---------------------------------------------------------------- ##
#      parse.dssp <- function(file, keepfiles = FALSE)               #
## ---------------------------------------------------------------- ##
#' Parse a DSSP File to Return a Dataframe
#' @description Parses a DSSP file to return a dataframe.
#' @usage parse.dssp(file, keepfiles = FALSE)
#' @param file input dssp file.
#' @param keepfiles logical, if TRUE the dataframe will be saved in the working directory and we will keep the dssp file.
#' @details If the argument 'keepfiles' is not set to TRUE, the dssp file used to get the parsed dataframe will be removed.
#' @return Returns a dataframe providing data for 'acc', 'ss', 'phi' and 'psi' for each residues from the structure.
#' @author Juan Carlos Aledo
#' @examples \dontrun{compute.dssp('3cwm'); parse.dssp('3cwm.dssp')}
#' @references Touw et al (2015) Nucl. Ac. Res. 43(Database issue): D364-D368 (PMID: 25352545).
#' @seealso download.dssp(), compute.dssp(), mkdssp() and acc.dssp()
#' @export

parse.dssp <- function(file, keepfiles = FALSE){
  ## --------------- Reading the dssp file ------------------ ##
  con <- file(file, 'r')

  counter <- 0
  resnum <- c()
  respdb <- c()
  chain <- c()
  aa <- c()
  ss <- c()
  sasa <- c()
  phi <- c()
  psi <- c()

  while(TRUE){
    line <- readLines(con, n = 1)
    counter <- counter + 1

    if (counter == 1){
      l <- strsplit(line, split = "")[[1]]
      l <- paste(l, collapse = "")
      if ("have bz2" %in% l){
        first_valid_line <- 29 # dssp file coming from the API
      } else {
        first_valid_line <- 28 # dssp file coming from the sync
      }
    }

    if (counter > first_valid_line & length(line) != 0){
      a <- strsplit(line, split = "")[[1]]
      resnum <- c(resnum, paste(a[1:5], collapse = ""))
      respdb <- c(respdb, paste(a[6:10], collapse = ""))
      chain <- c(chain, paste(a[11:12], collapse = ""))
      aa <- c(aa, paste(a[13:14], collapse = ""))
      ss <- c(ss, paste(a[15:17], collapse = ""))
      sasa <- c(sasa, paste(a[36:38], collapse = ""))
      phi <- c(phi, paste(a[104:109], collapse = ""))
      psi <- c(psi, paste(a[110:115], collapse = ""))
    }
    if (length(line) == 0){
      break
    }
  }
  close(con)

  ## ------ Setting the variable types ------------- ##
  resnum <- as.numeric(resnum)
  respdb <- as.numeric(respdb)
  chain <- gsub(" ", "", chain)
  aa <- gsub(" ", "", aa)
  ss <- gsub("   ", "C", ss)
  ss <- gsub(" ", "", ss)

  ## -------- Building the dataframe ---------------- ##
  df <- as.data.frame(matrix(c(resnum, respdb, chain, aa,
                               ss, sasa, phi, psi), ncol = 8),
                      stringsAsFactors = FALSE)

  colnames(df) <- c('resnum', 'respdb', 'chain', 'aa', 'ss',
                    'sasa', 'phi', 'psi')

  df$resnum <- as.numeric(df$resnum)
  df$respdb <- as.numeric(df$respdb)
  df$sasa <- as.numeric(df$sasa)
  df$phi <- as.numeric(df$phi)
  df$psi <- as.numeric(df$psi)

  if (keepfiles == TRUE){
    save(df, file = paste(file, ".Rda", sep = ""))
  } else {
    file.remove(file)
  }

  ## --------------- Remove empty lines between chains ------------- ##
  badlines <- c()
  for (i in 1:nrow(df)){
    if (df$aa[i] == '!' | df$aa[i] == 'X'){
      badlines <- c(badlines, i)
    }
  }
  if (length(badlines) != 0){
    df <- df[-badlines,]
    df$resnum <- 1:nrow(df)
  }
  return(df)
}

## ---------------------------------------------------------------- ##
#      compute.dssp <- function(pdb, destfile = './')                #
## ---------------------------------------------------------------- ##
#' Compute and Return a DSSP File
#' @description Computes and returns a DSSP file.
#' @usage compute.dssp(pdb, destfile = './')
#' @param pdb is either a PDB id, or the path to a pdb file.
#' @param destfile a character string with the path where the DSSP file is going to be saved.
#' @details A drawback of this function is that it depends on DSSP's server and in occasions it can take a long time to process the request.
#' @return an online computed dssp file that is saved at the indicated location.
#' @author Juan Carlos Aledo
#' @examples \dontrun{compute.dssp(pdb = '3cwm', destfile = './')}
#' @references Touw et al (2015) Nucl. Ac. Res. 43(Database issue): D364-D368 (PMID: 25352545).
#' @seealso download.dssp(), parse.dssp(), mkdssp() and acc.dssp()
#' @importFrom httr content
#' @importFrom httr GET
#' @importFrom httr POST
#' @importFrom httr upload_file
#' @importFrom bio3d get.pdb
#' @export

compute.dssp <- function(pdb, destfile = './'){

  del <- FALSE
  if (nchar(pdb) == 4){ # when input is a PDB ID
    mypdb <- suppressWarnings(bio3d::get.pdb(pdb)) # avoids warning: 'pdb exists. Skipping download'
    file <- paste("./", pdb, ".pdb", sep = "")
    del <- TRUE
  } else {
    file <- pdb
  }

  url_create <- 'https://www3.cmbi.umcn.nl/xssp/api/create/pdb_file/dssp/'
  # url_create <- 'http://www.cmbi.umcn.nl/xssp/api/create/pdb_file/dssp/'
  body <- list(file_ = httr::upload_file(file))
  response_create <- httr::POST(url_create, body = body)

  if (response_create$status_code %in% c(200, 202)){
    job_id <- httr::content(response_create)
  } else {
    stop(print(response_create$status_code))
  }

  url <- 'http://www.cmbi.umcn.nl/xssp/api/status/pdb_file/dssp/'
  url_status <- paste(url, job_id, sep="")
  response_status <- httr::GET(url_status)

  if (response_status$status_code == 200){
    job_status <- httr::content(response_status)
    ready <- FALSE
    attempts <- 0
    while(!ready & attempts < 3){
      print(attempts)
      attempts <- attempts + 1
      if (job_status == 'SUCCESS'){
        ready = TRUE
      } else if (job_status %in% c('FAILURE', 'REVOKED')){
        stop(print(job_status))
      } else {
        Sys.sleep(5)
      }
    }

    if (del){ # if a pdb file was downloaded now is deleted
      file.remove(file)
    }
    if (job_status != 'SUCCESS'){
      stop("After three attempts the server didn't responde")
    }

    t <- strsplit(file, split = "\\/")[[1]] # name and location of the saved file
    # pdb_id <- substring(t[length(t)], 1,4)
    pdb_id <- strsplit(t[length(t)], split = "\\.")[[1]][1]
    destfile = paste(destfile, pdb_id, ".dssp", sep = "")

    url <- 'http://www.cmbi.umcn.nl/xssp/api/result/pdb_file/dssp/'
    url_results <- paste(url, job_id, sep="")
    response_results <- httr::GET(url_results)

    if (response_results$status_code == 200){
      a <- httr::content(response_results)
      cat(a$result, file = destfile)
      return(paste("Work done!. See file at: ", destfile, sep = ""))
    } else {
      stop(print(response_results$status_code))
    }

  } else {
    stop("Response_status fails")
  }
}

## ---------------------------------------------------------------- ##
#         mkdssp <- function(pdb, method, exefile = "dssp")          #
## ---------------------------------------------------------------- ##
#' Compute DSSP File Using an In-House Version of the DSSP Software
#' @description Computes the DSSP file using an in-house version of the DSSP software.
#' @usage mkdssp(pdb, method = 'ptm', exefile = "dssp")
#' @param pdb is either a 4-character identifier of the PDB structure, or the path to a pdb file.
#' @param method a character string specifying the desired method to get the dssp dataframe; it should be one of 'ptm' or 'bio3d'.
#' @param exefile  file path to the DSSP executable on your system (i.e. how is DSSP invoked).
#' @details The structure of the output data depends on the method chosen, but it will always contain the DSSP-related data.
#' @return Returns either a dataframe containing the information extracted from the dssp file (method ptm), or a list with that information (method bio3d).
#' @author Juan Carlos Aledo
#' @examples \dontrun{mkdssp('3cwm', method = 'ptm')}
#' @references Touw et al (2015) Nucl. Ac. Res. 43(Database issue): D364-D368 (PMID: 25352545).
#' @seealso download.dssp(), parse.dssp(), compute.dssp() and acc.dssp()
#' @importFrom bio3d dssp
#' @importFrom bio3d read.pdb
#' @importFrom bio3d get.pdb
#' @export

mkdssp <- function(pdb, method = 'ptm', exefile = "dssp"){

  ## --- pdb id or path to its file
  del <- FALSE
  if (nchar(pdb) == 4){ # when input is a PDB ID
    suppressWarnings(bio3d::get.pdb(pdb)) # avoids warning: 'pdb exists. Skipping download'
    file <- paste("./", pdb, ".pdb", sep = "")
    del <- TRUE
  } else {
    file <- pdb
  }

  ## ---- Path to DSSP executable
  exefile <- .get.exepath(exefile)
  success <- .test.exefile(exefile)

  if (!success){
    stop(paste("Launching external program 'dssp' (or 'mkdssp') failed\n",
               "  make sure '", exefile, "' is in your search path", sep = ""))
  }

  ## ------ Method to compute DSSP
  if(method == 'ptm'){

    mydssp <- system(paste(exefile, " -i ", file, " -o ./temp.dssp", sep = ""))
    mydssp <- parse.dssp('./temp.dssp')

  } else if (method == 'bio3d'){
    mydssp <- suppressWarnings(bio3d::dssp(read.pdb(file) , exefile = exefile)) # avoids warning: Non-protein residues detected in input PDB.
  }

  if (del){ # if a pdb file was downloaded now is deleted
    file.remove(file)
  }
  return(mydssp)
}
