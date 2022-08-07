## ---------------------------------------------------------------- ##
#      compute.dssp <- function(pdb, destfile = './')                #
## ---------------------------------------------------------------- ##
#' Compute and Return a DSSP File
#' @description Computes and returns a DSSP file.
#' @usage compute.dssp(pdb, destfile = './')
#' @param pdb is either a PDB id, or the path to a pdb file.
#' @param destfile a character string with the path where the DSSP file is going to be saved.
#' @details A drawback of this function is that it depends on DSSP's server and in occasions it can take a long time to process the request.
#' @return An online computed dssp file that is saved at the indicated location.
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
  .Deprecated(new = "mkdssp", old = "compute.dssp")
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
  body <- tryCatch(
    {
      list(file_ = httr::upload_file(file))
    },
    error = function(cond){
      return(NULL)
    }
  )
  if (is.null(body)){
    message("Sorry, no file was found to be uploaded")
    return(NULL)
  }

  response_create <- tryCatch(
    {
      httr::POST(url_create, body = body)
    },
    error = function(cond){
      message(cond)
      return(NULL)
    }
  )
  if (is.null(response_create)){
    message("Sorry, pdb file could be posted")
    return(NULL)
  }

  if (response_create$status_code %in% c(200, 202)){
    job_id <- httr::content(response_create)
  }

  url <- 'http://www.cmbi.umcn.nl/xssp/api/status/pdb_file/dssp/'
  url_status <- paste(url, job_id, sep="")

  response_status <- tryCatch(
    {
      httr::GET(url_status)
    },
    error = function(cond){
      message(cond)
      return(NULL)
    }
  )

  if (is.null(response_status)){
    message("Sorry, no response from the DSSP server")
    return(NULL)
  }

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
      message("After three attempts the server didn't responde")
      return(NULL)
    }

    t <- strsplit(file, split = "\\/")[[1]] # name and location of the saved file
    # pdb_id <- substring(t[length(t)], 1,4)
    pdb_id <- strsplit(t[length(t)], split = "\\.")[[1]][1]
    destfile = paste(destfile, pdb_id, ".dssp", sep = "")

    url <- 'http://www.cmbi.umcn.nl/xssp/api/result/pdb_file/dssp/'
    url_results <- paste(url, job_id, sep="")

    response_results <- tryCatch(
      {
        httr::GET(url_results)
      },
      error = function(cond){
        message(cond)
        return(NULL)
      }
    )
    if (is.null(response_results)){
      return(NULL)
    }

    if (response_results$status_code == 200){
      a <- httr::content(response_results)
      cat(a$result, file = destfile)
      return(paste("Work done!. See file at: ", destfile, sep = ""))
    } else {
      message(paste("Sorry, ", response_results$status_code))
      return(NULL)
    }

  } else {
    message("Response_status fails")
    return(NULL)
  }
}
