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
  con <- tryCatch(
    {
      file(file, 'r')
    },
    error = function(cond){
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(con) | grepl("no fue posible", con) | grepl("No such file", con)){
    message("Sorry, no connection could be established")
    return(NULL)
  }

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

