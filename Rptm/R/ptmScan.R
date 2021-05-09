## ------- ptmScan.R -------------- ##
#                                    #
#   p.scan                           #
#   ac.scan                          #
#   me.scan                          #
#   ub.scan                          #
#   su.scan                          #
#   gl.scan                          #
#   sni.scan                         #
#   ni.scan                          #
#   ptm.scan                         #
#   reg.scan                         #
#   dis.scan                         #
#                                    #
## -------------------------------- ##

## ---------------------------------------------------------------- ##
#           p.scan <- function(up_id, db = 'all')                    #
## ---------------------------------------------------------------- ##
#' Scan a Protein in Search of  Phosphosites
#' @description Scans the indicated protein in search of phosphosites.
#' @usage p.scan(up_id, db = 'all')
#' @param up_id a character string corresponding to the UniProt ID.
#' @param db the database where to search. It should be one among 'PSP', 'dbPTM','dbPAF', 'PhosPhAt', 'Phospho.ELM', 'all'.
#' @details If db = 'all' has been selected, it may happen that the same residue appears in several rows if it is present in different databases.
#' @return Returns a dataframe where each row corresponds to a phosphorylatable residue.
#' @author Juan Carlos Aledo
#' @examples \dontrun{p.scan('P01009', db = 'PSP')}
#' @references Hornbeck et al. Nucleic Acids Res. 2019 47:D433-D441, (PMID: 30445427).
#' @references Huang et al. Nucleic Acids Res. 2019 47:D298-D308, (PMID: 30418626).
#' @references Ullah et al. Sci. Rep. 2016 6:23534, (PMID: 27010073).
#' @references Durek et al. Nucleic Acids Res.2010 38:D828-D834, (PMID: 19880383).
#' @references Dinkel et al. Nucleic Acids Res. 2011 39:D261-D567 (PMID: 21062810).
#' @seealso meto.scan(), ac.scan(), me.scan(), ub.scan(), su.scan(), gl.scan(), sni.scan(), ni.scan(), ptm.scan(), reg.scan(), dis.scan()
#' @export

p.scan <- function(up_id, db = 'all'){

  if (! db %in% c('PSP', 'dbPTM','dbPAF', 'PhosPhAt', 'Phospho.ELM', 'all')){
    stop("Please, select an appropiated database!")
  }

  p_db <- NULL
  baseUrl <- "https://github.com/jcaledo/p_db/blob/master/"
  call <- paste(baseUrl, "p_db_", up_id, ".Rda?raw=true", sep = "")

  tryCatch(
    {
      load(url(call))
    },
    error = function(cond){
      message("Sorry, no result could be retrieved")
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(p_db)){
    message("Sorry, no result could be retrieved")
    return(NULL)
  }

  if (db != 'all'){
    p_db <- p_db[which(p_db$database == db),]
  }

  return(p_db)
}


## ---------------------------------------------------------------- ##
#               ac.scan <- function(up_id, db = 'all')               #
## ---------------------------------------------------------------- ##
#' Scan a Protein in Search of Acetylation Sites
#' @description Scans the indicated protein in search of acetylation sites.
#' @usage ac.scan(up_id, db = 'all')
#' @param up_id a character string corresponding to the UniProt ID.
#' @param db the database where to search. It should be one among 'PSP', 'dbPTM', 'all'.
#' @details If db = 'all' has been selected, it may happen that the same residue appears in several rows if it is present in different databases.
#' @return Returns a dataframe where each row corresponds to an acetylable residue.
#' @author Juan Carlos Aledo
#' @examples \dontrun{ac.scan('P01009', db = 'PSP')}
#' @references Hornbeck et al. Nucleic Acids Res. 2019 47:D433-D441, (PMID: 30445427).
#' @references Huang et al. Nucleic Acids Res. 2019 47:D298-D308, (PMID: 30418626).
#' @seealso meto.scan(), p.scan(), me.scan(), ub.scan(), su.scan(), gl.scan(), sni.scan(), ni.scan(), ptm.scan(), reg.scan(), dis.scan()
#' @export

ac.scan <- function(up_id, db = 'all'){

  if (! db %in% c('PSP', 'dbPTM','dbPAF', 'PhosPhAt', 'Phospho.ELM', 'all')){
    stop("Please, select an appropiated database!")
  }

  ac_db <- NULL
  baseUrl <- "https://github.com/jcaledo/ac_db/blob/master/"
  call <- paste(baseUrl, "ac_db_", up_id, ".Rda?raw=true", sep = "")

  tryCatch(
    {
      load(url(call))
    },
    error = function(cond){
      message("Sorry, no result could be retrieved")
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(ac_db)){
    message("Sorry, no result could be retrieved")
    return(NULL)
  }

  if (db != 'all'){
    ac_db <- ac_db[which(ac_db$database == db),]
  }

  return(ac_db)
}

## ---------------------------------------------------------------- ##
#               me.scan <- function(up_id, db = 'all')               #
## ---------------------------------------------------------------- ##
#' Scan a Protein in Search of Methylation Sites
#' @description Scans the indicated protein in search of methylation sites.
#' @usage me.scan(up_id, db = 'all')
#' @param up_id a character string corresponding to the UniProt ID.
#' @param db the database where to search. It should be one among 'PSP', 'dbPTM', 'all'.
#' @details If db = 'all' has been selected, it may happen that the same residue appears in several rows if it is present in different databases.
#' @return Returns a dataframe where each row corresponds to a modifiable residue.
#' @author Juan Carlos Aledo
#' @references Hornbeck et al. Nucleic Acids Res. 2019 47:D433-D441, (PMID: 30445427).
#' @references Huang et al. Nucleic Acids Res. 2019 47:D298-D308, (PMID: 30418626).
#' @examples \dontrun{me.scan('Q16695', db = 'PSP')}
#' @seealso meto.scan(), ac.scan(), p.scan(), ub.scan(), su.scan(), gl.scan(), sni.scan(), ni.scan(), ptm.scan(), reg.scan(), dis.scan()
#' @export

me.scan <- function(up_id, db = 'all'){

  if (! db %in% c('PSP', 'dbPTM','dbPAF', 'PhosPhAt', 'Phospho.ELM', 'all')){
    stop("Please, select an appropiated database!")
  }

  me_db <- NULL
  baseUrl <- "https://github.com/jcaledo/me_db/blob/master/"
  call <- paste(baseUrl, "me_db_", up_id, ".Rda?raw=true", sep = "")

  tryCatch(
    {
      load(url(call))
    },
    error = function(cond){
      message("Sorry, no result could be retrieved")
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(me_db)){
    message("Sorry, no result could be retrieved")
    return(NULL)
  }

  if (db != 'all'){
    me_db <- me_db[which(me_db$database == db),]
  }

  return(me_db)
}

## ---------------------------------------------------------------- ##
#               ub.scan <- function(up_id, db = 'all')               #
## ---------------------------------------------------------------- ##
#' Scan a Protein in Search of Ubiquitination Sites
#' @description Scans the indicated protein in search of ubiquitination sites.
#' @usage ub.scan(up_id, db = 'all')
#' @param up_id a character string corresponding to the UniProt ID.
#' @param db the database where to search. It should be one among 'PSP', 'dbPTM', 'all'.
#' @details If db = 'all' has been selected, it may happen that the same residue appears in several rows if it is present in different databases.
#' @return Returns a dataframe where each row corresponds to a modifiable residue.
#' @author Juan Carlos Aledo
#' @examples \dontrun{ub.scan('Q16695', db = 'PSP')}
#' @references Hornbeck et al. Nucleic Acids Res. 2019 47:D433-D441, (PMID: 30445427).
#' @references Huang et al. Nucleic Acids Res. 2019 47:D298-D308, (PMID: 30418626).
#' @seealso meto.scan(), ac.scan(), me.scan(), p.scan(), su.scan(), gl.scan(), sni.scan(), ni.scan(), ptm.scan(), reg.scan(), dis.scan()
#' @export

ub.scan <- function(up_id, db = 'all'){

  if (! db %in% c('PSP', 'dbPTM','dbPAF', 'PhosPhAt', 'Phospho.ELM', 'all')){
    stop("Please, select an appropiated database!")
  }

  ub_db <- NULL
  baseUrl <- "https://github.com/jcaledo/ub_db/blob/master/"
  call <- paste(baseUrl, "ub_db_", up_id, ".Rda?raw=true", sep = "")

  tryCatch(
    {
      load(url(call))
    },
    error = function(cond){
      message("Sorry, no result could be retrieved")
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(ub_db)){
    message("Sorry, no result could be retrieved")
    return(NULL)
  }

  if (db != 'all'){
    ub_db <- ub_db[which(ub_db$database == db),]
  }

  return(ub_db)
}

## ---------------------------------------------------------------- ##
#               su.scan <- function(up_id, db = 'all')               #
## ---------------------------------------------------------------- ##
#' Scan a Protein in Search of Sumoylation Sites
#' @description Scans the indicated protein in search of sumoylation sites.
#' @usage su.scan(up_id, db = 'all')
#' @param up_id a character string corresponding to the UniProt ID.
#' @param db the database where to search. It should be one among 'PSP', 'dbPTM', 'all'.
#' @details If db = 'all' has been selected, it may happen that the same residue appears in several rows if it is present in different databases.
#' @return Returns a dataframe where each row corresponds to a modifiable residue.
#' @author Juan Carlos Aledo
#' @examples \dontrun{su.scan('Q16695', db = 'PSP')}
#' @references Hornbeck et al. Nucleic Acids Res. 2019 47:D433-D441, (PMID: 30445427).
#' @references Huang et al. Nucleic Acids Res. 2019 47:D298-D308, (PMID: 30418626).
#' @seealso meto.scan(), ac.scan(), me.scan(), ub.scan(), p.scan(), gl.scan(), sni.scan(), ni.scan(), ptm.scan(), reg.scan(), dis.scan()
#' @export

su.scan <- function(up_id, db = 'all'){

  if (! db %in% c('PSP', 'dbPTM','dbPAF', 'PhosPhAt', 'Phospho.ELM', 'all')){
    stop("Please, select an appropiated database!")
  }

  su_db <- NULL
  baseUrl <- "https://github.com/jcaledo/su_db/blob/master/"
  call <- paste(baseUrl, "su_db_", up_id, ".Rda?raw=true", sep = "")

  tryCatch(
    {
      load(url(call))
    },
    error = function(cond){
      message("Sorry, no result could be retrieved")
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(su_db)){
    message("Sorry, no result could be retrieved")
    return(NULL)
  }

  if (db != 'all'){
    su_db <- su_db[which(su_db$database == db),]
  }

  return(su_db)
}

## ---------------------------------------------------------------- ##
#               gl.scan <- function(up_id, db = 'all')               #
## ---------------------------------------------------------------- ##
#' Scan a Protein in Search of OGlcNAc Sites
#' @description Scans the indicated protein in search of glycosylation sites.
#' @usage gl.scan(up_id, db = 'all')
#' @param up_id a character string corresponding to the UniProt ID.
#' @param db the database where to search. It should be one among 'PSP', 'dbPTM', 'all'.
#' @details If db = 'all' has been selected, it may happen that the same residue appears in several rows if it is present in different databases.
#' @return Returns a dataframe where each row corresponds to a modifiable residue.
#' @author Juan Carlos Aledo
#' @examples \dontrun{gl.scan('P08670', db = 'PSP')}
#' @references Hornbeck et al. Nucleic Acids Res. 2019 47:D433-D441, (PMID: 30445427).
#' @references Huang et al. Nucleic Acids Res. 2019 47:D298-D308, (PMID: 30418626).
#' @seealso meto.scan(), ac.scan(), me.scan(), ub.scan(), su.scan(), p.scan(), sni.scan(), ni.scan(), ptm.scan(), reg.scan(), dis.scan()
#' @export

gl.scan <- function(up_id, db = 'all'){

  if (! db %in% c('PSP', 'dbPTM','dbPAF', 'PhosPhAt', 'Phospho.ELM', 'all')){
    stop("Please, select an appropiated database!")
  }

  gl_db <- NULL
  baseUrl <- "https://github.com/jcaledo/gl_db/blob/master/"
  call <- paste(baseUrl, "gl_db_", up_id, ".Rda?raw=true", sep = "")

  tryCatch(
    {
      load(url(call))
    },
    error = function(cond){
      message("Sorry, no result could be retrieved")
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(gl_db)){
    message("Sorry, no result could be retrieved")
    return(NULL)
  }

  if (db != 'all'){
    gl_db <- gl_db[which(gl_db$database == db),]
  }

  return(gl_db)
}

## ---------------------------------------------------------------- ##
#            sni.scan  <- function(up_id", db = 'all')               #
## ---------------------------------------------------------------- ##
#' Scan a Protein in Search of S-nitrosylation Sites
#' @description Scans the indicated protein in search of S-nitrosylation sites.
#' @usage sni.scan(up_id, db = 'all')
#' @param up_id a character string corresponding to the UniProt ID.
#' @param db the database where to search. It should be one among 'PSP', 'dbPTM', 'all'.
#' @return Returns a dataframe where each row corresponds to a modifiable residue.
#' @author Juan Carlos Aledo
#' @examples \dontrun{sni.scan('P01009')}
#' @references Huang et al. Nucleic Acids Res. 2019 47:D298-D308, (PMID: 30418626).
#' @seealso meto.scan(), ac.scan(), me.scan(), ub.scan(), su.scan(), gl.scan(), p.scan(), ni.scan(), ptm.scan(), reg.scan(), dis.scan()
#' @export

sni.scan <- function(up_id, db = 'all'){

  if (! db %in% c('PSP', 'dbPTM','dbPAF', 'PhosPhAt', 'Phospho.ELM', 'all')){
    stop("Please, select an appropiated database!")
  }

  sni_db <- NULL
  baseUrl <- "https://github.com/jcaledo/sni_db/blob/master/"
  call <- paste(baseUrl, "sni_db_", up_id, ".Rda?raw=true", sep = "")

  tryCatch(
    {
      load(url(call))
    },
    error = function(cond){
      message("Sorry, no result could be retrieved")
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(sni_db)){
    message("Sorry, no result could be retrieved")
    return(NULL)
  }

  if (db != 'all'){
    sni_db <- sni_db[which(sni_db$database == db),]
  }

  return(sni_db)
}

## ---------------------------------------------------------------- ##
#            ni.scan  <- function(up_id, db = 'all')                 #
## ---------------------------------------------------------------- ##
#' Scan a Protein in Search of Nitration Sites
#' @description Scans the indicated protein in search of nitration sites.
#' @usage ni.scan(up_id, db = 'all')
#' @param up_id a character string corresponding to the UniProt ID.
#' @param db the database where to search. It should be one among 'PSP', 'dbPTM', 'all'.
#' @return Returns a dataframe where each row corresponds to a modified residue.
#' @author Juan Carlos Aledo
#' @examples \dontrun{ni.scan('P05202')}
#' @references Huang et al. Nucleic Acids Res. 2019 47:D298-D308, (PMID: 30418626).
#' @seealso meto.scan(), ac.scan(), me.scan(), ub.scan(), su.scan(), gl.scan(), sni.scan(), p.scan(), ptm.scan(), reg.scan(), dis.scan()
#' @export

ni.scan <- function(up_id, db = 'all'){

  if (! db %in% c('PSP', 'dbPTM','dbPAF', 'PhosPhAt', 'Phospho.ELM', 'all')){
    stop("Please, select an appropiated database!")
  }

  ni_db <- NULL
  baseUrl <- "https://github.com/jcaledo/ni_db/blob/master/"
  call <- paste(baseUrl, "ni_db_", up_id, ".Rda?raw=true", sep = "")

  tryCatch(
    {
      load(url(call))
    },
    error = function(cond){
      message("Sorry, no result could be retrieved")
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(ni_db)){
    message("Sorry, no result could be retrieved")
    return(NULL)
  }

  if (db != 'all'){
    ni_db <- ni_db[which(ni_db$database == db),]
  }

  return(ni_db)
}


## ---------------------------------------------------------------- ##
#           ptm.scan <- function(up_id, renumerate = TRUE)           #
## ---------------------------------------------------------------- ##
#' Scan a Protein in Search of PTM Sites
#' @description Scans the indicated protein in search of PTM sites.
#' @usage ptm.scan(up_id, renumerate = TRUE)
#' @param up_id a character string corresponding to the UniProt ID.
#' @param renumerate logical, when TRUE the  sequence numeration of MetO sites is that given by Uniprot, which may not coincide with that from MetOSite.
#' @details The numerations of the sequences given by UniProt and MetOSite may or may not match. Sometimes one of the sequences corresponds to the precursor protein and the other to the processed mature protein.
#' @return Returns a dataframe where each row corresponds to a residue, and the columns inform about the modifications.
#' @author Juan Carlos Aledo
#' @examples \dontrun{ptm.scan('P01009', renumerate = TRUE)}
#' @references Hornbeck et al. Nucleic Acids Res. 2019 47:D433-D441, (PMID: 30445427).
#' @references Huang et al. Nucleic Acids Res. 2019 47:D298-D308, (PMID: 30418626).
#' @references Ullah et al. Sci. Rep. 2016 6:23534, (PMID: 27010073).
#' @references Durek et al. Nucleic Acids Res.2010 38:D828-D834, (PMID: 19880383).
#' @references Dinkel et al. Nucleic Acids Res. 2011 39:D261-D567 (PMID: 21062810).
#' @seealso meto.scan(), ac.scan(), me.scan(), ub.scan(), su.scan(), gl.scan(), sni.scan(), ni.scan(), p.scan(), reg.scan(), dis.scan()
#' @export

ptm.scan <- function(up_id, renumerate = TRUE){
  seq <- tryCatch(
    {
      get.seq(up_id, as.string = FALSE)
    },
    error = function(cond){
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(seq)){
    message("Sorry, get.seq failed")
    return(NULL)
  } else {
    seq <- seq[[1]]
  }

  output <- as.data.frame(matrix(rep(NA, length(seq)*14), ncol = 14))
  names(output) <- c('id','n', 'aa', 'meto', 'p', 'ac', 'me', 'ub', 'su', 'gl', 'sni', 'ni', 'reg', 'dis')
  output$id <- up_id
  output$n <- 1:nrow(output)
  output$aa <- seq

  ## ----- Sulfoxidation -------- ##
  meto <- tryCatch(
    {
      meto.scan(up_id)
    },
    error = function(cond){
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(meto)){
    message("Sorry, meto.scan failed")
    return(NULL)
  } else {
    meto <- meto[[1]]
  }

  if (!is.null(nrow(meto))){ # if there some meto site
    for (i in 1:nrow(meto)){
      t <- as.numeric(meto$met_pos[i])
      if (renumerate){
        t <- tryCatch(
          {
            renum(up_id, t, from = 'metosite', to = 'uniprot')
          },
          error = function(cond){
            return(NA)
          },
          warning = function(w) conditionMessage(w)
        )
      }
      output$meto[t] <- TRUE
    }
  }


  ## ----- Phosphorylation -------- ##
  p <- p.scan(up_id)
  if (!is.null(p)){
    if (nrow(p) != 0){ # if there some site
          u <- unique(p$modification)
          for (i in 1:length(u)){
            t <- strsplit(u[i], split = "-")[[1]][-2]
            t <- as.numeric(substring(t, 2))
            output$p[t] <- TRUE
          }
        }
  }

  ## ----- Acetylation -------- ##
  ac <- ac.scan(up_id)

  if (!is.null(ac)){
    if (nrow(ac) != 0){ # if there some site
      u <- unique(ac$modification)
      for (i in 1:length(u)){
        t <- strsplit(u[i], split = "-")[[1]][-2]
        t <- as.numeric(substring(t, 2))
        output$ac[t] <- TRUE
      }
    }
  }


  ## ----- Methylation -------- ##
  me <- me.scan(up_id)
  if (!is.null(me)){
    if (nrow(me) != 0){ # if there is some site
      u <- unique(me$modification)
      for (i in 1:length(u)){
        t <- strsplit(u[i], split = "-")[[1]][-2]
        t <- as.numeric(substring(t, 2))
        output$me[t] <- TRUE
      }
    }
  }


  ## ----- Ubiquitination -------- ##
  ub <- ub.scan(up_id)
  if (!is.null(ub)){
    if (nrow(ub) != 0){ # if there is some site
      u <- unique(ub$modification)
      for (i in 1:length(u)){
        t <- strsplit(u[i], split = "-")[[1]][-2]
        t <- as.numeric(substring(t, 2))
        output$ub[t] <- TRUE
      }
    }
  }


  ## ----- Sumoylation -------- ##
  su <- su.scan(up_id)
  if (!is.null(su)){
    if (nrow(su) != 0){ # if there is some site
      u <- unique(su$modification)
      for (i in 1:length(u)){
        t <- strsplit(u[i], split = "-")[[1]][-2]
        t <- as.numeric(substring(t, 2))
        output$su[t] <- TRUE
      }
    }
  }


  ## ----- OGlcNAc -------- ##
  gl <- gl.scan(up_id)
  if (!is.null(gl)){
    if (nrow(gl) != 0){ # if there is some site
      u <- unique(gl$modification)
      for (i in 1:length(u)){
        t <- strsplit(u[i], split = "-")[[1]][-2]
        t <- as.numeric(substring(t, 2))
        output$gl[t] <- TRUE
      }
    }
  }


  ## ----- S-nitrosylation -------- ##
  sni <- sni.scan(up_id)
  if (!is.null(sni)){
    if (nrow(sni) != 0){ # if there is some site
      u <- unique(sni$modification)
      for (i in 1:length(u)){
        t <- strsplit(u[i], split = "-")[[1]][-2]
        t <- as.numeric(substring(t, 2))
        output$sni[t] <- TRUE
      }
    }
  }


  ## ----- Nitration -------- ##
  ni <- ni.scan(up_id)
  if (!is.null(ni)){
    if (nrow(ni) != 0){ # if there is some site
      u <- unique(ni$modification)
      for (i in 1:length(u)){
        t <- strsplit(u[i], split = "-")[[1]][-2]
        t <- as.numeric(substring(t, 2))
        output$ni[t] <- TRUE
      }
    }
  }


  ## ----- Regulation -------- ##
  reg <- reg.scan(up_id)
  if (!is.null(reg)){
    if (nrow(reg) != 0){ # if there is some site
      u <- unique(reg$modification)
      for (i in 1:length(u)){
        t <- strsplit(u[i], split = "-")[[1]][-2]
        aa <- substr(t, 1, 1)
        t <- as.numeric(substring(t, 2))
        if (aa == 'M' & renumerate){
          t <- renum(up_id, t, from = 'metosite', to = 'uniprot')
        }
        output$reg[t] <- TRUE
      }
    }
  }


  ## ----- Disease -------- ##
  dis <- dis.scan(up_id)
  if (!is.null(dis)){
    if (nrow(dis) != 0){ # if there is any site
      u <- unique(dis$modification)
      for (i in 1:length(u)){
        t <- strsplit(u[i], split = "-")[[1]][-2]
        t <- as.numeric(substring(t, 2))
        output$dis[t] <- TRUE
      }
    }
  }

  output <- output[rowSums(is.na(output[, 4:14])) != 11,]

  if (nrow(output) == 0){
    message("Sorry, no modification sites were found for this protein")
    output <- NULL
  } else {
    o <- as.matrix(output)
    output$multi <- NA
    for (i in 1:nrow(o)){
      # output$multi[i] <- sum(as.logical(o[i,3:10]), na.rm = TRUE)
      output$multi[i] <- sum(as.logical(o[i,4:12]), na.rm = TRUE)
    }
  }

  if (!is.null(output)){
    attr(output, 'prot_id') <- up_id
    attr(output, 'prot_length') <- length(seq)
  }

  return(output)
}


## ---------------------------------------------------------------- ##
#                       reg.scan <- function(up_id)                  #
## ---------------------------------------------------------------- ##
#' Scan a Protein in Search of Regulatory PTM Sites
#' @description Scans the indicated protein in search of regulatory PTM sites.
#' @usage reg.scan(up_id)
#' @param up_id a character string corresponding to the UniProt ID.
#' @return Returns a dataframe where each row corresponds to a residue, and the columns inform about the regulatory modifications.
#' @author Juan Carlos Aledo
#' @examples \dontrun{reg.scan('P01009')}
#' @references Hornbeck et al. Nucleic Acids Res. 2019 47:D433-D441, (PMID: 30445427).
#' @seealso  meto.scan(), ac.scan(), me.scan(), ub.scan(), su.scan(), gl.scan(), sni.scan(), ni.scan(), ptm.scan(), p.scan(), dis.scan()
#' @export

reg.scan <- function(up_id){

  reg_db <- NULL
  baseUrl <- "https://github.com/jcaledo/reg_db/blob/master/"
  call <- paste(baseUrl, "reg_db_", up_id, ".Rda?raw=true", sep = "")

  tryCatch(
    {
      load(url(call))
    },
    error = function(cond){
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  closeAllConnections()

  if (is.null(reg_db)){
    message("Sorry, no regulatory sites found for this protein")
    return(NULL)
  } else {
    t <-  reg_db[which(reg_db$up_id == up_id),]
  }

  m <- meto.scan(up_id, report = 2)

  if (length(m$Metosites) != 0){
    tt <- m$Metosites[which(m$Metosites$reg_id > 2), ]

    if (nrow(tt) > 0){
      meto <- as.data.frame(matrix(rep(NA, 4*nrow(tt)), ncol = 4))
      names(meto) <- c('up_id','organism', 'modification', 'database')
      for (i in 1:nrow(tt)){
        meto$up_id[i] <- up_id
        meto$organism[i] <- m$prot_sp
        meto$modification[i] <- paste("M", tt$met_pos[i], "-ox", sep = "")
        meto$database[i] <- 'MetOSite'
      }
      if (nrow(t) > 0){
        t <- rbind(t, meto)
      } else {
        t <- meto
      }
    }
  }

  if (nrow(t) > 0){
    return(t)
  } else {
    message("Sorry, no modification sites were found for this protein")
    return(NULL)
  }
}

## ---------------------------------------------------------------- ##
#                     dis.scan <- function(up_id)                    #
## ---------------------------------------------------------------- ##
#' Scan a Protein in Search of Disease-Related PTM Sites
#' @description Scans the indicated protein in search of disease-related PTM sites.
#' @usage dis.scan(up_id)
#' @param up_id a character string corresponding to the UniProt ID.
#' @return Returns a dataframe where each row corresponds to a residue, and the columns inform about the disease-related modifications.
#' @author Juan Carlos Aledo
#' @examples \dontrun{dis.scan('P31749')}
#' @references Hornbeck et al. Nucleic Acids Res. 2019 47:D433-D441, (PMID: 30445427).
#' @seealso  meto.scan(), ac.scan(), me.scan(), ub.scan(), su.scan(), gl.scan(), sni.scan(), ni.scan(), ptm.scan(), reg.scan(), p.scan()
#' @export

dis.scan <- function(up_id){

  dis_db <- NULL
  baseUrl <- "https://github.com/jcaledo/dis_db/blob/master/"
  call <- paste(baseUrl, "dis_db_", up_id, ".Rda?raw=true", sep = "")

  tryCatch(
    {
      load(url(call))
    },
    error = function(cond){
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  closeAllConnections()

  if (is.null(dis_db)){
    message("Sorry, no modification sites were found for this protein")
    return(NULL)
  } else {
    t <-  dis_db[which(dis_db$up_id == up_id),]
    return(t)
  }
}
