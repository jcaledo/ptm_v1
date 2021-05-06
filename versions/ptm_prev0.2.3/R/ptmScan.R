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

  resp <- try(load(url(call)), silent = TRUE)
  if (inherits(resp, "try-error")) {
    text <- "Sorry, no modification sites were found for this protein"
    return(text)
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

  resp <- try(load(url(call)), silent = TRUE)
  if (inherits(resp, "try-error")) {
    text <- "Sorry, no modification sites were found for this protein"
    return(text)
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
#' @examples me.scan('Q16695', db = 'PSP')
#' @seealso meto.scan(), ac.scan(), p.scan(), ub.scan(), su.scan(), gl.scan(), sni.scan(), ni.scan(), ptm.scan(), reg.scan(), dis.scan()
#' @export

me.scan <- function(up_id, db = 'all'){

  if (! db %in% c('PSP', 'dbPTM','dbPAF', 'PhosPhAt', 'Phospho.ELM', 'all')){
    stop("Please, select an appropiated database!")
  }

  me_db <- NULL
  baseUrl <- "https://github.com/jcaledo/me_db/blob/master/"
  call <- paste(baseUrl, "me_db_", up_id, ".Rda?raw=true", sep = "")

  resp <- try(load(url(call)), silent = TRUE)
  if (inherits(resp, "try-error")) {
    text <- "Sorry, no modification sites were found for this protein"
    return(text)
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
#' @examples ub.scan('Q16695', db = 'PSP')
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

  resp <- try(load(url(call)), silent = TRUE)
  if (inherits(resp, "try-error")) {
    text <- "Sorry, no modification sites were found for this protein"
    return(text)
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
#' @examples su.scan('Q16695', db = 'PSP')
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

  resp <- try(load(url(call)), silent = TRUE)
  if (inherits(resp, "try-error")) {
    text <- "Sorry, no modification sites were found for this protein"
    return(text)
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
#' @examples gl.scan('P08670', db = 'PSP')
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

  resp <- try(load(url(call)), silent = TRUE)
  if (inherits(resp, "try-error")) {
    text <- "Sorry, no modification sites were found for this protein"
    return(text)
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
#' @examples sni.scan('P01009')
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

  resp <- try(load(url(call)), silent = TRUE)
  if (inherits(resp, "try-error")) {
    text <- "Sorry, no modification sites were found for this protein"
    return(text)
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
#' @examples ni.scan('P05202')
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

  resp <- try(load(url(call)), silent = TRUE)
  if (inherits(resp, "try-error")) {
    text <- "Sorry, no modification sites were found for this protein"
    return(text)
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
  seq <- get.seq(up_id, as.string = FALSE)[[1]]
  output <- as.data.frame(matrix(rep(NA, length(seq)*14), ncol = 14))
  names(output) <- c('id','n', 'aa', 'meto', 'p', 'ac', 'me', 'ub', 'su', 'gl', 'sni', 'ni', 'reg', 'dis')
  output$id <- up_id
  output$n <- 1:nrow(output)
  output$aa <- seq

  ## ----- Sulfoxidation -------- ##
  meto <- meto.scan(up_id)[[1]]
  if (!is.null(nrow(meto))){ # if there is any meto site
    for (i in 1:nrow(meto)){
      t <- as.numeric(meto$met_pos[i])
      if (renumerate){
        t <- renum(up_id, t, from = 'metosite', to = 'uniprot')
      }
      output$meto[t] <- TRUE
    }
  }


  ## ----- Phosphorylation -------- ##
  p <- p.scan(up_id)
  if (!grepl("Sorry", p[1])){
    if (nrow(p) != 0){ # if there is any site
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
  if (!grepl("Sorry", ac[1])){
    if (nrow(ac) != 0){ # if there is any site
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
  if (!grepl("Sorry", me[1])){
    if (nrow(me) != 0){ # if there is any site
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
  if (!grepl("Sorry", ub[1])){
    if (nrow(ub) != 0){ # if there is any site
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
  if (!grepl("Sorry", su[1])){
    if (nrow(su) != 0){ # if there is any site
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
  if (!grepl("Sorry", gl[1])){
    if (nrow(gl) != 0){ # if there is any site
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
  if (!grepl("Sorry", sni[1])){
    if (nrow(sni) != 0){ # if there is any site
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
  if (!grepl("Sorry", ni[1])){
    if (nrow(ni) != 0){ # if there is any site
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
  if (!grepl("Sorry", reg[1])){
    if (nrow(reg) != 0){ # if there is any site
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
  if (!grepl("Sorry", dis[1])){
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
    output <- "Sorry, no modification sites were found for this protein"
  } else {
    o <- as.matrix(output)
    output$multi <- NA
    for (i in 1:nrow(o)){
      # output$multi[i] <- sum(as.logical(o[i,3:10]), na.rm = TRUE)
      output$multi[i] <- sum(as.logical(o[i,4:12]), na.rm = TRUE)
    }
  }

  attr(output, 'prot_id') <- up_id
  attr(output, 'prot_length') <- length(seq)

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
#' @examples reg.scan('P01009')
#' @references Hornbeck et al. Nucleic Acids Res. 2019 47:D433-D441, (PMID: 30445427).
#' @seealso  meto.scan(), ac.scan(), me.scan(), ub.scan(), su.scan(), gl.scan(), sni.scan(), ni.scan(), ptm.scan(), p.scan(), dis.scan()
#' @export

reg.scan <- function(up_id){

  reg_db <- NULL
  baseUrl <- "https://github.com/jcaledo/reg_db/blob/master/"
  call <- paste(baseUrl, "reg_db_", up_id, ".Rda?raw=true", sep = "")
  resp <- try(load(url(call)), silent = TRUE)
  if (inherits(resp, "try-error")) {
    t <- data.frame()
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
    text <- "Sorry, no modification sites were found for this protein"
    return(text)
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
#' @examples dis.scan('P31749')
#' @references Hornbeck et al. Nucleic Acids Res. 2019 47:D433-D441, (PMID: 30445427).
#' @seealso  meto.scan(), ac.scan(), me.scan(), ub.scan(), su.scan(), gl.scan(), sni.scan(), ni.scan(), ptm.scan(), reg.scan(), p.scan()
#' @export

dis.scan <- function(up_id){

  dis_db <- NULL
  baseUrl <- "https://github.com/jcaledo/dis_db/blob/master/"
  call <- paste(baseUrl, "dis_db_", up_id, ".Rda?raw=true", sep = "")
  resp <- try(load(url(call)), silent = TRUE)
  if (inherits(resp, "try-error")) {
    text <- "Sorry, no modification sites were found for this protein"
    return(text)
  }

  t <-  dis_db[which(dis_db$up_id == up_id),]
  return(t)
}
