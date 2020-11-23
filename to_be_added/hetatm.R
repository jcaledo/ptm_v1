# hetatm.R
## ----------------------------------------------------------------- ##
#   search.hetatm <-function(file, destfile = './', rmH2O = TRUE)     #
## ----------------------------------------------------------------- ##
#' Search for Heteroatoms
#' @description Search for heteroatoms within a PDB file.
#' @usage search.hetatm(file, destfile = './', rmH2O = TRUE)
#' @param file a character string with the path to the PDB file of interst.
#' @param destfile a character string with the path where the split files will be saved.
#' @param rmH2O logical, when TRUE those atoms from water molecules will not be considered as heteroatoms.
#' @details The names of the heteroatoms are given as an attribute of the output object.
#' @return This function returns a dataframe with the minimal, maximal and averaged distances distances between the protein residues and the heteroatoms.
#' @author Juan Carlos Aledo
#' @examples search.hetatm('./1aaq.pdb')
#' @references Laurie & Jackson (2005) Bioinformatics 21:1908-1916
#' @seealso res.dist
#' @importFrom bio3d read.pdb
#' @importFrom bio3d write.pdb
#' @export

search.hetatm <- function(file, destfile = './', rmH2O = TRUE){

  all <- read.pdb(file)$atom
  atm <- all[which(all$type == 'ATOM'),]
  het <- all[which(all$type == 'HETATM'),]
  het_names <- unique(het$resid)

  if (rmH2O){
    het <- het[which(het$resid != 'HOH'),]
  }

  if (nrow(het) == 0){
    return("No heteroatoms found")
    stop
  }

  id <- stringr::str_extract(file, "\\w{4}\\.pdb")
  atm_id <- gsub("\\.", "_atm\\.", id)
  het_id <- gsub("\\.", "_het\\.", id)

  file_atm <- gsub("\\w{4}\\.pdb", atm_id, file)
  file_het <- gsub("\\w{4}\\.pdb", het_id, file)

  ## ----- Saving the split pdb files --------- ##
  atm_xyz <- c()
  for (i in 1:nrow(atm)){
    atm_xyz <- c(atm_xyz, atm$x[i], atm$y[i], atm$z[i])
  }

  write.pdb(file = file_atm, type = atm$type, eleno = atm$eleno,
            elety = atm$elety, resid = atm$resid, chain = atm$chain,
            resno = atm$resno, xyz = atm_xyz)


  het_xyz <- c()
  for (i in 1:nrow(het)){
    het_xyz <- c(het_xyz, het$x[i], het$y[i], het$z[i])
  }

  write.pdb(file = file_het, type = het$type, eleno = het$eleno,
            elety = het$elety, resid = het$resid, chain = het$chain,
            resno = het$resno, xyz = het_xyz)

  ## ------- Computing distances -------------- ##

  #-- Coordinates for the heteroatoms:
  het_xyz <- het[ , c(which(colnames(het) == 'x' |
                              colnames(het) == 'y' |
                              colnames(het) == 'z'))]
  # If chain id is absent substitute by eleno to avoid
  # repeated rows:
  het$chain[which(is.na(het$chain))] <- het$eleno
  rownames(het_xyz) <- paste(het$elety, het$resid,
                             het$resno, het$chain, sep = "-")

  #-- Regarding the apoprotein moiety (residue by residue):
  u <- unique(paste(atm$resid, atm$resno, atm$chain, sep = "-"))
  unico <- strsplit(u, split = "-")

  #-- Output dataframe's structure:
  output <- as.data.frame(matrix(rep(NA, length(unico)*4), ncol = 4))
  colnames(output) <- c('res','min', 'max', 'avg')
  output$res <- u
  #--
  for (i in 1:length(unico)){
    aa <- unico[[i]][1]
    at <- unico[[i]][2]
    chain <- unico[[i]][3]
    atm_t <- atm[which(atm$resid == aa & atm$resno == at & atm$chain == chain),]
    atm_txyz <- atm_t[ , c(which(colnames(atm_t) == 'x' |
                                 colnames(atm_t) == 'y' |
                                 colnames(atm_t) == 'z'))]
    res_het <- rbind(atm_txyz, het_xyz)
    d <- round(as.matrix(dist(res_het)), 2)
    d <- d[1:nrow(atm_t), (nrow(atm_t)+1):nrow(d)]
    output$min[i] <- min(d)
    output$max[i] <- max(d)
    output$avg[i] <- round(mean(d), 2)
  }
  attr(output, 'hetatm') <- het_names
  return(output)
}


## ----------------------------------------------------------------- ##
#   asdist.distribution <- function(file,  res = 'M')                 #
## ----------------------------------------------------------------- ##
#' Distance to the Active Site Distribution
#' @description Plots the distance (to the AS) distribution.
#' @usage asdist.distribution(file, res = 'M')
#' @param file a character string with the path to the PDB file of interst.
#' @param res a character indicating the amino acid to be studied.
#' @details Bla, bla, bla.
#' @return This function returns a plot of the distance to the AS for each residue emphasizing the points corresponding to the chosen residue.
#' @author Juan Carlos Aledo
#' @examples asdist.distribution('./1aaq.pdb')
#' @references bla, bla, bla
#' @seealso plot.ptm()
#' @importFrom bio3d read.pdb
#' @export

asdist.distribution <- function(file, res = 'M'){
  ## -------- Dataframe with distances to the Active Site (South Pole) ------- ##
  df <- search.hetatm(file)
  df$chain <- df$aa <- df$pos <- NA
  for (i in 1:nrow(df)){
    t <- strsplit(df$res[i], split = "-")[[1]]
    df$aa[i] <- bio3d::aa321(t[1])
    df$pos[i] <- t[2]
    df$chain[i] <- t[3]
    df$resno[i] <- i
  }
  ## ------------------------ Latitude --------------------------------------- ##
  max <- max(df$avg)
  min <- min(df$avg)
  equator <- round(mean(c(max, min)), 1)
  cancer <- round(mean(c(equator, max)), 1)
  capricorn <- round(mean(c(equator, min)), 1)

  df$latitude <- NA
  for (i in 1:nrow(df)){
    if (df$avg[i] >= cancer){
      df$latitude[i] <- 'artic'
    } else if (df$avg[i] < cancer & df$avg[i] >= equator){
      df$latitude[i] <- 'cancer'
    } else if (df$avg[i] < equator & df$avg[i] >= capricorn){
      df$latitude[i] <- 'capricorn'
    } else if (df$avg[i] < capricorn) {
      df$latitude[i] <- 'antartic'
    } else {
      stop("Something went awfully wrong!")
    }
  }

  df$latitude <- as.factor(df$latitude)
  # order (artic -> cancer -> capricorn -> antartic):
  df$latitude = factor(df$latitude, levels = levels(df$latitude)[c(2,3,4,1)])

  # Number or residues in the different regions:
  n_artic <- length(which(df$latitude == 'artic'))
  n_cancer <- length(which(df$latitude == 'cancer'))
  n_capricorn <- length(which(df$latitude == 'capricorn'))
  n_antartic <- length(which(df$latitude == 'antartic'))
  population <- paste("(", n_artic, ",", n_cancer, ",", n_capricorn, ",", n_antartic, ")", sep = "")

  ## ------------------------ Plotting Results --------------------------------- ##
  t <- strsplit(file, split = "\\/")[[1]]
  pdb_id <- strsplit(t[length(t)], split = "\\.")[[1]][1]
  title <- paste(pdb_id, res, population, sep = "-")
  dylab <- expression(paste("Distance to AS (", ring(A), ")", sep = ""))
  plot(1:nrow(df), df$avg, type = 'l', main = title, xlab = 'Residue Number', ylab = '')
  mtext(text = dylab, side = 2, line = 2)

  abline(equator, 0, lty = 2)
  abline(capricorn, 0, lty = 3)
  abline(cancer, 0, lty = 3)

  points(df$resno[which(df$aa == res & df$avg >= cancer)],
         df$avg[which(df$aa == res & df$avg >= cancer)], pch = 19, col = 'blue')

  points(df$resno[which(df$aa == res & df$avg < cancer & df$avg >= equator)],
         df$avg[which(df$aa == res & df$avg < cancer & df$avg >= equator)], pch = 19, col = 'green')

  points(df$resno[which(df$aa == res & df$avg < equator & df$avg >= capricorn)],
         df$avg[which(df$aa == res & df$avg < equator & df$avg >= capricorn)], pch = 19, col = 'orange')

  points(df$resno[which(df$aa == res & df$avg < capricorn)],
         df$avg[which(df$aa == res & df$avg < capricorn)], pch = 19, col = 'red')

}

