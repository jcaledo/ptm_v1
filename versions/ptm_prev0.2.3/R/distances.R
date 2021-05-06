## ------- distances.R ------------ ##
#                                    #
#   pairwise.dist                    #
#   res.dist                         #
#   dist2closest                     #
#   ball                             #
#                                    #
## -------------------------------- ##

## ---------------------------------------------------------------- ##
#         pairwise.dist <- function(a, b, squared = TRUE)            #
## ---------------------------------------------------------------- ##
#' Compute Euclidean Distances
#' @description Computes the pairwise distance matrix between two sets of points
#' @usage pairwise.dist(a, b, squared = TRUE)
#' @param a,b matrices (NxD) and (MxD), respectively, where each row represents a D-dimensional point.
#' @param squared return containing squared Euclidean distance
#' @return Euclidean distance matrix (NxM). An attribute "squared" set to the
#' value of param \code{squared} is provided.
#' @examples pairwise.dist(matrix(1:9, ncol = 3), matrix(9:1, ncol = 3))
#' @seealso res.dist(), dist2closest(), ball()
#' @export

pairwise.dist <- function(a, b, squared = TRUE){
  an <- apply(a, 1, function(x) crossprod(x, x))
  bn <- apply(b, 1, function(x) crossprod(x, x))

  m <- length(an)
  n <- length(bn)

  an_bn <- matrix(rep(an, n), nrow=m) + matrix(rep(bn, m), nrow=m, byrow=T)
  d2 <- an_bn - 2 * tcrossprod(a, b)

  if(!squared){
    d2 <- sqrt(d2)
  }
  attr(d2, "squared") <- squared
  return(d2)
}

## ----------------------------------------------------------------- ##
# res.dist <-function(pdb, rA, chainA, rB, chainB, backbone, hatoms)  #
## ----------------------------------------------------------------- ##
#' Compute Distances Between Residues
#' @description Computes the euclidean distance between two given residues
#' @usage res.dist(pdb,  rA, chainA, rB, chainB, backbone = FALSE, hatoms = FALSE)
#' @param pdb is either a PDB id, or the path to a pdb file.
#' @param rA an integer indicating the position of the first residue.
#' @param chainA a character indicating the chain to which belong the first residue.
#' @param rB an integer indicating the position of the second residue.
#' @param chainB a character indicating the chain to which belong the second residue.
#' @param backbone logical, when TRUE it means that we include those atoms belonging to the main chain (CA, N, O and C) beside all the side chain atoms.
#' @param hatoms logical, if TRUE we include all the hydrogen atoms in the computation as long as the PDB provides their coordinates.
#' @return This function returns a list of three elements, where each of these elements is, in turn, a list of three elements providing information regarding minimal, maximal and averaged distances.
#' @author Juan Carlos Aledo
#' @examples \dontrun{res.dist('1q8k', 51, 'A', 55, 'A', backbone = TRUE, hatoms = TRUE)}
#' @seealso pairwise.dist(), dist2closest(), ball()
#' @importFrom bio3d read.pdb
#' @importFrom stats dist
#' @export

res.dist <- function(pdb, rA, chainA, rB, chainB, backbone = FALSE, hatoms = FALSE){

  # ------------ Data Preparation and Computation of Distances
  mypdb <- suppressWarnings(bio3d::read.pdb(pdb)) # avoids warning: 'pdb exists. Skipping download'

  atoA <- mypdb$atom[which(mypdb$atom$resno == rA & mypdb$atom$chain == chainA),]
  rownames(atoA) <- paste(atoA$elety, atoA$resid, atoA$resno, atoA$chain, sep = "-")
  rA <- paste(atoA$resid[1], atoA$resno[1], atoA$chain[1], sep = "-")
  atoA <- atoA[, c(which(colnames(atoA) == 'x' |
                           colnames(atoA) == 'y' |
                           colnames(atoA) == 'z'))]

  atoB <- mypdb$atom[which(mypdb$atom$resno == rB & mypdb$atom$chain == chainB),]
  rownames(atoB) <- paste(atoB$elety, atoB$resid, atoB$resno, atoB$chain, sep = "-")
  rB <- paste(atoB$resid[1], atoB$resno[1], atoB$chain[1], sep = "-")
  atoB <- atoB[, c(which(colnames(atoB) == 'x' |
                           colnames(atoB) == 'y' |
                           colnames(atoB) == 'z'))]
  atoms <- rbind(atoA, atoB)
  d <- as.matrix(round(dist(atoms), 2))
  d <- d[1:nrow(atoA), (nrow(atoA)+1):nrow(d)]

  #--------------- Main/Side Chain
  if (backbone == FALSE){
    # Atoms belonging to the backbone are removed
    residueA <- strsplit(rA, split = '-')[[1]][1]
    residueB <- strsplit(rB, split = '-')[[1]][1]
    if (residueA == 'GLY' & residueB != 'GLY'){ # Fist is Gly
      backboneA <- grep('N-[A-Z]|C-[A-Z]|O-[A-Z]', rownames(d))
      backboneB <- grep('CA-[A-Z]|N-[A-Z]|C-[A-Z]|O-[A-Z]', colnames(d))
    } else if (residueA != 'GLY' & residueB == 'GLY'){ # Secon is Gly
      backboneA <- grep('CA-[A-Z]|N-[A-Z]|C-[A-Z]|O-[A-Z]', rownames(d))
      backboneB <- grep('N-[A-Z]|C-[A-Z]|O-[A-Z]', colnames(d))
    } else if (residueA == 'GLY' & residueB == 'GLY'){ # Both are Gly
      backboneA <- grep('N-[A-Z]|C-[A-Z]|O-[A-Z]', rownames(d))
      backboneB <- grep('N-[A-Z]|C-[A-Z]|O-[A-Z]', colnames(d))
    } else { # Any of them is Gly
      backboneA <- grep('CA-[A-Z]|N-[A-Z]|C-[A-Z]|O-[A-Z]', rownames(d))
      backboneB <- grep('CA-[A-Z]|N-[A-Z]|C-[A-Z]|O-[A-Z]', colnames(d))
    }

    # When one (or both) of the residues is Ala, the
    # distance matrix (d) boils down to a vector and
    # that can cause problems with the row/column names.…
    if (dim(d)[1] == 5 & dim(d)[2] == 5){ # Both are Ala
      d <- d[-c(backboneA), -c(backboneB)]
      d <- as.matrix(d)
      rownames(d) <- paste("CB-", rA, sep = "")
      colnames(d) <- paste("CB-", rB, sep = "")
    } else if (dim(d)[1] == 5 & dim(d)[2] > 5){ # Fist is Ala
      d <- d[-c(backboneA), -c(backboneB)]
      d <- t(as.matrix(d))
      rownames(d) <- paste("CB-", rA, sep = "")
    } else if (dim(d)[1] > 5 & dim(d)[2] == 5){# Second is Ala
      d <- d[-c(backboneA), -c(backboneB)]
      d <- as.matrix(d)
      colnames(d) <- paste("CB-", rB, sep = "")
    } else { # Any is Ala
      d <- d[-c(backboneA), -c(backboneB)]
    }
  }

  #-------------- Hydrogen Atoms
  # Sometimes the pdb contains hydrogen atoms' coordinates
  if (hatoms == FALSE){
    filas <- rownames(d)[which(substr(rownames(d), 1,1) != "H")]
    columnas <- colnames(d)[which(substr(colnames(d), 1,1) != "H")]
    d <- as.matrix(d[which(rownames(d) %in%  filas),
                     which(colnames(d) %in% columnas)])
    if (dim(d)[1] == length(filas)){
      rownames(d) <- filas
      colnames(d) <- columnas
    } else if (dim(d)[2] == length(filas)){
      d <- t(d)
      rownames(d) <- filas
      colnames(d) <- columnas
    } else {
      stop("Problem with the distance matrix dimension")
    }
  }

  #--------------- Selecting Min, Max and Mean Distances
  min.d <- min(d)
  fila.min <- which(d == min(d), arr.ind = TRUE)[1]
  columna.min <- which(d == min(d), arr.ind = TRUE)[2]
  A <- rownames(d)[fila.min]
  B <- colnames(d)[columna.min]

  mylist <- list(list(min.d, A, B, "minimal distance"))

  max.d <- max(d)
  fila.max <- which(d == max(d), arr.ind = TRUE)[1]
  columna.max <- which(d == max(d), arr.ind=TRUE)[2]
  A <- rownames(d)[fila.max]
  B <- colnames(d)[columna.max]

  mylist[[2]] <- list(max.d, A, B, 'maximum distance')
  mylist[[3]] <- list(round(mean(d), 2), rA, rB, 'average distance')

  attr(mylist, "backbone") <- backbone
  return(mylist)
}


## ---------------------------------------------------------------- ##
#          dist2closest <- function(pdb, res, aa, backbone)          #
## ---------------------------------------------------------------- ##
#' Search and Compute the Distance to the Closest Aa
#' @description Computes the distance to the closest amino acid of the indicated type
#' @usage dist2closest(pdb, chain, res, aa = "M", backbone = FALSE)
#' @param pdb is either a PDB id, or the path to a pdb file.
#' @param chain chain ID.
#' @param res position of the residue of interest.
#' @param aa amino acid of interest.
#' @param backbone logical, when TRUE it means that we include those atoms belonging to the main chain (CA, N, O and C) beside all the side chain atoms.
#' @details The identity of the closest Aa is given as an attribute.
#' @return Numerical value indicating the distance in ångströms (Å).
#' @examples \dontrun{dist2closest(pdb = '1q8k', chain = 'A', res = 222, aa = 'S')}
#' @seealso res.dist(), pairwise.dist(), ball()
#' @importFrom bio3d get.pdb
#' @importFrom bio3d read.pdb
#' @importFrom bio3d aa123
#' @export

dist2closest <- function(pdb, chain, res, aa = 'M', backbone = FALSE){

  del <- FALSE
  if (nchar(pdb) == 4){ # when input is a PDB ID
    mypdb <- bio3d::get.pdb(pdb)
    pdb <- paste("./", pdb, ".pdb", sep = "")
    del <- TRUE
  }

  mypdb <- bio3d::read.pdb(pdb)$atom
  mypdb <- mypdb[which(mypdb$elety == 'CA'),]

  aa.pos <- mypdb$resno[which(mypdb$resid == aa123(aa))]
  aa.ch <- mypdb$chain[which(mypdb$resid == aa123(aa))]

  if (length(aa.pos) == 0){
    return(paste("Sorry, ", aa, " residues have not been found in this protein", sep = ""))
  } else if (length(aa.pos) == 1 & res %in% aa.pos){
    attr(m, 'closest residue') <- 'itself'
    return(0)
  }

  rA <- res
  chainA <- chain

  distances <- c()
  targets <- c()
  for (i in 1:length(aa.pos)){
    rB <- aa.pos[i]
    chainB <- aa.ch[i]
    t <- res.dist(pdb, rA, chainA, rB, chainB, backbone = backbone)
    distances <- c(distances, t[[1]][[1]])
    targets <- c(targets, paste(t[[1]][[3]], "<-->", t[[1]][[2]]))
  }

  m <- min(distances[which(distances != 0)]) # the min distance
  t <- targets[which(distances == min(distances[distances != 0]))] # the amino acid

  if (del){
    file.remove(pdb)
  }

  attr(m, 'closest residue') <- t

  return(m)
}


## ---------------------------------------------------------------- ##
#          ball <- function(pdb, chain, res, r, backbone)            #
## ---------------------------------------------------------------- ##
#' Search for Atoms Close to a Given Atom
#' @description Finds the atoms within a sphere with the indicated center and radius
#' @usage ball(pdb, chain, res, r, backbone = FALSE)
#' @param pdb is either a PDB id, or the path to a pdb file.
#' @param chain a character indicating the chain to which the residue belongs
#' @param res position in the primary structure of the residue of interest, which will be use as center of the sphere.
#' @param r radius in ångströms of the sphere.
#' @param backbone logical, when TRUE it means that we include those atoms belonging to the main chain (CA, N, O and C) beside all the side chain atoms.
#' @details The position indicated by res must be occupied by one of the following amino acids (otherwise the function will return an error message): Met (SD atom), Cys (SG atom), Glu (centroid of OE1 and OE2 atoms), Asp (centroid of OD1 and OD2 atoms), His (centroid of ND1 and NE2 atoms), Lys (NZ atom), Arg (centroid of NE, NH1 and NH2 atoms), Phe (centroid of CG, CD1, CD2, CE1, CE2 and CZ atoms), Tyr  (centroid of CG, CD1, CD2, CE1, CE2 and CZ atoms), Trp (centroid of the indol rings, ring-1: CG, CD1, CD2, NE1, CE2; ring-2: CD2, CE3, CZ2, CZ3, CE2). All the atoms belonging to the central residue are excluded from the results.
#' @return A dataframe with the atoms identification and their distances to the central atom.
#' @examples \dontrun{ball(pdb = '6e7f', chain = 'A', res = 181, r = 6, backbone = TRUE)}
#' @seealso res.dist(), pairwise.dist(), dist2closest()
#' @importFrom bio3d get.pdb
#' @importFrom bio3d read.pdb
#' @importFrom bio3d aa123
#' @export

ball <- function(pdb, chain, res, r, backbone = FALSE){

  mypdb <- suppressWarnings(bio3d::read.pdb(pdb)$atom) # avoids warning: 'pdb exists. Skipping download'

  ## --------------- Atoms to be used in the centroid computation
  residue <- mypdb$resid[which(mypdb$resno == res & mypdb$chain == chain)][1]
  if (residue == "MET"){
    atom <- mypdb[which(mypdb$elety == "SD" & mypdb$resno == res & mypdb$chain == chain), ]
  } else if (residue == "CYS"){
    atom <- mypdb[which(mypdb$elety == "SG" & mypdb$resno == res & mypdb$chain == chain), ]
  } else if (residue == "GLU"){
    atom <- mypdb[which(mypdb$elety %in% c('OE1', 'OE2') & mypdb$resno == res & mypdb$chain == chain), ]
  } else if (residue == "ASP"){
    atom <- mypdb[which(mypdb$elety %in% c('OD1', 'OD2') & mypdb$resno == res & mypdb$chain == chain), ]
  } else if (residue == "HIS"){
    atom <- mypdb[which(mypdb$elety %in% c('ND1', 'NE2') & mypdb$resno == res & mypdb$chain == chain), ]
  } else if (residue == "ARG"){
    Gdm<- c('NE', 'NH1', 'NH2')
    atom <- mypdb[which(mypdb$elety %in% Gdm & mypdb$resno == res & mypdb$chain == chain), ]
  } else if (residue == "LYS"){
    atom <- mypdb[which(mypdb$elety == "NZ" & mypdb$resno == res & mypdb$chain == chain), ]
  } else if (residue == "PHE" | residue == "TYR"){
    ben <- c('CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ')
    atom <- mypdb[which(mypdb$elety %in% ben & mypdb$resno == res & mypdb$chain == chain), ]
  } else if (residue == "TRP"){
    ring1 <- c('CG', 'CD1', 'CD2', 'NE1', 'CE2')
    ring2 <- c('CD2', 'CE3', 'CZ2', 'CH2', 'CZ3', 'CE2')
    atom1 <- mypdb[which(mypdb$elety %in% ring1 & mypdb$resno == res & mypdb$chain == chain), ]
    atom2 <- mypdb[which(mypdb$elety %in% ring2& mypdb$resno == res & mypdb$chain == chain), ]
  } else if (residue == "SER"){
    atom <- mypdb[which(mypdb$elety == "OG" & mypdb$resno == res & mypdb$chain == chain), ]
  } else if (residue == "THR"){
    atom <- mypdb[which(mypdb$elety == "OG1" & mypdb$resno == res & mypdb$chain == chain), ]
  } else {
    stop(paste("The residue ", residue, " is not a valid one!", sep =""))
  }

  ## ------------- Centroid computation
  if (residue == 'TRP'){
    centroid <- matrix(c(mean(atom1$x), mean(atom1$y), mean(atom1$z),
                         mean(atom2$x), mean(atom2$y), mean(atom2$z)), ncol = 3, byrow = TRUE)
  } else {
    centroid <- matrix(c(mean(atom$x), mean(atom$y), mean(atom$z)), ncol = 3)
  }

  ## -------------- Main/Side Chain
  r_mypdb <- mypdb[which(mypdb$resno != res), ] # atoms from the residue of interest are removed.
  if (backbone == FALSE){
    # Atoms belonging to the backbone are removed
    backbone_atoms <- grep('^CA$|^N$|^C$|^O$', r_mypdb$elety)
    r_mypdb <- r_mypdb[-backbone_atoms, ]
  }

  ## -------------- Matrices of points for pairwise.dist
  a <- matrix(c(r_mypdb$x, r_mypdb$y, r_mypdb$z), ncol = 3, byrow = FALSE)
  b <- centroid

  ## ------------ Computing pairwise distances
  output <- r_mypdb[ , colnames(r_mypdb) %in% c('elety', 'resid', 'chain',
                                                'resno', 'x', 'y', 'z', 'distance')]

  d <- pairwise.dist(a, b, squared = FALSE)
  if (residue == "TRP"){
    output$distance1 <- d[ , 1]
    output$distance2 <- d[ , 2]
    output <- output[which(output$distance1 <= r | output$distance2 <= r ), ]
  } else {
    output$distance <- d
    output <- output[which(output$distance <= r), ]
  }

  attr(output, 'structure') <- pdb
  attr(output, 'residue') <- residue
  attr(output, 'position') <- res
  attr(output, 'coordinates') <- centroid

  return(output)
}
