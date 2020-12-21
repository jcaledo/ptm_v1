## --------------- accdpx.R ----------------- ##
#                                              #
#       acc.dssp                               #
#       get.area                               #
#       dpx                                    #
#       atom.dpx                               #
#       res.dpx                                #
#       stru.part                              #
#                                              #
## ------------------------------------------ ##

## ---------------------------------------------------------------- ##
#     acc.dssp <- function(pdb, dssp = 'compute', aa = 'all')        #
## ---------------------------------------------------------------- ##
#' Compute Residue Accessibility and SASA
#' @description Computes the accessibility as well as the SASA for each reside from the indicated protein.
#' @usage acc.dssp(pdb, dssp = 'compute', aa = 'all')
#' @param pdb is either a PDB id, or the path to a pdb file.
#' @param dssp string indicating the preferred method to obtain the dssp file. It must be either 'compute' or 'mkdssp'.
#' @param aa one letter code for the amino acid of interest, or 'all' for all the protein residues.
#' @details For the given PDB the function obtains its corresponding DSSP file using the chosen method. The argument dssp allows two alternative methods:
#' 'compute' (it calls to the function compute.dssp(), which in turn uses an API to run DSSP at the CMBI);
#' 'mkdssp' (if you have installed DSSP on your system and in the search path for executables).
#' @return A dataframe where each row is an individual residue of the selected protein. The variables computed, among others, are:
#' (i) the secondary structure (ss) element to which the residue belongs,
#' (ii) the solvent accessible surface area (sasa) of each residue in square angstrom (Å²), and
#' (iii) the accessibility (acc) computed as the percent of the sasa that the residue X would have in the tripeptide GXG with the polypeptide skeleton in an extended conformation and the side chain in the conformation most frequently observed in proteins.
#' @author Juan Carlos Aledo
#' @examples \dontrun{acc.dssp('3cwm')
#' acc.dssp(pdb = '3cwm', aa = 'M')}
#' @references Miller et al (1987) J. Mol. Biol. 196: 641-656 (PMID: 3681970).
#' @references Touw et al (2015) Nucl. Ac. Res. 43(Database issue): D364-D368 (PMID: 25352545).
#' @seealso compute.dssp(), atom.dpx(), res.dpx(), str.part()
#' @importFrom bio3d get.pdb
#' @importFrom bio3d pdbsplit
#' @importFrom bio3d read.pdb
#' @export

acc.dssp <- function(pdb, dssp = 'compute', aa = 'all'){

  ## --------------- Download and split the PDB -------------- ##
  del <- FALSE
  if (nchar(pdb) == 4){
    del <- TRUE
    suppressWarnings(bio3d::get.pdb(pdb)) # avoids warning: 'pdb exists. Skipping download'
    file <- paste("./", pdb, ".pdb", sep = "")
    id <- pdb
  } else {
    file <- pdb
    t <- strsplit(file, split = "/")[[1]]
    id <- strsplit(t[length(t)], split = "\\.")[[1]][1]
  }
  chains <- pdb.chain(file, keepfiles = TRUE)

  ## ------------ Get the whole protein dssp file ------------ ##
  if (dssp == 'compute'){
    compute.dssp(pdb)
    dssp_file <- paste('./', id, '.dssp', sep = "")
    df <- parse.dssp(dssp_file, keepfiles = FALSE)

  } else if (dssp == 'mkdssp'){
    df <- mkdssp(id, method = 'ptm')
    # system(paste("mkdssp -i ", file, " -o ", dssp_file, sep =""))
  } else {
    stop("A proper dssp method should be provided!")
  }

  ## -------- Starting the dataframe construction ------------ ##
  df <- df[,1:6] # keep only relevant variables
  colnames(df)[6] <- 'sasa_complex'
  df$sasa_chain <- NA
  df$delta_sasa <- NA
  df$acc_complex <- NA
  df$acc_chain <- NA
  df$delta_acc <- NA

  ## --------- Get dssp files for individual chains  --------- ##
  chains <- unique(df$chain) # Sometimes we have to remove non-protein chains
  chain_files <- paste('./split_chain/', id, '_', chains,  '.pdb', sep = "")
  dssp_files <- paste('./split_chain/', id, '_', chains, '.dssp', sep = "")
  if (dssp == 'compute'){
    chain_files <- chain_files[file.exists(chain_files)]
    lapply(chain_files, function(x) compute.dssp(x, destfile = './split_chain/'))
    dssp_files <- dssp_files[file.exists(dssp_files)]
    df$sasa_chain <- unlist(lapply(dssp_files, function(x) parse.dssp(x, keepfiles = FALSE)$sasa))
  } else if (dssp == 'mkdssp'){
    df$sasa_chain <- unlist(lapply(chain_files, function(x) mkdssp(x, method = "ptm")$sasa))
  }
  df$delta_sasa <- df$sasa_chain - df$sasa_complex

  ## ---------------- Compute accessibility ------------------ ##
  df$acc_complex[which(df$aa == "A")] <-
    round(df$sasa_complex[which(df$aa == "A")]/113, 3)
  df$acc_chain[which(df$aa == "A")] <-
    round(df$sasa_chain[which(df$aa == "A")]/113, 3)

  df$acc_complex[which(df$aa == "R")] <-
    round(df$sasa_complex[which(df$aa == "R")]/241, 3)
  df$acc_chain[which(df$aa == "R")] <-
    round(df$sasa_chain[which(df$aa == "R")]/241, 3)

  df$acc_complex[which(df$aa == "N")] <-
    round(df$sasa_complex[which(df$aa == "N")]/158, 3)
  df$acc_chain[which(df$aa == "N")] <-
    round(df$sasa_chain[which(df$aa == "N")]/158, 3)

  df$acc_complex[which(df$aa == "D")] <-
    round(df$sasa_complex[which(df$aa == "D")]/151, 3)
  df$acc_chain[which(df$aa == "D")] <-
    round(df$sasa_chain[which(df$aa == "D")]/151, 3)

  df$acc_complex[which(df$aa == "C")] <-
    round(df$sasa_complex[which(df$aa == "C")]/140, 3)
  df$acc_chain[which(df$aa == "C")] <-
    round(df$sasa_chain[which(df$aa == "C")]/140, 3)
  # Because dssp convention, lower case letters (a, b, ...)
  # are introduced for bridged cysteines:
  df$acc_complex[which(df$aa %in% letters)] <-
    round(df$sasa_complex[which(df$aa %in% letters)]/140, 3)
  df$acc_chain[which(df$aa %in% letters)] <-
    round(df$sasa_chain[which(df$aa %in% letters)]/140, 3)

  df$acc_complex[which(df$aa == "Q")] <-
    round(df$sasa_complex[which(df$aa == "Q")]/189, 3)
  df$acc_chain[which(df$aa == "Q")] <-
    round(df$sasa_chain[which(df$aa == "Q")]/189, 3)

  df$acc_complex[which(df$aa == "E")] <-
    round(df$sasa_complex[which(df$aa == "E")]/183, 3)
  df$acc_chain[which(df$aa == "E")] <-
    round(df$sasa_chain[which(df$aa == "E")]/183, 3)

  df$acc_complex[which(df$aa == "G")] <-
    round(df$sasa_complex[which(df$aa == "G")]/85, 3)
  df$acc_chain[which(df$aa == "G")] <-
    round(df$sasa_chain[which(df$aa == "G")]/85, 3)

  df$acc_complex[which(df$aa == "H")] <-
    round(df$sasa_complex[which(df$aa == "H")]/194, 3)
  df$acc_chain[which(df$aa == "H")] <-
    round(df$sasa_chain[which(df$aa == "H")]/194, 3)

  df$acc_complex[which(df$aa == "I")] <-
    round(df$sasa_complex[which(df$aa == "I")]/182, 3)
  df$acc_chain[which(df$aa == "I")] <-
    round(df$sasa_chain[which(df$aa == "I")]/182, 3)

  df$acc_complex[which(df$aa == "L")] <-
    round(df$sasa_complex[which(df$aa == "L")]/180, 3)
  df$acc_chain[which(df$aa == "L")] <-
    round(df$sasa_chain[which(df$aa == "L")]/180, 3)

  df$acc_complex[which(df$aa == "K")] <-
    round(df$sasa_complex[which(df$aa == "K")]/211, 3)
  df$acc_chain[which(df$aa == "K")] <-
    round(df$sasa_chain[which(df$aa == "K")]/211, 3)

  df$acc_complex[which(df$aa == "M")] <-
    round(df$sasa_complex[which(df$aa == "M")]/204, 3)
  df$acc_chain[which(df$aa == "M")] <-
    round(df$sasa_chain[which(df$aa == "M")]/204, 3)

  df$acc_complex[which(df$aa == "F")] <-
    round(df$sasa_complex[which(df$aa == "F")]/218, 3)
  df$acc_chain[which(df$aa == "F")] <-
    round(df$sasa_chain[which(df$aa == "F")]/218, 3)

  df$acc_complex[which(df$aa == "P")] <-
    round(df$sasa_complex[which(df$aa == "P")]/143, 3)
  df$acc_chain[which(df$aa == "P")] <-
    round(df$sasa_chain[which(df$aa == "P")]/143, 3)

  df$acc_complex[which(df$aa == "S")] <-
    round(df$sasa_complex[which(df$aa == "S")]/122, 3)
  df$acc_chain[which(df$aa == "S")] <-
    round(df$sasa_chain[which(df$aa == "S")]/122, 3)

  df$acc_complex[which(df$aa == "T")] <-
    round(df$sasa_complex[which(df$aa == "T")]/146, 3)
  df$acc_chain[which(df$aa == "T")] <-
    round(df$sasa_chain[which(df$aa == "T")]/146, 3)

  df$acc_complex[which(df$aa == "W")] <-
    round(df$sasa_complex[which(df$aa == "W")]/259, 3)
  df$acc_chain[which(df$aa == "W")] <-
    round(df$sasa_chain[which(df$aa == "W")]/259, 3)

  df$acc_complex[which(df$aa == "Y")] <-
    round(df$sasa_complex[which(df$aa == "Y")]/229, 3)
  df$acc_chain[which(df$aa == "Y")] <-
    round(df$sasa_chain[which(df$aa == "Y")]/229, 3)

  df$acc_complex[which(df$aa == "V")] <-
    round(df$sasa_complex[which(df$aa == "V")]/160, 3)
  df$acc_chain[which(df$aa == "V")] <-
    round(df$sasa_chain[which(df$aa == "V")]/160, 3)

  df$delta_acc <- df$acc_chain - df$acc_complex

  ## ---------------- Cleaning and Tidying ------------------- ##
  unlink("./temp_split", recursive = TRUE)
  aai <- as.character(aai$aa)

  if (del & file.exists(file)){
    file.remove(file)
  }

  ## ------------------ Returning Results -------------------- ##
  if (aa != 'all' & aa %in% aai){
    df <- df[which(df$aa == aa),]
    return(df)
  } else if (aa == 'all'){
    return(df)
  } else {
    stop("A proper amino acid must be selected!")
  }
}


## ---------------------------------------------------------------- ##
#         get.area <- function(pdb, keepfiles = FALSE)               #
## ---------------------------------------------------------------- ##
#' Atomic Solvation Energies.
#' @description Computes online surface energies using the Getarea server.
#' @usage get.area(pdb, keepfiles = FALSE)
#' @param pdb is either a PDB id, or the path to a pdb file.
#' @param keepfiles logical, if TRUE the dataframe will be saved in the working directory and we will keep the getarea txt file.
#' @details If the option keepfiles is set as TRUE, then txt and Rda files are saved in the working directory.
#' @return This function returns a dataframe containing the requested information.
#' @author Juan Carlos Aledo
#' @examples \dontrun{get.area('3cwm')}
#' @references Fraczkiewicz, R. and Braun, W. (1998) J. Comp. Chem., 19, 319-333.
#' @seealso compute.dssp(), atom.dpx(), res.dpx(), acc.dssp(), str.part()
#' @importFrom bio3d aa321
#' @importFrom RCurl fileUpload
#' @importFrom RCurl postForm
#' @export

get.area <- function(pdb, keepfiles = FALSE){

  del <- FALSE
  if (nchar(pdb) == 4){
    del <- TRUE
    suppressWarnings(get.pdb(pdb)) # avoids warning: 'pdb exists. Skipping download'
    file <- paste("./", pdb, ".pdb", sep = "")
    id <- pdb
  } else {
    file <- pdb
    t <- strsplit(file, split = "/")[[1]]
    id <- strsplit(t[length(t)], split = "\\.")[[1]][1]
  }

  output_file = gsub("\\.pdb", "_getarea.txt", file)
  url <- "http://curie.utmb.edu/cgi-bin/getarea.cgi"

  result <-  postForm(url,
                     "water" = "1.4",
                     "gradient" = "n",
                     "name" = "test",
                     "email" = 'metosite2018@gmail.com',
                     "Method" = "4",
                     "PDBfile" = fileUpload(file))

  writeLines(gsub("[</pre></td>|<td><pre>]","", result), con = output_file)

  con <- file(output_file, 'r')
  counter <- 0
  eleno <- c()
  elety <- c()
  alt <- c()
  resid <- c()
  resno <- c()
  areaenergy <- c()

  while(TRUE){
    counter <- counter + 1
    line <- readLines(con, n = 1)

    if (counter > 17 & length(line) != 0){
      a <- strsplit(line, split = "")[[1]]
      eleno <- c(eleno, paste(a[1:6], collapse = ""))
      elety <- c(elety, paste(a[7:10], collapse = ""))
      alt <- c(alt, paste(a[11:12], collapse = ""))
      resid <- c(resid, paste(a[13:15], collapse = ""))
      resno <- c(resno, paste(a[16:23], collapse = ""))
      areaenergy <- c(areaenergy, paste(a[24:30], collapse = ""))
    } else if (length(line) == 0) {
      break
    }
  }
  close(con)

  ## ------ Setting the variable types ------------- ##
  eleno <- eleno[1:(length(eleno)-11)]
  elety <- gsub(" ", "", elety[1:(length(elety)-11)])
  alt <- alt[1:(length(alt)-11)]
  resid <- aa321(resid[1:(length(resid)-11)])
  resno <- gsub(" ", "", resno[1:(length(resno)-11)])
  areaenergy <- gsub(" ", "", areaenergy[1:(length(areaenergy)-11)])

  ## -------- Building the dataframe ---------------- ##
  df <- as.data.frame(matrix(c(alt, eleno, elety, resid, resno, areaenergy),
                             ncol = 6), stringsAsFactors = FALSE)

  colnames(df) <- c('alt','eleno', 'elety', 'resid', 'resno', 'areaenergy')

  df$alt <- gsub(" ", "", as.character(df$alt))
  df$eleno <- as.numeric(df$eleno)
  df$resno <- as.numeric(df$resno)
  df$areaenergy <- as.numeric(df$areaenergy)

  ## --- When PDB has ALT records, taking A only --- ##
  df <- df[which(df$alt != "B"),]
  df <- df[,-1]

  if (keepfiles == TRUE){
    rda_file <- gsub(".txt", ".Rda", output_file)
    save(df, file = rda_file)
  } else {
    file.remove(output_file)
    if (del){
      file.remove(file)
    }
  }
  return(df)
}

## ---------------------------------------------------------------- ##
#                   dpx <- function(pdb)                             #
## ---------------------------------------------------------------- ##
#' Atom Depth Analysis
#' @description Computes the depth from the surface for each protein's atom.
#' @usage dpx(pdb)
#' @param pdb is either a PDB id, or the path to a pdb file.
#' @details This function computes the depth, defined as the distance in angstroms between the target atom and the closest atom on the protein surface.
#' @return A dataframe with the computed depths.
#' @author Juan Carlos Aledo
#' @references Pintar et al. 2003. Bioinformatics 19:313-314 (PMID: 12538266)
#' @seealso compute.dssp(), atom.dpx(), res.dpx(), acc.dssp(), str.part()
#' @importFrom bio3d read.pdb
#' @export

dpx <- function(pdb){

  ## --------------- Generate a dataframe ---------------- ##
  atom <- get.area(pdb)
  mypdb <- suppressWarnings(read.pdb(pdb)) # avoids warning: 'pdb exists. Skipping download'

  mypdb <- mypdb$atom[which(mypdb$atom$type == 'ATOM'),]
  atom$x <- mypdb$x
  atom$y <- mypdb$y
  atom$z <- mypdb$z
  atom$dpx <- NA
  atom$dpx[which(atom$areaenergy != 0)] = 0

  exposed <- atom[which(atom$areaenergy != 0), c(6:8)]
  exposed <- as.matrix(exposed) # xyz coordinates of exposed atoms
  buried <- as.matrix(atom[is.na(atom$dpx), 6:8]) # xyz coordinates of buried atoms
  dpx <- round(pairwise.dist(buried, exposed, FALSE),2)
  dpx_buried <- apply(dpx, 1, min) # distance to the closest exposed atom
  atom$dpx[is.na(atom$dpx)] <- dpx_buried

  return(atom)
}


## ---------------------------------------------------------------- ##
#                atom.dpx <- function(pdb)                           #
## ---------------------------------------------------------------- ##
#' Atom Depth Analysis
#' @description Computes the depth from the surface for each protein's atom.
#' @usage atom.dpx(pdb)
#' @param pdb is either a PDB id, or the path to a pdb file.
#' @details This function computes the depth, defined as the distance in angstroms between the target atom and the closest atom on the protein surface. When the protein is composed of several subunits, the calculations are made for both, the atom being part of the complex, and the atom being only part of the polypeptide chain to which it belongs.
#' @return A dataframe with the computed depths.
#' @author Juan Carlos Aledo
#' @examples \dontrun{atom.dpx('1cll')}
#' @references Pintar et al. 2003. Bioinformatics 19:313-314 (PMID: 12538266)
#' @seealso res.dpx(), acc.dssp(), str.part()
#' @importFrom bio3d read.pdb
#' @export

atom.dpx <- function(pdb){
  del <- FALSE
  if (nchar(pdb) == 4){
    del <- TRUE
    id <- pdb
  } else {
    t <- strsplit(pdb, split = "/")[[1]]
    id <- strsplit(t[length(t)], split = "\\.")[[1]][1]
  }

  ## ------------ For the whole structure ---------------- ##
  atom <- dpx(pdb)
  energies <- atom$areaenergy
  dpx <- atom$dpx
  atom <- atom[,1:4]
  atom$chain <- NA
  atom$areaenergy <- energies
  atom$dpx_complex <- dpx
  atom$dpx_chain <- NA
  atom$delta_dpx <- NA

  ## ----------- For multi-chains structure -------------- ##
  chains <- suppressWarnings(pdb.chain(pdb, keepfiles = TRUE))
  counter_a <- 0
  for (i in seq_len(length(chains))){
    t <- paste('./split_chain/', id, '_', chains[i], '.pdb', sep = "")
    atom_chain <- dpx(t)
    counter_b <- counter_a + nrow(atom_chain)
    atom$chain[(counter_a + 1) : counter_b] <- chains[i]
    atom$dpx_chain[(counter_a + 1) : counter_b] <- atom_chain$dpx
    counter_a <- counter_b
  }
  atom$delta_dpx <- atom$dpx_complex - atom$dpx_chain

  ## ---------------- Cleaning and Tidying ------------------- ##
  unlink("./split_chain", recursive = TRUE)
  if (del){
    file.remove(paste("./", id, ".pdb", sep = ""))
  }

  return(atom)
}


## ---------------------------------------------------------------- ##
#             res.dpx <- function(pdb, aa = 'all')                   #
## ---------------------------------------------------------------- ##
#' Residue Depth Analysis
#' @description Computes the depth from the surface for each protein's residue.
#' @usage res.dpx(pdb, aa = 'all')
#' @param pdb is either a PDB id, or the path to a pdb file.
#' @param aa one letter code for the amino acid of interest, or 'all' for all the protein residues.
#' @details This function computes the depth, defined as the distance in angstroms between the residue and the closest atom on the protein surface.
#' @return A dataframe with the computed depths.
#' @author Juan Carlos Aledo
#' @examples \dontrun{res.dpx('1cll')}
#' @references Pintar et al. 2003. Bioinformatics 19:313-314 (PMID: 12538266)
#' @seealso atom.dpx(), acc.dssp(), str.part()
#' @export

res.dpx <- function(pdb, aa = 'all'){

  atom <- atom.dpx(pdb)
  atom$res <- paste(atom$resid, atom$resno, atom$chain, sep ="-")
  res <- unique(atom$res)

  ## ------------------ Dataframe structure -------------------- ##
  df <- as.data.frame(matrix(rep(NA, length(res) * 10), ncol = 10))
  colnames(df) <- c('res','resno', 'resid', 'chain', 'min_dpx_complex',
                    'avg_dpx_complex', 'max_dpx_complex',
                    'min_dpx_chain', 'avg_dpx_chain', 'max_dpx_chain')
  df$res <- res
  for (i in 1:nrow(df)){
    t <- strsplit(res[i], split = '-')[[1]]
    df$resno[i] <- as.numeric(t[2])
    df$resid[i] <- t[1]
    df$chain[i] <- t[3]
  }

  ## -------------------- Depth computation -------------------- ##
  df$min_dpx_complex <- df$min_dpx_chain <- NA
  df$max_dpx_complex <- df$max_dpx_chain <- NA
  df$avg_dpx_complex <- df$avg_dpx_chain <- NA

  for (i in 1:length(res)){
    temp <- atom[which(atom$res == res[i]),]

    df$min_dpx_chain[i] <- min(temp$dpx_chain)
    df$avg_dpx_chain[i] <- round(mean(temp$dpx_chain), 2)
    df$max_dpx_chain[i] <- max(temp$dpx_chain)

    df$min_dpx_complex[i] <- min(temp$dpx_complex)
    df$avg_dpx_complex[i] <- round(mean(temp$dpx_complex), 2)
    df$max_dpx_complex[i] <- max(temp$dpx_complex)
  }

  ## ------------------ Returning Results -------------------- ##
  if (aa != 'all' & aa %in% aai$aa){
    df <- df[which(df$resid == aa),]
    return(df)
  } else if (aa == 'all'){
    return(df)
  } else {
    stop("A proper amino acid must be provided!")
  }
}

## ---------------------------------------------------------------- ##
#           stru.part <- function(pdb, cutoff = 0.25)                #
## ---------------------------------------------------------------- ##
#' Partition of Structural Regions
#' @description Carries out a partition of the structural regions of a given protein.
#' @usage stru.part(pdb, cutoff = 0.25)
#' @param pdb is either a PDB id, or the path to a pdb file
#' @param cutoff accessibility below which a residue is considered to be buried.
#' @details The accessibilities of a residue computed in the complex (ACCc) and in the monomer (ACCm) allow to distinguish four structural regions as follows.
#'
#' Interior: ACCc < cutoff & (ACCm - ACCc) = 0.
#'
#' Surface: ACCc > cutoff &  (ACCm - ACCc) = 0.
#'
#' Support: ACCm < cutoff & (ACCm - ACCc) > 0.
#'
#' Rim: ACCc > cutoff & (ACCm - ACCc) > 0.
#'
#' Core: ACCm > cutoff & ACCc < cutoff.
#' @return A dataframe where each residue is assigned to one of the four structural groups considered.
#' @author Juan Carlos Aledo
#' @examples \dontrun{stru.part('1u8f')}
#' @references Levy (2010) J. Mol. Biol. 403: 660-670 (PMID: 20868694).
#' @seealso  atom.dpx(), res.dpx(), acc.dssp()
#' @export

stru.part <- function(pdb, cutoff = 0.25){

  ## --- Getting Accessibilities
  t <- acc.dssp(pdb, dssp = 'compute', aa = 'all')
  output <- t[, c(1:5, 9:11)]
  names(output) <- c(names(t)[1:5], "ACCc", "ACCm", "DACC")
  output$str <- NA

  ## --- Sorting Structural Regions
  for (i in 1:nrow(output)){
    if (output$DACC[i] == 0 & output$ACCc[i] < cutoff){
      output$str[i] <- "interior"
    } else if (output$DACC[i] == 0 & output$ACCc[i] > cutoff){
      output$str[i] <- 'surface'
    } else if (output$DACC[i] > 0 & output$ACCm[i] < cutoff){
      output$str[i] <- 'support'
    } else if (output$DACC[i] > 0 & output$ACCc[i] > cutoff){
      output$str[i] <- 'rim'
    } else {
      output$str[i] <- 'core'
    }
  }
  return(output)
}


