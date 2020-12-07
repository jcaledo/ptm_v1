## ------------- aa.R ------------- ##
#                                    #
#     aa.at                          #
#     is.at                          #
#     aa.comp                        #
#     renum.pdb                      #
#     renum.meto                     #
#     renum                          #
#                                    #
## -------------------------------- ##

## --------------------------------------------------------------- ##
#             aa.at(at, target, uniprot = TRUE)                     #
## --------------------------------------------------------------- ##
#' Residue Found at the Requested Position
#' @description Returns the residue found at the requested position.
#' @usage aa.at(at, target, uniprot = TRUE)
#' @param at the position in the primary structure of the protein.
#' @param target a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param uniprot logical, if TRUE the argument 'target' should be an ID.
#' @details Please, note that when uniprot is set to FALSE, target can be the string returned by a suitable function, such as get.seq or other.
#' @return Returns a single character representing the residue found at the indicated position in the indicated protein.
#' @author Juan Carlos Aledo
#' @examples aa.at(28, 'P01009')
#' aa.at(at = 80, target = get.seq('P00004', 'metosite'), uniprot = FALSE)
#' @seealso is.at(), renum.pdb(), renum.meto(), renum(), aa.comp()
#' @importFrom bio3d read.fasta
#' @export

aa.at <- function(at, target, uniprot = TRUE){
  if (uniprot == TRUE){
    target <- get.seq(target, as.string = FALSE)[[1]]
  } else {
    if (length(target) == 1){
      target <- strsplit(target, split="")[[1]]
    }
  }
  if (at %in% 1:length(target)){
    return(target[at])
  } else {
    return(paste(at , " isn't a valid position for this protein", sep=""))
  }
}


## --------------------------------------------------------------- ##
#             is.at(at, target, aa = 'M', uniprot = TRUE)           #
## --------------------------------------------------------------- ##
#' Check Residue a Fixed Position
#' @description Checks if a given amino acid is at a given position.
#' @usage  is.at(at, target, aa = 'M', uniprot = TRUE)
#' @param at the position in the primary structure of the protein.
#' @param target a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param aa the amino acid of interest.
#' @param uniprot logical, if TRUE the argument 'target' should be an ID.
#' @details Please, note that when uniprot is set to FALSE, target can be the string returned by a suitable function, such as get.seq or other.
#' @return Returns a boolean. Either the residue is present at that position or not.
#' @author Juan Carlos Aledo
#' @examples is.at(28, 'P01009', 'Q')
#' is.at(at = 80, target = get.seq('P00004', 'metosite'), uniprot = FALSE)
#' @seealso aa.at(), renum.pdb(), renum.meto(), renum(), aa.comp()
#' @export

is.at <- function(at, target, aa = 'M', uniprot = TRUE){
  if (uniprot == TRUE){
    target <- get.seq(target)[[1]]
  }
  return(at %in% gregexpr(aa, target)[[1]])
}


## --------------------------------------------------------------- ##
#             aa.comp(target, uniprot = TRUE)                       #
## --------------------------------------------------------------- ##
#' Amino Acid Composition
#' @description Returns a table with the amino acid composition of the target protein.
#' @usage aa.comp(target, uniprot = TRUE)
#' @param target a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param uniprot logical, if TRUE the argument 'target' should be an ID.
#' @return Returns a dataframe with the absolute frequency of each type of residue found in the target peptide.
#' @author Juan Carlos Aledo
#' @examples aa.comp('MPSSVSWGILLLAGLCCLVPVSLAEDPQGDAAQK', uniprot = FALSE)
#' aa.comp('P01009')
#' @seealso is.at(), renum.pdb(), renum.meto(), renum(), aa.at()
#' @export

aa.comp <- function(target, uniprot = TRUE){

  if (uniprot){
    seq <- get.seq(id = target) # Protein sequence
    id <- target
  } else {
    seq <- target
    id <- 'user-provided sequence'
  }

  output <- data.frame(aa = aai$aa, frequency = NA)

  for (aa in output$aa){
    t <- gregexpr(aa, seq)[[1]]
    if (t[1] == -1){
      output$frequency[which(output$aa == aa)] <- 0
    } else {
      output$frequency[which(output$aa == aa)] <- length(t)
    }
  }
  attr(output, 'seq') <- target
  return(output)
}


## --------------------------------------------------------------- ##
#                renum.pdb(pdb, chain, uniprot)                     #
## --------------------------------------------------------------- ##
#' Renumerate Residue Position
#' @description Renumerates residue position of a PDB sequence to match the corresponding UniProt sequence.
#' @usage  renum.pdb(pdb, chain, uniprot)
#' @param pdb the PDB ID or the path to a pdb file.
#' @param chain the chain of interest.
#' @param uniprot the UniProt ID.
#' @return Returns a dataframe containing the re-numerated sequence.
#' @author Juan Carlos Aledo
#' @examples renum.pdb(pdb = '121P', chain = 'A', uniprot = 'P01112')
#' @seealso is.at(), aa.at(), renum.meto(), renum(), aa.compo()
#' @importFrom bio3d read.pdb
#' @importFrom bio3d aa321
#' @export

renum.pdb <- function(pdb, chain, uniprot){

  ## ---------------------- Protein from PDB --------------- ##
  prot_pdb <- suppressWarnings(bio3d::read.pdb(pdb)$atom)
  prot_pdb <- prot_pdb[which(prot_pdb$elety == 'CA'),]
  seq_pdb <- bio3d::aa321(prot_pdb$resid[which(prot_pdb$elety == 'CA' &
                                   prot_pdb$chain == chain)])

  ## ----------------- Protein from UniProt --------------- ##
  seq_uni <- get.seq(uniprot, as.string = FALSE)[[1]]


  ## ----------------- Proteins Alignement --------------- ##
  sequences <- c(paste(seq_uni, collapse = ""),
                 paste(seq_pdb, collapse = ""))
  ids <- c('seq_uni', 'seq_pdb')
  aln <- msa(sequences, ids)

  ## ------------------- Renumerating ------------------- ##
  valn <- as.data.frame(matrix(rep(NA, dim(aln$ali)[2]*6), ncol = 6)) # vertical alignment
  names(valn) <- c('aln_pos', 'uni_pos', 'uniprot','pdb', 'pdb_pos', 'pdb_renum')
  valn$aln_pos <- 1:dim(aln$ali)[2]
  valn$uniprot <- aln$ali[1,]
  valn$pdb <- aln$ali[2,]

  k <- 1 # UniProt Index
  j <- 1 # PDB Index
  for (i in 1:nrow(valn)){
    if (valn$uniprot[i] != '-'){
      valn$uni_pos[i] <- k
      k <- k + 1
    }
    if (valn$pdb[i] != '-'){
      valn$pdb_pos[i] <- j
      j <- j + 1
    }
    if (valn$uniprot[i] == valn$pdb[i]){
      valn$pdb_renum[i] <- valn$uni_pos[i]
    }
  }
  return(valn)
}


## --------------------------------------------------------------- ##
#                     renum.meto(uniprot)                           #
## --------------------------------------------------------------- ##
#' Renumerate Residue Position
#' @description Renumerates residue position of a MetOSite sequence to match the corresponding UniProt sequence
#' @usage  renum.meto(uniprot)
#' @param uniprot the UniProt ID.
#' @return Returns a dataframe containing the re-numerated sequence.
#' @author Juan Carlos Aledo
#' @examples renum.meto('P01009')
#' @seealso is.at(), aa.at(), renum.pdb(), renum(), aa.comp()
#' @export

renum.meto <- function(uniprot){

  ## ----------------- Proteins Alignement --------------- ##
  seq_uni <- get.seq(id = uniprot)
  seq_meto <- get.seq(id = uniprot, db = 'metosite')
  sequences <- c(seq_uni, seq_meto)
  names(sequences) <- c("uniprot", "metosite")
  aln <- msa(sequences)

  ## ------------------- Renumerating ------------------- ##
  valn <- as.data.frame(matrix(rep(NA, dim(aln$ali)[2]*6), ncol = 6)) # vertical alignment
  names(valn) <- c('aln_pos', 'uni_pos', 'uniprot','meto', 'meto_pos', 'meto_renum')
  valn$aln_pos <- 1:length(aln$ali[1,])
  valn$uniprot <- aln$ali[1,]
  valn$meto <- aln$ali[2,]

  k <- 1 # UniProt Index
  j <- 1 # PDB Index
  for (i in 1:nrow(valn)){
    if (valn$uniprot[i] != '-'){
      valn$uni_pos[i] <- k
      k <- k + 1
    }
    if (valn$meto[i] != '-'){
      valn$meto_pos[i] <- j
      j <- j + 1
    }
    if (valn$uniprot[i] == valn$meto[i]){
      valn$meto_renum[i] <- valn$uni_pos[i]
    }
  }

  return(valn)
}


## --------------------------------------------------------------- ##
#                     renum(up_id, pos, from, to, ...)              #
## --------------------------------------------------------------- ##
#' Renumerate Residue Position
#' @description Renumerates residue position.
#' @usage  renum(up_id, pos, from, to, ...)
#' @param up_id the UniProt ID.
#' @param pos position in the initial sequence.
#' @param from origin of the initial sequence, it should be one among 'uniprot', 'metosite' and 'pdb'.
#' @param to target sequence, it should be one among 'uniprot', 'metosite' and 'pdb'.
#' @param ...  additional arguments (PDB ID and chain) when 'pdb' is either origin or destination.
#' @details Either the origin sequence or the target sequence should be uniprot. Nevertheless, the conversion pdb -> metosite, for instance, can be achieved through the path: pdb -> uniprot -> metosite. If 'pdb' is selected, then the PDB ID and the involved chain must be provided, in that order.
#' @return Returns the final position.
#' @author Juan Carlos Aledo
#' @examples renum(up_id = 'P01009', pos = 351, from = 'metosite', to = 'uniprot')
#' @examples \dontrun{renum(up_id = 'P01009', pos = 60, from = 'uniprot',
#'                          to = 'pdb', pdb = '1ATU', chain = 'A')}
#' @seealso is.at(), aa.at(), renum.pdb(), renum.meto(), aa.comp()
#' @export

renum <- function(up_id, pos, from, to, ...){

  z <- list(...)
  if (length(z) != 0){
    pdb <- z[[1]][1]
    chain <- z[[2]][1]
  }

  if (from == 'uniprot' & to == 'pdb'){
     t <- renum.pdb(pdb = pdb, chain = chain, uniprot = up_id)
     rpos <- t$pdb_pos[which(t$uni_pos == pos)]
  } else if (from == 'pdb' & to == 'uniprot'){
     t <- renum.pdb(pdb = pdb, chain = chain, uniprot = up_id)
     rpos <- t$uni_pos[which(t$pdb_pos == pos)]
  } else if (from == 'uniprot' & to == 'metosite'){
     t <- renum.meto(up_id)
     rpos <- t$meto_pos[which(t$uni_pos == pos)]
  } else if(from == 'metosite' & to == 'uniprot'){
     t <- renum.meto(up_id)
     rpos <- t$uni_pos[which(t$meto_pos == pos)]
  } else {
    stop("Provide proper values for the parameters!")
  }

  return(rpos)
}
