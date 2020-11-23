# contacts.R

## ---------------------------------------------------------------- ##
#           contacts <- function(pdb,  threshold = 7)                #
## ---------------------------------------------------------------- ##
#' Search for Residue Contacts
#' @description Searchs for intra- and inter-chain residues contacts.
#' @usage contact(pdb, threshold = 7)
#' @param pdb either the path to the PDB file of interest or the 4-letters identifier.
#' @details For each residue this function cumputes the distances to all the remaining residues. When that distance is equal or lower to the threshold, it will be computed as a contact.
#' @value The function returns a dataframe where each residue from the protein represents a row. The number of intra-chain contacts is provided as well as the number of total contants (intra- and inter-chains).
#' @author Juan Carlos Aledo
#' @examples contacts('2gls')
#' contacts('./pdb/1u8f.pdb')
#' @references Fulano et al (xxxx) J. B. Chem. xxx:xx-xx.
#' @seealso res.dist() and pairwaise.distance()
#' @importFrom bio3d aa321
#' @importFrom bio3d get.pdb
#' @importFrom bio3d read.pdb
#' @export

contacts <- function(pdb, threshold = 7, rawdata = FALSE){

  ## ------ Check whether is necessary to download the pdb file ----- ##
  if (gregexpr("\\.pdb", pdb) == -1){
    pdb_id <- pdb
    bio3d::get.pdb(pdb_id, ".")
    path <- './'
    pdb <- paste(path, pdb_id, '.pdb', sep = "")
  } else {
    pdb_id <- strsplit(pdb, split = "/")[[1]]
    path <- paste(pdb_id[-length(pdb_id)], collapse = "/")
    pdb_id <- strsplit(pdb_id[length(pdb_id)], split = "\\.")[[1]][1]
  }
  ## -------------- Splitting the pdb into chains ------------------- ##
  chains_id <- pdb.chain(pdb)

  ## -------------------- Output dataframe -------------------------- ##
  aa <- c('A','R','N','D','C','Q','E','G','H','I',
          'L','K','M','F','P','S','T','W','Y','V')
  x <- bio3d::read.pdb(pdb)
  x <- x$atom # from the whole protein (which may be an oligomeric one)
  x <- x[which(x$type == "ATOM" & aa321(x$resid) %in% aa),] # remove ligands if present

  Calpha <- x[which(x$elety == 'CA' &  aa321(x$resid) %in% aa), ]
  output <- as.data.frame(matrix(rep(NA, nrow(Calpha)* 5), ncol = 5))
  names(output) <- c('residue', 'chain','intra_contacts', 'inter_contacts','total_contacts')
  output$residue <- Calpha$resid
  output$chain <- Calpha$chain

  ## --------------- Computing distances whole protein ----------------------- ##
  x$factor <- paste(x$resid, x$chain, x$resno, sep = "-")
  cumulative <- c()
  f <- unique(x$factor)

  for (i in 1:length(f)){
    cumulative <- c(cumulative, f[i])
    t <- as.matrix(x[which(x$factor == f[i]), colnames(x) %in% c('x', 'y', 'z')])
    tc <- as.matrix(x[which(! x$factor  %in% cumulative), colnames(x) %in% c('x', 'y', 'z')])
    d <- pairwise.dist(t, tc, squared = FALSE)
    output$total_contacts[i] <- length(which(d < threshold))
  }

  ## ---------------- Computing distances single chain ----------------------- ##
  chains_id <- unique(x$chain) # To avoid chains that are not protein (RNA, etc)
  counter <- 1
  for (ch in chains_id){
    x <- bio3d::read.pdb(paste(path, "/temp/", pdb_id, "_", ch, ".pdb", sep = ""))
    x <- x$atom # from the chain
    x <- x[which(x$type == "ATOM" & aa321(x$resid) %in% aa),] # remove ligands if present
    x$factor <- paste(x$resid, x$chain, x$resno, sep = "-")

    cumulative <- c()
    f <- unique(x$factor)

    for (i in 1:length(f)){
      cumulative <- c(cumulative, f[i])
      t <- as.matrix(x[which(x$factor == f[i]), colnames(x) %in% c('x', 'y', 'z')])
      tc <- as.matrix(x[which(! x$factor  %in% cumulative), colnames(x) %in% c('x', 'y', 'z')])
      d <- pairwise.dist(t, tc, squared = FALSE)
      output$intra_contacts[counter] <- length(which(d < threshold))
      counter <- counter + 1
    }
  }

  output$inter_contacts <- output$total_contacts - output$intra_contacts
  return(output)
}
