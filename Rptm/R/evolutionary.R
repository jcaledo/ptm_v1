## ----------- envolutionary.R ------------ ##
#                                            #
#     msa                                    #
#     custom.aln                             #
#     list.hom                               #
#     parse.hssp                             #
#     get.hssp                               #
#     shannon                                #
#     site.type                              #
#                                            #
## ---------------------------------------- ##

## ----------------------------------------------------------------- ##
#   msa <- function(sequences, ids = names(sequences),                #
#                     sfile = FALSE, inhouse = FALSE)                 #
## ----------------------------------------------------------------- ##
#' Multiple Sequence Alignment
#' @description Aligns multiple protein sequences.
#' @usage msa(sequences, ids = names(sequences), sfile = FALSE, inhouse = FALSE)
#' @param sequences vector containing the sequences.
#' @param ids vector containing the sequences' ids.
#' @param sfile path to the file where the fasta alignment should be saved, if any.
#' @param inhouse logical, if TRUE the in-house MUSCLE software is used. It must be installed on your system and in the search path for executables.
#' @return Returns a list of four elements. The first one ($seq) provides the sequences analyzed, the second element ($ids) returns the identifiers, the third element ($aln) provides the alignment in fasta format and the fourth element ($ali) gives the alignment in matrix format.
#' @examples \dontrun{msa(sequences = sapply(c("P19446", "P40925", "P40926"), ptm::get.seq),
#'  ids = c("wmelon", "cyt", "mit"))}
#' @references Edgar RC. Nucl. Ac. Res. 2004 32:1792-1797.
#' @references Edgar RC. BMC Bioinformatics 5(1):113.
#' @seealso custom.aln(), list.hom(), parse.hssp(), get.hssp(), shannon()
#' @importFrom Biostrings AAStringSet
#' @importFrom muscle muscle
#' @importFrom bio3d seqbind
#' @importFrom bio3d seqaln
#' @importFrom bio3d write.fasta
#' @export

msa <- function(sequences, ids = names(sequences), sfile = FALSE, inhouse = FALSE){

  ## --- Checking the inputs
  if (length(sequences) < 2){
    stop("At least two sequences are required!")
  } else if (length(sequences) != length(ids)){
    stop("The number of sequences and sequences' ids doesn't match!")
  }

  ## --- Using the program MUSCLE, which must be in the user system path
  if (inhouse){
    seqs <- lapply(sequences, function(x) strsplit(x, split = "")[[1]])
    sqs <- bio3d::seqbind(seqs[[1]], seqs[[2]], blank = "-")
    c <- 2
    while (c < length(sequences)){
      c <- c + 1
      sqs <- bio3d::seqbind(sqs, seqs[[c]], blank = "-")
    }
    aln <- bio3d::seqaln(sqs, id = ids, exefile = "muscle")
    aln$seq <- sequences

    if (sfile != FALSE){
      bio3d::write.fasta(aln, file = paste("./alignment ", Sys.time(), sep = ""))
    }
    return(aln)

  } else {
    ## --- Using the Biostring and muscle R packages
    if (requireNamespace('Biostrings', quietly = TRUE)){
      seq <- Biostrings::AAStringSet(sequences)
    } else {
      stop("You must install the package Biostrings in order to use this function")
    }

    if (requireNamespace('muscle', quietly = TRUE)){
      aln1 <- muscle::muscle(seq)
    } else {
      stop("You must install the package muscle in order to use this function")
    }

    aln <- list()
    aln$seq <- sequences
    aln$ids <- ids
    aln$aln <- as.character(aln1)
    l <- sapply(aln$aln, function (x) strsplit(x, split = ""))
    aln$ali <- matrix(unlist(l), nrow = length(sequences), byrow = TRUE)

    if (sfile != FALSE){
      for (i in 1:length(aln$aln)){
        t <- paste(">", aln$ids[i], sep = "")
        cat(t, file = sfile, append = TRUE)
        tt <- paste("\n", aln$aln[i], "\n", sep = "" )
        cat(tt, file = sfile, append = TRUE)
      }
    }
    return(aln)
  }
}


## ---------------------------------------------------------------- ##
#  custom.aln <- function(target, species, molecule, sfile = FALSE)  #
## ---------------------------------------------------------------- ##
#' Download and Align Orthologous Sequences
#' @description Downloads orthologous sequences and carries out their alignment
#' @usage custom.aln(target, species, molecule = 'protein', sfile = FALSE)
#' @param target the KEGG identifier of the protein of interest.
#' @param species a character vector containing the KEGG code for the species of interest.
#' @param molecule a character string specifying the nature of the sequence; it must be one of 'dna', 'protein'.
#' @param sfile logical, if TRUE the alignment in fasta format is saved in the current directory.
#' @details We can build the list of species or, alternatively, we can choose between four pre-established options: 'vertebrates', 'plants','one-hundred' and 'two-hundred'. The first will use the following seven species: human (hsa), chimp (ptr), gorilla (ggo), rat (rno), cow (bta), chicken (gga), western clawed frog (xtr) and zebrafish (dre). The second, A. thaliana (ara), A. lyrata (aly), B. oleracea (boe), G. max (gmax), S. lycopersicum (sly), O. sativa (osa) and C. reinhardtii (cre). The third and fourth options will use orthologous sequences from one hundred and two hundred different species, respectively.
#' @return Returns a list of class "fasta" with three components: 'ali' (an alignment character matrix with a row per sequence and a column per residue), 'id' (sequence identifiers) and 'call' (the matched call).
#' @author Juan Carlos Aledo
#' @examples \dontrun{custom.aln('hsa:4069', species = c('pps', 'pon', 'mcc', 'ssc'))}
#' \dontrun{custom.aln('cge:100773737', 'vertebrates', molecule = 'dna' )}
#' @references Edgar RC. Nucl. Ac. Res. 2004 32:1792-1797.
#' @references Edgar RC. BMC Bioinformatics 5(1):113.
#' @seealso msa(), list.hom(), parse.hssp(), get.hssp(), shannon()
#' @importFrom bio3d seqbind
#' @importFrom bio3d seqaln
#' @export

custom.aln <- function(target, species, molecule = 'protein',  sfile = FALSE){

  ## -------- Getting the KEGG orthologous identifiers ------------##
  orthDF <- list.hom(target)

  ## Cleaning the orthDF (species shouldn't be repeated):
  wrongLines <- which(orthDF$species %in% names(which(table(orthDF$species) > 1)))
  if (length(wrongLines != 0)){
    orthDF <- orthDF[-wrongLines,]
  }
  ## Cleaning the orthDF (species should have proper IDs):
  orthDF <- orthDF[which(nchar(orthDF$species) < 5),]

  ## -------- Setting the rigth subset of species:id ----------- ##
  if (species[1] == 'vertebrates'){
    id <- c('hsa', 'ptr', 'ggo', 'rno', 'bta', 'gga','xtr', 'dre')
    orthDF <- orthDF[which(orthDF$species %in% id),]
    n <- nrow(orthDF)
  } else if (species[1] == 'plants'){
    id <- c('ara', 'aly', 'boe', 'gmx', 'sly', 'osa', 'cre')
    orthDF <- orthDF[which(orthDF$species %in% id),]
    n <- nrow(orthDF)
  } else if (species[1] == 'one-hundred'){
    n <- 100
  } else if (species[1] == 'two-hundred'){
    n <- 200
  } else {
    id <- species
    orthDF <- orthDF[which(orthDF$species %in% id),]
    n <- nrow(orthDF)
  }

  ## ---------------- Downloading the sequences ------------------##
  if (molecule == 'dna'){
    db <- 'kegg-nt'
    protein <- FALSE
  } else {
    db <- 'kegg-aa'
    protein <- TRUE
  }
  ref_seq <- ptm::get.seq(target, db = db, as.string = FALSE)
  id <- c('ref')
  i <- 1
  while (length(id) <= n & i <= nrow(orthDF)){
    t <- ptm::get.seq(paste(orthDF$species[i], orthDF$entry[i], sep = ":"),
                      db = db, as.string = FALSE)
    if (!grepl("Sorry", t)){
      if (i == 1){
        sequences <- seqbind(ref_seq[[1]], t[[1]], blank = "-")
      } else {
        sequences <- seqbind(sequences, t[[1]])
      }
      id <- c(id, orthDF$species[i])
    }
    i <- i + 1
  }
  if (length(id) < n){
    status <- paste("Only ", length(id), " species could be aligned!")
  } else {
    status <- "Success"
  }

  ## ----------------- Carrying out the alignment ------------------ ##
  fa_file <- paste(gsub(":", '_', target), 'aln.fa', sep = "_")
  aln <- bio3d::seqaln(sequences, id <- id, exefile = 'muscle',
                       outfile = fa_file,
                       protein = protein)
  if (! sfile){
    if (file.exists(fa_file)){
      file.remove(fa_file)
    }
  }
  attr(aln, 'status') <- status
  return(aln)
}


## ---------------------------------------------------------------- ##
#         list.hom <- function(target, homology = 'o')               #
## ---------------------------------------------------------------- ##
#' Search Homologous Entries
#' @description Searches for homologous entries.
#' @usage list.hom(target, homology = 'o')
#' @param target the KEGG identifier of the protein of interest.
#' @param homology one letter indicating the type of homology. It should be either 'o' (orthologs) or 'p' (paralogs).
#' @details This application rests on the KEGG Sequence Similarity Database, which contains the information about amino acid sequence similarities among all protein-coding genes in the complete genomes, as well as the addendum and virus categories, of the GENES database.
#' @return Returns a dataframe with the requested entries.
#' @author Juan Carlos Aledo
#' @examples \dontrun{list.hom('hsa:2744', hom = 'p')}
#' @references Kanehisa et al (2017) Nucl. Ac. Res. 33:D353-D361.
#' @seealso msa(), custom.aln(), parse.hssp(), get.hssp(), shannon()
#' @importFrom httr GET
#' @importFrom httr http_error
#' @importFrom httr status_code
#' @importFrom httr content
#' @importFrom XML htmlTreeParse
#' @importFrom XML xpathApply
#' @importFrom XML xmlValue
#' @export

list.hom <- function(target, homology = 'o'){

  ## --------------------------- Orthologs versus Paralogs --------------------------------- ##
  if (homology == 'o'){
    url <- paste('https://www.kegg.jp/ssdb-bin/ssdb_best?org_gene=', target, sep = "")
  } else if (homology == 'p'){
    url <- paste('https://www.kegg.jp/ssdb-bin/ssdb_paralog?org_gene=', target, sep = "")
  } else {
    stop("The parameter 'hom' should be either 'o' (orthologs) or 'p' (paralogs)")
  }

  ## -------------------- Getting the KEGG homologous identifiers ---------------------------##
  resp <- httr::GET(url)
  if (httr::http_error(resp)){
    stop(paste("The request could not be properly solved. The status code was: ", status_code(resp), sep = ""))
  }

  cont <- httr::content(resp, as = 'text')
  xml_html <- XML::htmlTreeParse(cont, useInternalNodes = TRUE)
  u <- XML::xpathApply(xml_html, "//body//text()[not(ancestor::script)][not(ancestor::style)][not(ancestor::noscript)]", xmlValue)

  pos <- c() # numbers of the lines containing 'species:id'
  start_at <- which(grepl("Entry", u)) + 1
  for (i in start_at:length(u)){
    t <- u[i]
    if (grepl("[a-z]+:\\w", t)){
      pos <- c(pos, i)
    }
  }

  sp <- c()
  N <- length(pos)
  if (N > 0){
    for (i in 1:N){
      sp <-  c(sp, strsplit(u[[pos[i]]], split = ":")[[1]][1])
      artif <- which(nchar(sp) > 4) # species id longer than 4 character are artifacts
    }
  } else {
    return("No homologous found")
  }

  if (length(artif) > 0){
    pos <- pos[- artif]
    N <- length(pos)
  }

  hom <- as.data.frame(matrix(rep(NA, N*7), ncol = 7))
  names(hom) <- c("species","entry", "name", "ko", "len", "identity", "overlap")

  for (i in 1:N){
    hom$species[i] <- strsplit(u[[pos[i]]], split = ":")[[1]][1]
    hom$entry[i] <- strsplit(u[[pos[i]]], split = ":")[[1]][2]
    hom$name[i] <- gsub("^\\s+|\\s+$", "", u[[pos[i] + 1]])
    hom$ko[i] <- gsub("^\\s+|\\s+$", "", u[[pos[i] + 2]])

    if (grepl("->", hom$name[i])){ # lacks KO and the name is an unsplitted string. Need to be amended
      hom$len[i] <- NA          # They will be completed later
      hom$identity[i] <- NA
      hom$overlap[i] <- NA
    } else {
      tr <- strsplit(u[[pos[i] + 3]], split = "\\s")[[1]]
      tr <- tr[which(tr != "")]
      hom$len[i] <- as.numeric(tr[1])
      hom$identity[i] <- as.numeric(tr[6])
      hom$overlap[i] <- as.numeric(tr[7])
    }
  }

  amend <- which(grepl("->", hom$name)) # lacks KO and the name is an unsplitted string. Need to be amended
  for (i in seq_len(length(amend))){
    t <- hom$name[amend[i]]
    t <- strsplit(t, split = "\\s{2,}")[[1]]
    hom$name[amend[i]] <- sub("^\\s+|\\s+$", "",t[1])
    hom$ko[amend[i]] <- NA
    hom$len[amend[i]] <- as.numeric(t[2])
    hom$identity[amend[i]] <- as.numeric(t[6])
    hom$overlap[amend[i]] <- as.numeric(t[7])
  }

  attr(hom, 'target') <- target
  attr(hom, 'homlogy') <- homology
  return(hom)
}


## ------------------------------------------------------------------------------- ##
#             parse.hssp <- function(file, keepfiles = FALSE)                       #
## ------------------------------------------------------------------------------- ##
#' Parse a HSSP File to Return Dataframes
#' @description Parses a HSSP file to return dataframes.
#' @usage parse.hssp(file, keepfiles = FALSE)
#' @param file input hssp file.
#' @param keepfiles logical, if TRUE the dataframes will be saved in the working directory and we will keep the hssp file.
#' @details If the argument 'keepfiles' is not set to TRUE, the hssp file used to get the parsed dataframe will be removed. Otherwise, 4 dataframes will be saved:
#' \itemize{
#' \item{id_seq_list.Rda:}  {This block of information holds the metadata per sequence, and some alignment statistic. See https://swift.cmbi.umcn.nl/gv/hssp for a detailed description of the information that can be find in this block.}
#' \item{id_aln.Rda:}  {This dataframe contains the alignment itself (each sequence is a column). Additional information such as secondary structure, SASA, etc., is also found in this block.}
#' \item{id_profile.Rda:}  {This dataframe holds per amino acid type its percentage in the list of residues observed at that position. In addition, this dataframe also informs about the entropy at each position, as well as the number of sequences spanning this position (NOOC).}
#' \item{id_insertions.Rda:} {A dataframe with information regarding those sequences that contain insertions. See https://swift.cmbi.umcn.nl/gv/hssp for further details.}
#' }
#' @return Returns a dataframe corresponding to the profile.Rda described above.
#' @author Juan Carlos Aledo
#' @examples \dontrun{parse.hssp(file = './1u8f.hssp')}
#' @references Touw et al (2015) Nucl. Ac. Res. 43:D364-368.
#' @seealso msa(), custom.aln(), list.hom(), get.hssp(), shannon()
#' @export

parse.hssp <- function(file, keepfiles = FALSE){

  ## ------- Check file argument --------- ##
  if (!file.exists(file)){
    stop("Please provide a proper file")
  } else if (grepl("hssp", file)){
    t <- strsplit(file, split = "/")[[1]]
    id <- strsplit(t[length(t)], split = "\\.")[[1]][1]
  } else {
    stop("Please, provide a file of type hssp")
  }

  con <- file(file, 'r')
  l <- readLines(con)
  close(con)

  ## -------------------------------- The Sequence List --------------------------------------- ##
  # ID: EMBL/SWISSPROT identifier of the aligned (homologous) protein.
  # IDE: Percentage of residue identity of the alignment.
  # WSIM: Weighted similarity of the alignment.
  # IFIR: First residue of the alignment in the test sequence.
  # ILAS: Last residue of the alignment in the test sequence.
  # JFIR: First residue of the alignment in the alignend protein.
  # JLAS: Last residue of the alignment in the alignend protein.
  # LALI: Length of the alignment excluding insertions and deletions.
  # NGAP: Number of insertions and deletions in the alignment.
  # LGAP: Total length of all insertions and deletions.
  # LSEQ2: Length of the entire sequence of the aligned protein.
  # ACCNUM: SwissProt accession number.
  # PROTEIN: One-line description of aligned protein.

  start_at <- which(grepl("  NR.    ID         STRID   %IDE %WSIM", l)) + 1 # first line for headers
  end_at <- which(grepl("## ALIGNMENTS ",  l))[1] - 1 # last line for next block's headers
  number_lines <- end_at - start_at + 1
  seq_list <- data.frame(NR = 1:number_lines,
                         ID = rep(NA, number_lines),
                         IDE = rep(NA, number_lines),
                         WSIM = rep(NA, number_lines),
                         IFIR = rep(NA, number_lines),
                         ILAS = rep(NA, number_lines),
                         JFIR = rep(NA, number_lines),
                         JLAS = rep(NA, number_lines),
                         LALI = rep(NA, number_lines),
                         NGAP = rep(NA, number_lines),
                         LGAP = rep(NA, number_lines),
                         LESEQ2 = rep(NA, number_lines),
                         ACCNUM = rep(NA, number_lines),
                         PROTEIN =  rep(NA, number_lines)
  )

  for (i in 1:number_lines){
    t <- strsplit(l[i + start_at -1], split = "\\s+")[[1]]
    seq_list$ID[i] <- t[4]
    seq_list$IDE[i] <- as.numeric(t[5])
    seq_list$WSIM[i] <- as.numeric(t[6])
    seq_list$IFIR[i] <- as.numeric(t[7])
    seq_list$ILAS[i] <- as.numeric(t[8])
    seq_list$JFIR[i] <- as.numeric(t[9])
    seq_list$JLAS[i] <- as.numeric(t[10])
    seq_list$LALI[i] <- as.numeric(t[11])
    seq_list$NGAP[i] <- as.numeric(t[12])
    seq_list$LGAP[i] <- as.numeric(t[13])
    seq_list$LESEQ2[i] <- as.numeric(t[14])
    seq_list$ACCNUM[i] <- t[15]
    seq_list$PROTEIN[i] <- paste(t[16:length(t)], collapse = " ")
  }

  if (keepfiles){
    save(seq_list, file = paste(id, "_seq_list.Rda", sep = ""))
  }
  ## ----------------------------------- The Actual Alignment --------------------------------------- ##
  # SeqNo: Sequence residue number.
  # PDBNo: PDB residue number.
  # Chain: Chain ID.
  # AA: Amino acid.
  # SS: Secondary structure.
  # ACC: Solvent exposure as in DSSP.
  # NOCC: Number of aligned sequences spanning this position (including the test sequence).
  # VAR: Sequence variability on a scale of 0-100 as derived from the NALIGN (number of sequences aligned) alignments.

  head <- which(grepl("## ALIGNMENTS", l))
  nres <- head[2] - head[1] -2 # two header lines
  aln <- data.frame(SeqNo = rep(NA, nres),
                    PDBNo =rep(NA, nres),
                    Chain = rep(NA, nres),
                    AA = rep(NA, nres),
                    SS = rep(NA, nres),
                    ACC = rep(NA, nres),
                    NOCC = rep(NA, nres),
                    VAR = rep(NA, nres))


  # Structural data:
  c <- 0
  for (j in (head[1]+2):(head[2]-1)){ # as many lines as residues
    c <- c + 1
    t <- l[j] # each line is 121 char long
    aln$SeqNo[c] <- as.numeric(gsub("^\\s+|\\s+$", "", substr(t, start =  1, stop = 6)))
    aln$PDBNo[c] <- as.numeric(gsub("^\\s+|\\s+$", "", substr(t, start =  7, stop = 11)))
    aln$Chain[c] <- gsub("^\\s+|\\s+$", "", substr(t, start =  12, stop = 13))
    aln$AA[c] <- gsub("^\\s+|\\s+$", "", substr(t, start =  14, stop = 16))
    aln$SS[c] <- gsub("^\\s+|\\s+$", "", substr(t, start =  17, stop = 19))
    aln$ACC[c] <- as.numeric(gsub("^\\s+|\\s+$", "", substr(t, start =  37, stop = 40)))
    aln$NOCC[c] <- as.numeric(gsub("^\\s+|\\s+$", "", substr(t, start =  42, stop = 46)))
    aln$VAR[c] <- as.numeric(gsub("^\\s+|\\s+$", "", substr(t, start = 47, stop = 50)))
  }

  # Alignment:
  for (i in 1:length(head)){ # Iterate through blocks (of 70 sequences in the alignment)

    c <- 0
    block <- as.data.frame(matrix(rep(NA, nres * 70), ncol = 70))

    for (j in (head[i]+2):(head[i]+1+nres)){ # Within each block (of 70 sequences in the alignment)
      c <- c + 1
      t <- l[j] # each line is 121 char long
      a <- gsub(" ", "-", strsplit(substring(t, 52,121), split = "")[[1]])
      if (identical(a, character(0))){ # empty line (see block 1, line 50 of 5xge)
        block[c,] <- rep(NA, 70)
      } else {
        block[c,] <-  gsub(" ", "-", strsplit(substring(t, 52,121), split = "")[[1]])
      }
    }

    aln <- cbind(aln, block)
  }

  columns_to_be_remove <- which(apply(aln, 2, function(x) all(x == "-"))) # columns with all rows "-"
  if (length(columns_to_be_remove) > 0){
    aln <- aln[, -columns_to_be_remove] # remove
  }
  names(aln) <- c(names(aln)[1:8], seq_list$ACCNUM) # name each seq in the alignmen by its accnum

  if (keepfiles){
    save(aln, file = paste(id, "_aln.Rda", sep = ""))
  }

  ## -------------------------------------- The Profile ------------------------------------------ ##
  # Relative Frequency of an Amino Acid Type at Each Position (scaled to 100).
  # Asx and Glx are in their acid/amide form in proportion to their database frequencies.
  # SeqNo: Sequence residue number.
  # PDBNo: PDB residue number.
  # NOCC: Number of aligned sequences spanning this position (including the test sequence).
  # NDEL: Number of sequences with a deletion in the test protein at this position.
  # NINS: Number of sequences with an insertion in the test protein at this position.
  # ENTROPY: Entropy measure of sequence variability at this position.
  # RELENT: Relative entropy, i.e.  entropy normalized to the range 0-100.
  # WEIGHT: Conservation weight.


  at <- which(grepl("## SEQUENCE PROFILE AND ENTROPY", l)) + 1
  n <- strsplit(l[at], split = "\\s+")[[1]][-1] # remove the first and the two last
  n <- n[-length(n)]
  n <- n[-length(n)]


  profile <- as.data.frame(matrix(rep(NA, nres * length(n)), ncol = length(n)))
  names(profile) <- n

  c <- 0
  for (i in (at+1):(at+nres)){
    c <- c + 1
    t <- strsplit(l[i], split = "\\s+")[[1]][-1]
    profile$SeqNo[c] <- as.numeric(t[1])
    profile$PDBNo[c] <- as.numeric(t[2])
    profile$V[c] <- as.numeric(t[4])
    profile$L[c] <- as.numeric(t[5])
    profile$I[c] <- as.numeric(t[6])
    profile$M[c] <- as.numeric(t[7])
    profile$F[c] <- as.numeric(t[8])
    profile$W[c] <- as.numeric(t[9])
    profile$Y[c] <- as.numeric(t[10])
    profile$G[c] <- as.numeric(t[11])
    profile$A[c] <- as.numeric(t[12])
    profile$P[c] <- as.numeric(t[13])
    profile$S[c] <- as.numeric(t[14])
    profile$T[c] <- as.numeric(t[15])
    profile$C[c] <- as.numeric(t[16])
    profile$H[c] <- as.numeric(t[17])
    profile$R[c] <- as.numeric(t[18])
    profile$K[c] <- as.numeric(t[19])
    profile$Q[c] <- as.numeric(t[20])
    profile$E[c] <- as.numeric(t[21])
    profile$N[c] <- as.numeric(t[22])
    profile$D[c] <- as.numeric(t[23])
    profile$NOCC[c] <- as.numeric(t[24])
    profile$NDEL[c] <- as.numeric(t[25])
    profile$NINS[c] <- as.numeric(t[26])
    profile$ENTROPY[c] <- as.numeric(t[27])
    profile$RELENT[c] <- as.numeric(t[28])
    profile$WEIGHT[c] <- as.numeric(t[29])
  }

  profile <- profile[, -c(29:32)] # remove the last four columns
  profile <- profile[, c(1:9, 11,10,12:15, 19,21,16:18,20,22:28)] # ordered by physicochemical properties

  if (keepfiles){
    save(profile, file = paste(id, "_profile.Rda", sep = ""))
  }
  ## -------------------------------------- The Insection List------------------------------------------ ##

  at <- which(grepl("## INSERTION LIST", l)) + 1
  end <-  which(grepl("^//$", l))
  nlines <- (end-1)-(at+1)

  inser <- as.data.frame(matrix(rep(NA, nlines * 5), ncol = 5))
  names(inser) <- c("AliNo", "IPOS", "JPOS", "Len", "Sequence")

  w <- FALSE
  if (nlines < 1){
    inser <- "No insertions recorded"
    warning("No insertions recorded")
    w <- TRUE
  } else {
    c <- 0
    for (i in (at+1):(end-1)){
      c <- c + 1
      t <- strsplit(l[i], split = "\\s+")[[1]][-1]
      inser[c,] <- t
    }
  }

  if (keepfiles & !w){
    save(inser, file = paste(id, "_insertions.Rda", sep = ""))
  }

  return(profile)
}


## ------------------------------------------------------------------------------- ##
#               get.hssp <- function(pdb, path)                                     #
## ------------------------------------------------------------------------------- ##
#' Get a HSSP File
#' @description Gets a HSSP file of the requested structure.
#' @usage get.hssp(pdb, path, keepfiles = TRUE)
#' @param pdb the 4-letter identifier of the PDB file.
#' @param path character string providing the path to the in-house HSSP database.
#' @param keepfiles logical, if TRUE the dataframes will be saved in the working directory and we will keep the hssp file.
#' @details In order to use this function, you need to obtain a local copy of the HSSB database. This function will obtain and parse the requested HSSP file. When the argument ‘keepfiles’ is set to TRUE, the get.hssp() function will build and save (in the working directory) the following 4 dataframes:
#' \itemize{
#' \item{id_seq_list.Rda:}  {This block of information holds the metadata per sequence, and some alignment statistic. See https://swift.cmbi.umcn.nl/gv/hssp for a detailed description of the information that can be find in this block.}
#' \item{id_aln.Rda}  {This dataframe contains the alignment itself (each sequence is a column). Additional information such as secondary structure, SASA, etc., is also found in this block.}
#' \item{id_profile.Rda}  {This dataframe holds per amino acid type its percentage in the list of residues observed at that position. In addition, this dataframe also informs about the entropy at each position, as well as the number of sequences spanning this position (NOOC).}
#' \item{id_insertions.Rda} {A dataframe with information regarding those sequences that contain insertions. See https://swift.cmbi.umcn.nl/gv/hssp for further details.}
#' }
#' @return Returns a dataframe corresponding to the profile.Rda described above.
#' @author Juan Carlos Aledo
#' @examples \dontrun{get.hssp(file = './1u8f.hssp')}
#' @references Touw et al (2015) Nucl. Ac. Res. 43:D364-368.
#' @seealso msa(), custom.aln(), list.hom(), parse.hssp(), shannon()
#' @export

get.hssp <- function(pdb, path = "/Users/juancarlosaledo/Dropbox/local_HSSP/", keepfiles = TRUE){

  if (nchar(pdb) != 4){
    stop("A proper PDB ID should be provided")
  } else {
    file <- paste(path, pdb, ".hssp.bz2", sep = "")
    id <- pdb
  }

  if (file.exists(file)){
    con <- file(file, "r")
    l <- readLines(con)
    close(con)
  } else {
    profile <- "Sorry, the requested file couldn't be found!"
    return(profile)
  }

  ## -------------------------------- The Sequence List --------------------------------------- ##
  # ID: EMBL/SWISSPROT identifier of the aligned (homologous) protein.
  # IDE: Percentage of residue identity of the alignment.
  # WSIM: Weighted similarity of the alignment.
  # IFIR: First residue of the alignment in the test sequence.
  # ILAS: Last residue of the alignment in the test sequence.
  # JFIR: First residue of the alignment in the alignend protein.
  # JLAS: Last residue of the alignment in the alignend protein.
  # LALI: Length of the alignment excluding insertions and deletions.
  # NGAP: Number of insertions and deletions in the alignment.
  # LGAP: Total length of all insertions and deletions.
  # LSEQ2: Length of the entire sequence of the aligned protein.
  # ACCNUM: SwissProt accession number.
  # PROTEIN: One-line description of aligned protein.

  start_at <- which(grepl("  NR.    ID         STRID   %IDE %WSIM", l)) + 1 # first line for headers
  end_at <- which(grepl("## ALIGNMENTS ",  l))[1] - 1 # last line for next block's headers
  number_lines <- end_at - start_at + 1
  seq_list <- data.frame(NR = 1:number_lines,
                         ID = rep(NA, number_lines),
                         IDE = rep(NA, number_lines),
                         WSIM = rep(NA, number_lines),
                         IFIR = rep(NA, number_lines),
                         ILAS = rep(NA, number_lines),
                         JFIR = rep(NA, number_lines),
                         JLAS = rep(NA, number_lines),
                         LALI = rep(NA, number_lines),
                         NGAP = rep(NA, number_lines),
                         LGAP = rep(NA, number_lines),
                         LESEQ2 = rep(NA, number_lines),
                         ACCNUM = rep(NA, number_lines),
                         PROTEIN =  rep(NA, number_lines)
  )

  for (i in 1:number_lines){
    t <- strsplit(l[i + start_at -1], split = "\\s+")[[1]]
    seq_list$ID[i] <- t[4]
    seq_list$IDE[i] <- as.numeric(t[5])
    seq_list$WSIM[i] <- as.numeric(t[6])
    seq_list$IFIR[i] <- as.numeric(t[7])
    seq_list$ILAS[i] <- as.numeric(t[8])
    seq_list$JFIR[i] <- as.numeric(t[9])
    seq_list$JLAS[i] <- as.numeric(t[10])
    seq_list$LALI[i] <- as.numeric(t[11])
    seq_list$NGAP[i] <- as.numeric(t[12])
    seq_list$LGAP[i] <- as.numeric(t[13])
    seq_list$LESEQ2[i] <- as.numeric(t[14])
    seq_list$ACCNUM[i] <- t[15]
    seq_list$PROTEIN[i] <- paste(t[16:length(t)], collapse = " ")
  }

  if (keepfiles){
    save(seq_list, file = paste(id, "_seq_list.Rda", sep = ""))
  }
  ## ----------------------------------- The Actual Alignment --------------------------------------- ##
  # SeqNo: Sequence residue number.
  # PDBNo: PDB residue number.
  # Chain: Chain ID.
  # AA: Amino acid.
  # SS: Secondary structure.
  # ACC: Solvent exposure as in DSSP.
  # NOCC: Number of aligned sequences spanning this position (including the test sequence).
  # VAR: Sequence variability on a scale of 0-100 as derived from the NALIGN (number of sequences aligned) alignments.

  head <- which(grepl("## ALIGNMENTS", l))
  nres <- head[2] - head[1] -2 # two header lines
  aln <- data.frame(SeqNo = rep(NA, nres),
                    PDBNo =rep(NA, nres),
                    Chain = rep(NA, nres),
                    AA = rep(NA, nres),
                    SS = rep(NA, nres),
                    ACC = rep(NA, nres),
                    NOCC = rep(NA, nres),
                    VAR = rep(NA, nres))


  # Structural data:
  c <- 0
  for (j in (head[1]+2):(head[2]-1)){ # as many lines as residues
    c <- c + 1
    t <- l[j] # each line is 121 char long
    aln$SeqNo[c] <- as.numeric(gsub("^\\s+|\\s+$", "", substr(t, start =  1, stop = 6)))
    aln$PDBNo[c] <- as.numeric(gsub("^\\s+|\\s+$", "", substr(t, start =  7, stop = 11)))
    aln$Chain[c] <- gsub("^\\s+|\\s+$", "", substr(t, start =  12, stop = 13))
    aln$AA[c] <- gsub("^\\s+|\\s+$", "", substr(t, start =  14, stop = 16))
    aln$SS[c] <- gsub("^\\s+|\\s+$", "", substr(t, start =  17, stop = 19))
    aln$ACC[c] <- as.numeric(gsub("^\\s+|\\s+$", "", substr(t, start =  37, stop = 40)))
    aln$NOCC[c] <- as.numeric(gsub("^\\s+|\\s+$", "", substr(t, start =  42, stop = 46)))
    aln$VAR[c] <- as.numeric(gsub("^\\s+|\\s+$", "", substr(t, start = 47, stop = 50)))
  }

  # Alignment:
  for (i in 1:length(head)){ # Iterate through blocks (of 70 sequences in the alignment)

    c <- 0
    block <- as.data.frame(matrix(rep(NA, nres * 70), ncol = 70))

    for (j in (head[i]+2):(head[i]+1+nres)){ # Within each block (of 70 sequences in the alingment)
      c <- c + 1
      t <- l[j] # each line is 121 char long
      a <- gsub(" ", "-", strsplit(substring(t, 52,121), split = "")[[1]])
      if (identical(a, character(0))){ # empty line (see block 1, line 50 of 5xge)
        block[c,] <- rep(NA, 70)
      } else {
        block[c,] <-  gsub(" ", "-", strsplit(substring(t, 52,121), split = "")[[1]])
      }
    }

    aln <- cbind(aln, block)
  }

  columns_to_be_remove <- which(apply(aln, 2, function(x) all(x == "-"))) # columns with all rows "-"
  if (length(columns_to_be_remove) > 0){
    aln <- aln[, -columns_to_be_remove] # remove
  }
  names(aln) <- c(names(aln)[1:8], seq_list$ACCNUM) # name each seq in the alignmen by its accnum

  if (keepfiles){
    save(aln, file = paste(id, "_aln.Rda", sep = ""))
  }

  ## -------------------------------------- The Profile ------------------------------------------ ##
  # Relative Frequency of an Amino Acid Type at Each Position (scaled to 100).
  # Asx and Glx are in their acid/amide form in proportion to their database frequencies.
  # SeqNo: Sequence residue number.
  # PDBNo: PDB residue number.
  # NOCC: Number of aligned sequences spanning this position (including the test sequence).
  # NDEL: Number of sequences with a deletion in the test protein at this position.
  # NINS: Number of sequences with an insertion in the test protein at this position.
  # ENTROPY: Entropy measure of sequence variability at this position.
  # RELENT: Relative entropy, i.e.  entropy normalized to the range 0-100.
  # WEIGHT: Conservation weight.


  at <- which(grepl("## SEQUENCE PROFILE AND ENTROPY", l)) + 1
  n <- strsplit(l[at], split = "\\s+")[[1]][-1] # remove the first and the two last
  n <- n[-length(n)]
  n <- n[-length(n)]


  profile <- as.data.frame(matrix(rep(NA, nres * length(n)), ncol = length(n)))
  names(profile) <- n

  c <- 0
  for (i in (at+1):(at+nres)){
    c <- c + 1
    t <- strsplit(l[i], split = "\\s+")[[1]][-1]
    profile$SeqNo[c] <- as.numeric(t[1])
    profile$PDBNo[c] <- as.numeric(t[2])
    profile$V[c] <- as.numeric(t[4])
    profile$L[c] <- as.numeric(t[5])
    profile$I[c] <- as.numeric(t[6])
    profile$M[c] <- as.numeric(t[7])
    profile$F[c] <- as.numeric(t[8])
    profile$W[c] <- as.numeric(t[9])
    profile$Y[c] <- as.numeric(t[10])
    profile$G[c] <- as.numeric(t[11])
    profile$A[c] <- as.numeric(t[12])
    profile$P[c] <- as.numeric(t[13])
    profile$S[c] <- as.numeric(t[14])
    profile$T[c] <- as.numeric(t[15])
    profile$C[c] <- as.numeric(t[16])
    profile$H[c] <- as.numeric(t[17])
    profile$R[c] <- as.numeric(t[18])
    profile$K[c] <- as.numeric(t[19])
    profile$Q[c] <- as.numeric(t[20])
    profile$E[c] <- as.numeric(t[21])
    profile$N[c] <- as.numeric(t[22])
    profile$D[c] <- as.numeric(t[23])
    profile$NOCC[c] <- as.numeric(t[24])
    profile$NDEL[c] <- as.numeric(t[25])
    profile$NINS[c] <- as.numeric(t[26])
    profile$ENTROPY[c] <- as.numeric(t[27])
    profile$RELENT[c] <- as.numeric(t[28])
    profile$WEIGHT[c] <- as.numeric(t[29])
  }

  profile <- profile[, -c(29:32)] # remove the last four columns
  profile <- profile[, c(1:9, 11,10,12:15, 19,21,16:18,20,22:28)] # ordered by physicochemical properties

  if (keepfiles){
    save(profile, file = paste(id, "_profile.Rda", sep = ""))
  }
  ## -------------------------------------- The Insection List------------------------------------------ ##

  at <- which(grepl("## INSERTION LIST", l)) + 1
  end <-  which(grepl("^//$", l))
  nlines <- (end-1)-(at+1)

  inser <- as.data.frame(matrix(rep(NA, nlines * 5), ncol = 5))
  names(inser) <- c("AliNo", "IPOS", "JPOS", "Len", "Sequence")

  w <- FALSE
  if (nlines < 1){
    inser <- "No insertions recorded"
    warning("No insertions recorded")
    w <- TRUE
  } else {
    c <- 0
    for (i in (at+1):(end-1)){
      c <- c + 1
      t <- strsplit(l[i], split = "\\s+")[[1]][-1]
      inser[c,] <- t
    }
  }

  if (keepfiles & !w){
    save(inser, file = paste(id, "_insertions.Rda", sep = ""))
  }

  return(profile)
}


## ---------------------------------------------------------------- ##
#   shannon <- function(target, species, base = 2, alphabet = 21)    #
## ---------------------------------------------------------------- ##
#' Compute Shannon Entropy
#' @description Computes Shannon's entropies for both amino acid and codon positions
#' @usage shannon(target, species, base = 2, alphabet = 21)
#' @param target the KEGG identifier of the protein of interest.
#' @param species a character vector containing the KEGG code for the species of interest.
#' @param base integer that must take an allowed value: 2, 4 or 21 (see details).
#' @param alphabet a numeric value that can be either 21 or 4 (see details).
#' @details To compute the entropy at a given position in an amino acid alignment, we can consider an alphabet of 21 (20 amino acids plus the gap symbol, "-"), or alternatively an alphabet of 4 symbols when the amino acids are grouped according to their properties: charged (E,D,H,R,K), hydrophobic (A,L,I,V,M,F,W,Y), polar (S,T,C,Q,N) and special (G,P,-). The logarithm base used to compute the Shannon entropy can be set to 2 when we want to interpret the entropy in terms of bits, or it can be chosen to be  4 or 21 in order to get relative entropy values (ranging from 0 to 1) when we are using 4 or 21 letters alphabets, respectively. In the case of the codon entropy, we always use an alphabet of 62 codons. The logarithm base can be set to 2, in any other case it will take the value 62. Regarding the species included in the alingment, we can build the list of species or, alternatively, we can choose between pre-established options: 'vertebrates', 'plants', 'one-hundred', 'two-hundred'. The first will use the following seven species: human (hsa), chimp (ptr), gorilla (ggo), rat (rno), cow (bta), chicken (gga), western clawed frog (xtr) and zebrafish (dre). The second, A. thaliana (ara), A. lyrata (aly), B. oleracea (boe), G. max (gmax), S. lycopersicum (sly), O. sativa (osa) and C. reinhardtii (cre). The third and fourth options will use orthologous sequences from one hundred and two hundred different species, respectively.
#' @return Returns the computed entropy for each position of the molecular sequence of reference (both protein and cDNA sequences are considered).
#' @author Juan Carlos Aledo
#' @examples \dontrun{shannon('hsa:4069', 'vertebrates')}
#' @seealso msa(), custom.aln(), parse.hssp(), get.hssp()
#' @importFrom seqinr translate
#' @importFrom seqinr s2c
#' @export

shannon <- function(target, species, base = 2, alphabet = 21){

  ## ----------------------------- Checking arguments ---------------------------- ##
  if (! base %in% c(2, 4, 21)){
    base <- 2
    warn <- "The base has been set to 2"
  }
  if( ! alphabet %in% c(4, 21)){
    alphabet <- 21
    warn <- "The alphabet has been set to 21"
  }

  ## ---------------------------- Getting the alignment -------------------------- ##
  aln_aa <- custom.aln(target, species)
  aln_nt <- custom.aln(target, species, molecule = 'dna')


  ## ----------------------------- Reference Sequence ---------------------------- ##
  # The reference protein sequence:
  ref_aa <- aln_aa$ali[which(aln_aa$id == "ref"),]
  names(ref_aa) <- 1:length(ref_aa) # numerated according to the alignment
  g <- ref_aa[which(ref_aa == "-")] # gaps in ref seq
  ref_aa <- ref_aa[which(ref_aa != "-")] # removing gaps

  # The reference dna sequence:
  ref_nt <- aln_nt$ali[which(aln_nt$id == "ref"),]
  # split into codons:
  ref_nt <- paste(ref_nt[which(ref_nt != "-")], collapse = "")
  ref_nt <- strsplit(gsub('(.{3})', '\\1 ', ref_nt), split = " ")[[1]]
  if (seqinr::translate(seqinr::s2c(ref_nt[length(ref_nt)])) == "*"){
    ref_nt <- ref_nt[-length(ref_nt)] # remove stop codon
  }

  ## ------------------------ Working Sequence-Dataframes ------------------------ ##
  nr <- length(aln_aa$id)
  nc <- length(ref_aa)
  aa <- matrix(rep(NA, nr*nc), ncol = nc) # --- For amino acids
  rownames(aa) <- aln_aa$id
  colnames(aa) <- names(ref_aa)
  nt <- aa # ------- For nucleotides

  aa[1,] <- ref_aa
  nt[1,] <- ref_nt

  ## ------------------------  Output Dataframes -------------------------------- ##
  output <- as.data.frame(matrix(rep(NA, 5 * ncol(aa)), ncol = 5))
  names(output) <- c('n', 'aa', 'codon', 'Haa', 'Hcodon')
  output$n <- 1:nrow(output)
  output$aa <- ref_aa
  output$codon <- ref_nt

  ## ---------------- Cicling through orthologous sequences --------------------- ##
  mismatchings <- c()
  pos <- as.numeric(names(ref_aa))
  for (i in 2:nr){
    prot <- aln_aa$ali[i,]
    ngap <- which(prot != "-")
    dna <- aln_nt$ali[i,]
    dna <- paste(dna[which(dna != '-')], collapse = "")
    dna <- strsplit(gsub('(.{3})', '\\1 ', dna), split = " ")[[1]]

    if (nchar(dna[length(dna)]) == 3){
      if (seqinr::translate(seqinr::s2c(dna[length(dna)])) == "*"){
        dna <- dna[-length(dna)] # remove stop codon when needed
      }
    }

    count <- 1
    nombres <- c()
    for (j in 1:length(prot)){
      if (prot[j] == "-"){
        nombres <- c(nombres, 'g')
      } else {
        nombres <- c(nombres, dna[count])
        count <- count + 1
      }
    }
    prot <- prot[pos]
    nombres <- nombres[pos]

    ## --- Quality Control --- ##
    secuenciaProteina <- paste(prot[which(prot != '-')], collapse = "")
    cDNA <- paste(nombres[which(nombres != "g")], collapse = "")
    cDNA_translated <- paste(seqinr::translate(seqinr::s2c(cDNA)), collapse = "")
    if (secuenciaProteina == cDNA_translated){
      aa[i,] <- prot
      nt[i,] <- nombres
    } else {
      mismatchings <- c(mismatchings, aln_aa$id[i])
    }
  }

  ## ----------- Computing Shannon's Entropies for Codons ------------------------ ##
  tablas <- apply(nt, 2, table)
  tablas <- lapply(tablas, as.data.frame)

  # For each codon in the sequence of reference, we have $`pos`. Example:
  # $ `87`
  #    Var1 Freq
  # 1  AGC    1
  # 2  GCC    1
  # 3  GGC    2

  if (base != 2){
    base_codon <- 62
  } else {
    base_codon <- 2
  }
  H <- function(x) { # For each position we compute the Shannon's entropy
    -sum( ( x[,2]/(sum( x[,2] ) ) )*log(x[,2]/(sum(x[,2])), base= base_codon) )
  }

  entropy <- sapply(tablas, H) # Vector with the Shannon entropies for each position
  entropy <- sapply(entropy, function(x) round(x, digits=3))
  output$Hcodon <- entropy

  ## ----------- Computing Shannon's Entropies for Amino Acids ------------------ ##
  if (alphabet == 4){
    charged <- c('E', 'D', 'H', 'K', 'R')
    hydrophobic <- c('L', 'I', 'V', 'M', 'F', 'W', 'Y', 'A')
    polar <- c('S', 'T', 'C', 'Q', 'N')
    special <- c('G','P','-')

    for (i in seq_len(nrow(output))){
      S <- 0
      X <- c()
      X <- c(X, sum(aa[,i] %in% charged)/sum(!is.na(aa[,i]))) # NA doesn't count
      X <- c(X, sum(aa[,i] %in% hydrophobic)/sum(!is.na(aa[,i])))
      X <- c(X, sum(aa[,i] %in% polar)/sum(!is.na(aa[,i])))
      X <- c(X, sum(aa[,i] %in% special)/sum(!is.na(aa[,i])))
      X <- X[which(X != 0)]
      S <- sum(X*log(1/X, base = base))
      output$Haa[i] <- round(S, 3)
    }

  } else{
    tablas <- apply(aa, 2, table)
    tablas <- lapply(tablas, as.data.frame)

    # For each amino acid in the sequence of reference, we have $`pos`. Example:
    # $`87`
    #     Var1 Freq
    # 1    A    1
    # 2    G    2
    # 3    S    1

    H <- function(x) { # For each position we compute the Shannon's entropy
      -sum( ( x[,2]/(sum( x[,2] ) ) )*log(x[,2]/(sum(x[,2])), base=base) )
    }

    entropy <- sapply(tablas, H) # Vector with the Shannon entropies for each position
    entropy <- sapply(entropy, function(x) round(x, digits=3))
    output$Haa <- entropy
  }

  ## ----------------  Output ---------------------- ##
  if (is.null(mismatchings)){
    attr(output, 'mismatchings') <- 'Any'
  } else {
    attr(output, 'mismatchings') <- mismatchings
  }
  attr(output, 'base') <- base
  attr(output, 'alphabet') <- alphabet
  return(output)
}


## ---------------------------------------------------------------- ##
#          site.type <- function(target, species, th = 0.25)         #
## ---------------------------------------------------------------- ##
#' Compute Shannon Entropy and Sort out the Sites
#' @description Computes Shannon's entropies and performs a partition of the sites set.
#' @usage site.type(target, species, th = 0.25)
#' @param target the KEGG identifier of the protein of interest.
#' @param species a character vector containing the KEGG code for the species of interest.
#' @param th value between 0 and 1 indicating the percentile driving the site partition.
#' @details Each site can be classified according to their entropies into the following categories: invariant, pseudo-invariant, constrained, conservative, unconstrained, drastic.
#' @return Returns a dataframe including the category of each site according to its variability.
#' @author Juan Carlos Aledo
#' @examples \dontrun{site.type('hsa:4069', 'vertebrates')}
#' @seealso msa(), custom.aln(), parse.hssp(), get.hssp(), shannon()
#' @importFrom seqinr translate
#' @importFrom seqinr s2c
#' @importFrom stats quantile
#' @export

site.type <- function(target, species, th = 0.25){

  ## ---------------------------- Getting the alignment -------------------------- ##
  aln_aa <- custom.aln(target, species)

  ## ----------------------------- Reference Sequence ---------------------------- ##
  # The reference protein sequence:
  ref_aa <- aln_aa$ali[which(aln_aa$id == "ref"),]
  names(ref_aa) <- 1:length(ref_aa) # numerated according to the alignment
  g <- ref_aa[which(ref_aa == "-")] # gaps in ref seq
  ref_aa <- ref_aa[which(ref_aa != "-")] # removing gaps

  ## ------------------------ Working Sequence-Dataframes ------------------------ ##
  nr <- length(aln_aa$id)
  nc <- length(ref_aa)
  aa <- matrix(rep(NA, nr*nc), ncol = nc) # --- For amino acids
  rownames(aa) <- aln_aa$id
  colnames(aa) <- names(ref_aa)
  aa[1,] <- ref_aa

  ## ------------------------  Output Dataframes -------------------------------- ##
  output <- as.data.frame(matrix(rep(NA, 5 * ncol(aa)), ncol = 5))
  names(output) <- c('n', 'aa', 'H21', 'H4', 'site')
  output$n <- 1:nrow(output)
  output$aa <- ref_aa

  ## ---------------- Cicling through orthologous sequences --------------------- ##
  pos <- as.numeric(names(ref_aa))
  for (i in 2:nr){
    prot <- aln_aa$ali[i,]
    ngap <- which(prot != "-")
    prot <- prot[pos]
    aa[i,] <- prot
  }

  ## ----------- Computing H4 Shannon's Entropies  ------------------ ##
  charged <- c('E', 'D', 'H', 'K', 'R')
  hydrophobic <- c('L', 'I', 'V', 'M', 'F', 'W', 'Y', 'A')
  polar <- c('S', 'T', 'C', 'Q', 'N')
  special <- c('G','P','-')

  for (i in seq_len(nrow(output))){
    S <- 0
    X <- c()
    X <- c(X, sum(aa[,i] %in% charged)/sum(!is.na(aa[,i]))) # NA doesn't count
    X <- c(X, sum(aa[,i] %in% hydrophobic)/sum(!is.na(aa[,i])))
    X <- c(X, sum(aa[,i] %in% polar)/sum(!is.na(aa[,i])))
    X <- c(X, sum(aa[,i] %in% special)/sum(!is.na(aa[,i])))
    X <- X[which(X != 0)]
    S <- sum(X*log(1/X, base = 2))
    output$H4[i] <- round(S, 3)
  }

  ## ----------- Computing H21 Shannon's Entropies  ------------------ ##
  tablas <- apply(aa, 2, table)
  tablas <- lapply(tablas, as.data.frame)

    # For each amino acid in the sequence of reference, we have $`pos`. Example:
    # $`87`
    #     Var1 Freq
    # 1    A    1
    # 2    G    2
    # 3    S    1

  H <- function(x) { # For each position we compute the Shannon's entropy
      -sum( ( x[,2]/(sum( x[,2] ) ) )*log(x[,2]/(sum(x[,2])), base=2) )
  }

  entropy <- sapply(tablas, H) # Vector with the Shannon entropies for each position
  entropy <- sapply(entropy, function(x) round(x, digits=3))
  output$H21 <- entropy

  ## ------------------  Sorting sites  ----------------------- ##
  output$site[which(output$H21 < quantile(output$H21, th) & output$H4 < quantile(output$H4, th))] <- 'constrained'
  output$site[which(output$H21 > quantile(output$H21, th) & output$H4 < quantile(output$H4, th))] <- 'conservative'
  output$site[which(output$H21 > quantile(output$H21, (1- th)) & output$H4 > quantile(output$H4, (1- th)))] <- 'unconstrained'
  output$site[which(output$H21 < quantile(output$H21, th) & output$H4 > quantile(output$H4, th))] <- 'drastic'

  output$site[which(output$H21 == 0)] <- 'invariant'
  output$site[which(output$H21 > 0 & output$H4 == 0)] <- 'pseudo-invariant'

  ## ----------------  Plotting results  ---------------------- ##
  plot(output$H21, output$H4, xlab = 'H21', ylab = 'H4')
  abline(v = quantile(output$H21, th), lty = 2)
  abline(v = quantile(output$H21, 1-th), lty = 2)
  abline(h = quantile(output$H4, th), lty = 2)
  abline(h = quantile(output$H4, 1-th), lty = 2)

  points(output$H21[which(output$site == 'invariant')], output$H4[which(output$site == 'invariant')], col = 'blue', pch = 19)
  points(output$H21[which(output$site == 'pseudo-invariant')], output$H4[which(output$site == 'pseudo-invariant')], col = 'purple', pch = 19)
  points(output$H21[which(output$site == 'radical')], output$H4[which(output$site == 'radical')], col = 'red', pch = 19)
  points(output$H21[which(output$site == 'constrained')], output$H4[which(output$site == 'constrained')], col = 'cyan', pch = 19)
  points(output$H21[which(output$site == 'conservative')], output$H4[which(output$site == 'conservative')], col = 'aquamarine', pch = 19)
  points(output$H21[which(output$site == 'unconstrained')], output$H4[which(output$site == 'unconstrained')], col = 'orange', pch = 19)

  attr(output, 'units') <- 'bits'
  return(output)
}
