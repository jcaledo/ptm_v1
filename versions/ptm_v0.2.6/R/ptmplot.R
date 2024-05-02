## ---------- ptmplot.R ---------- ##
#                                    #
#     ptm.plot                       #
#     find.aaindex                   #
#                                    #
## -------------------------------- ##

## ------------------------------------------------------------------------------- ##
#   ptm.plot <- function(up_id, pdb="", property, ptm, window, sdata=F, ...)  #
## ------------------------------------------------------------------------------- ##
#' Plot Values of a Property and PTM Sites Along the Protein Sequence
#' @description Represents the values of a property and show the PTM sites along a protein sequence.
#' @usage ptm.plot(up_id, pdb = "", property, ptm, window = 1, sdata = FALSE, ...)
#' @param up_id a character string for the UniProt ID of the protein of interest.
#' @param pdb Optional argument to indicate the PDB and chain to be used (i.e. '1u8f.O'). If we leave this argument empty, the function will make the election for us whenever possible.
#' @param property a character string indicating the property of interest. It should be one of 'sasa', 'acc', 'dpx', 'eiip', 'volume', 'polarizability', 'avg.hyd', 'pi.hel', 'a.hel', 'b.sheet', 'B.factor', or 'own'.
#' @param ptm a character vector indicating the PTMs of interest. It should be among: 'ac' (acetylation), 'me' (methylation), 'meto' (sulfoxidation), 'p' (phosphorylation), 'ni' (nitration), 'su' (sumoylation) or 'ub' (ubiquitination), 'gl' (glycosylation), 'sni' (S-nitrosylation),'reg' (regulatory), 'dis' (disease).
#' @param window positive integer indicating the window size for smoothing with a sliding window average (default: 1, i.e. no smoothing).
#' @param sdata logical, if TRUE save a Rda file with the relevant data in the current directory.
#' @param ... when the user want to use his/her own amino acid index, it can be passed as a named vector.
#' @details If the property 'own' is selected, a named vector with the own index for the 20 amino acids should be passed as argument.
#' Currently the supported properties are:
#' \itemize{
#' \item{sasa:}  {Solvent-accessible surface area (3D)}
#' \item{acc:} {Accessibility (3D)}
#' \item{dpx:} {Depth (3D)}
#' \item{volume:} {Normalized van der Waals volume (1D)}
#' \item{mutability:} {Relative mutability, Jones 1992, (1D)}
#' \item{helix:} {Average relative probability of helix, Kanehisa-Tsong 1980,(1D)}
#' \item{beta-sheet:} {Average relative probability of beta-sheet, Kanehisa-Tsong 1980, (1D)}
#' \item{pi-helix:} {Propensity of amino acids within pi-helices, Fodje-Al-Karadaghi 2002, (1D)}
#' \item{hydropathy:} {Hydropathy index, Kyte-Doolittle 1982, (1D)}
#' \item{avg.hyd:} {Normalized average hydrophobicity scales, Cid et al 1992, (1D)}
#' \item{hplc:} {Retention coefficient in HPLC at pH7.4, Meek 1980, (1D)}
#' \item{argos:} {Hydrophobicity index, Argos et al 1982, (1D)}
#' \item{eiip:} {Electron-ion interaction potential, Veljkovic et al 1985, (1D)}
#' \item{polarizability:} {Polarizability parameter, Charton-Charton 1982, (1D)}
#' }
#' For 3D properties such as sasa, acc or dpx, for which different values can be obtained depending on the quaternary structure, we first compute the property values for each residue in the whole protein and plotted them against the residue position. Then, the value for this property is computed in the isolated chain (a single polypeptide chain) and in a second plot, the differences between the values in the whole protein and the chain are plotted against the residue position.
#' @return This function returns either one or two plots related to the chosen property along the primary structure, as well as the computed data if sdata has been set to TRUE.
#' @author Juan Carlos Aledo
#' @examples \dontrun{ptm.plot('P04406', property = 'sasa', window = 10, ptm = 'meto')}
#' @examples \dontrun{ptm.ptm('P04406', property = 'dpx', ptm = c('meto', 'p'))}
#' @seealso find.aaindex()
#' @importFrom graphics par
#' @importFrom graphics mtext
#' @importFrom graphics points
#' @export


ptm.plot <- function(up_id, pdb = "", property, ptm , window = 1, sdata = FALSE, ...){

  aa <- aai$aa

  ## -- Check if already exists a directory named 'split_chain'
  if (file.exists("./split_chain")){
    borrar <- FALSE
  } else {
    borrar <- TRUE
  }

  ## ------------------------------------------------------------------- ##
  ## ---------------- Is there any PDB for this protein? --------------- ##
  ## ------------------------------------------------------------------- ##
  if (nchar(pdb) == 6){ # ----- the user provide the PDB and chain
    t <- strsplit(pdb, split = "\\.")[[1]]
    pdb_id <- t[1]
    pdb_chain <- t[2]
    exists.pdb <- TRUE
  } else if (pdb == ""){ # ---- the PDB is selected by the script
    pdb <- pdb.select(up_id)
    if (is.null(pdb)){
      exists.pdb <- FALSE
    } else {
      exists.pdb <- TRUE
      pdb_id <- pdb[[1]][1]
      pdb_chain <- pdb[[2]][1]
      pdb_coverage <- attributes(pdb)$coverage
    }
  } else {
    message("Wrong pdb input!")
    return(NULL)
  }


  ## ------------------------------------------------------------------- ##
  ## ------------------- Use of cache when possible -------------------- ##
  ## ------------------------------------------------------------------- ##
  file_scan <- paste("./plotptm_cache/scan_", up_id, ".Rda", sep = "")
  if (file.exists(file_scan)){
    load(file_scan)
  }

  if (exists.pdb){
    file_sse <- paste("./plotptm_cache/sse_", pdb_id, ".Rda", sep = "")
    if (file.exists(file_sse)){
      load(file_sse)
    }
  }


  ## ------------------------------------------------------------------- ##
  ## ------- Checking that a suitable property has been selected ------- ##
  ## ------------------------------------------------------------------- ##
  pdb_property <- c('sasa', 'acc', 'dpx')
  aa_property <- c('volume','mutability', 'helix', 'beta-sheet', 'pi-helix',
                  'hydropathy', 'avg.hyd', 'hplc', 'argos', 'eiip',
                  'polarizability')

  names(aa_property) <- c('FAUJ880103', 'JOND920102', 'KANM800101', 'KANM800102',
                       'FODM020101', 'KYTJ820101', 'CIDH920105', 'MEEJ800101',
                       'ARGP820101', 'VELV850101', 'CHAM820101')

  # evo_property <- c('entropy7.aa', 'entropy100.aa','entropy7.condon', 'entropy100.codon')
  # all_property <- c(pdb_property, aa_property, evo_property, 'own')

  all_property <- c(pdb_property, aa_property,  'own')

  if (! property %in% all_property){
    message("A proper property must be indicated")
    return(NULL)
  }

  if (!exists.pdb & property %in% pdb_property){
    message('This property cannot be computed because no PDB file could be found')
    return(NULL)
  }

  if (!is.numeric(window) || window < 1){
    message("'window' must be numeric and positive")
    return(NULL)
  }

    if (property == 'own'){
      z <- list(...)
      if (length(z[[1]]) == 0){
        message("No aa index has been provided")
        return(NULL)
      } else {
        index <- z[[1]]
      }

      if (sum(aa == names(index)) != 20){
        message("The provided aa index must be a named numeric vector")
        return(NULL)
      }
    }

  ## ------------------------------------------------------------------- ##
  ## ---------- Checking that a suitable ptm has been selected --------- ##
  ## ------------------------------------------------------------------- ##
  supported_ptm <- c('ac', 'me', 'meto', 'p', 'su', 'ub', 'gl', 'sni',
                     'ni', 'reg', 'dis', 'all')
  names(supported_ptm) <- c('green', 'aquamarine', 'red', 'orange', 'deeppink',
                            'darkgreen', 'purple', 'yellow', 'blue',
                            'deepskyblue4', 'darkseagreen3', 'black')

  should_be_empty <- setdiff(ptm, supported_ptm)

  if (length(should_be_empty) != 0){
    message("The supported PTMs are 'ac', 'me', 'meto', 'p', 'su', 'ub', 'gl', 'sni',
                     'ni', 'reg', 'dis'")
    return(NULL)
  }

  if (ptm[1] == 'all'){
    ptm <- supported_ptm
  }


  ## ------------------------------------------------------------------- ##
  ## -------------- Scanning the protein for PTM sites ----------------- ##
  ## ------------------------------------------------------------------- ##
  if (! 'scan' %in% ls()){
    scan <- suppressWarnings(ptm.scan(up_id))
    if (is.null(scan)){
      message("Sorry, ptm.scan failed")
    }
    if (sdata){
      dir.create("plotptm_cache", showWarnings = FALSE)
      save(scan, file = paste("./plotptm_cache/scan_", up_id, ".Rda", sep = ""))
    }
  }
  if (! is.null(scan)){
    ptmScan <- scan[, c(2,3, which(colnames(scan) %in% ptm))]
    modifications <- !is.na(ptmScan[,3:dim(ptmScan)[2]])
    ptmScan$multi <- apply(as.matrix(modifications), 1, sum)
    ptmScan <- ptmScan[which(ptmScan$multi != 0), ]
    if (nrow(ptmScan) > 0){
      exists.ptm <- TRUE
      for (i in 1:nrow(ptmScan)){
        if (ptmScan$multi[i] == 1){
          ptmScan$col[i] <- names(ptmScan[which(ptmScan[i,-1] == TRUE) + 1][1])
        } else {
          ptmScan$col[i] <- "black"
        }
      }
      for (i in 1:nrow(ptmScan)){
        t <- ptmScan$col[i]
        if (t == "black"){
          color = "black"
        } else {
          color <- names(supported_ptm)[which(supported_ptm == t)]
        }
        ptmScan$col[i] <- color
      }
    } else { # no PTM sites found
      exists.ptm <- FALSE
    }
  } else {
    exists.ptm <- FALSE
  }


  ## ------------------------------------------------------------------- ##
  ## -------- Computing sse and raw 3d properties when required -------- ##
  ## ------------------------------------------------------------------- ##
  if (exists.pdb){
    if (! 'sse' %in% ls()){
      sse <- acc.dssp(pdb_id)
      if (is.null(sse)){
        message("Sorry, acc.dssp failed")
        return(NULL)
      }
      if (sdata){
        dir.create("plotptm_cache", showWarnings = FALSE)
        save(sse, file = paste("./plotptm_cache/sse_", pdb_id, ".Rda", sep = ""))
      }
    }
    sse <- sse[which(sse$chain == pdb_chain), ]
    seq <- sse$aa # aa sequence
    seq <- gsub("[a|b|c|d|e|f|g|h|i|j|k|l|m|n|o|p|q|r]", "C", seq) # half cystines
    names(seq) <- sse$respdb # position of each aa

    sse$sse <- NA
    for (i in 1:nrow(sse)){
      if (sse$ss[i] %in% c('G', 'H', 'I')){
        sse$sse[i] <- 1 # helixes arecoded as 1
      } else if (sse$ss[i] %in% c('E', 'B')){
        sse$sse[i] <- -1 # strands are coded as -1
      } else {
        sse$sse[i] <- 0 # coils are code as 0
      }
    }
    helices <- sse[which(sse$sse == 1),]
    strands <- sse[which(sse$sse == -1),]
    coils <- sse[which(sse$sse == 0),]

    if (window >= nrow(sse)) {
      message("'window' must be smaller than the sequence length")
      return(NULL)
    }

    if (property == 'dpx'){
      dpx <- res.dpx(pdb_id)
      if (is.null(dpx)){
        message("Sorry, res.dpx failed")
        return(NULL)
      }
    }


  } else { # When there is no pdb to be used

    seq <- get.seq(up_id, as.string = FALSE)[[1]] # aa sequence from Uniprot
    names(seq) <- 1:length(seq) # aa position

    if (window >= length(aa)){
      message("'window' must be smaller than the sequence length")
      return(NULL)
    }
  }

  ## ------------------------------------------------------------------- ##
  ## ---------- Computing the property values sequence ----------------- ##
  ## ------------------------------------------------------------------- ##
  if (property == 'own'){
    property_seq_mon <- property_seq_com <- index[seq]
  } else if (property %in% aa_property){
    i <- names(aa_property)[which(aa_property == property)]
    property_seq_mon <- property_seq_com <- bio3d::aa.index[[i]]$I[seq]
  } else if (property == "sasa"){
    property_seq_mon <- sse$sasa_chain
    property_seq_com <- sse$sasa_complex
  } else if (property == "acc"){
    property_seq_mon <- sse$acc_chain
    property_seq_com <- sse$acc_complex
  } else if (property == "dpx"){
    dpx <- res.dpx(pdb_id)[which(dpx$chain == pdb_chain), ]
    property_seq_mon <- dpx$min_dpx_chain
    property_seq_com <- dpx$min_dpx_complex
  }

  ## ------------------------------------------------------------------- ##
  ## ------- Smoothing the computed property values sequence ----------- ##
  ## ------------------------------------------------------------------- ##

  smooth <- function(x, window){
    if (window == 1) {
      y <- x
    } else {
      n <- length(x)
      y <- rep(NA, n)
      w <- ceiling(window/2)
    }
    if ((window%%2) == 0) {

      from <- w
      to <- n - w
      y[from:to] <- sapply(from:to, function(i) mean(x[(i - w + 1):(i + w)], na.rm = TRUE))
      if (from - 1 > 0) {
        y[1:(from - 1)] <- sapply(1:(from - 1), function(i) mean(x[1:(i + w)], na.rm = TRUE))
      }
      y[(to + 1):n] <- sapply((to + 1):n, function(i) mean(x[(i - w + 1):n], na.rm = TRUE))

    } else {

      from <- w
      to <- n - (w - 1)
      y[from:to] <- sapply(from:to, function(i) mean(x[(i - w + 1):(i + w - 1)], na.rm = TRUE))

      y[1:(from - 1)] <- sapply(1:(from - 1), function(i) mean(x[1:(i + w - 1)], na.rm = TRUE))

      y[(to + 1):n] <- sapply((to + 1):n, function(i) mean(x[(i - w + 1):n], na.rm = TRUE))

    }
    y <- round(y, 2)
    return(y)
  }

  if (window > 1){
    s_property_seq_mon <- smooth(property_seq_mon, window = window)
    s_property_seq_com <- smooth(property_seq_com, window = window)
  } else {
    s_property_seq_mon <- property_seq_mon
    s_property_seq_com <- property_seq_com
  }


  ## ---------------------------------------------------------------- ##
  ## --------------------- Dependent Variable  ---------------------- ##
  ## ---------------------------------------------------------------- ##

  if (property == 'sasa'){
    y <- s_property_seq_com
    ylab <- expression(paste('SASA (', ring(A)^2, ')', sep = ""))
    dy <- s_property_seq_mon - s_property_seq_com
    dylab <- expression(paste(Delta, 'SASA (', ring(A)^2, ')', sep = ""))
  } else if (property == 'acc'){
    y <- s_property_seq_com
    ylab <- 'Accessibility'
    dy <- s_property_seq_mon - s_property_seq_com
    dylab <- expression(paste(Delta, "Accessibility"))
  } else if (property == 'dpx'){
    y <- s_property_seq_com
    ylab <- expression(paste('Depth (', ring(A), ')', sep = ""))
    dy <- s_property_seq_com - s_property_seq_mon
    dylab <- expression(paste(Delta,'Depth (', ring(A), ')', sep = ""))
  } else if (property == 'eiip'){
    y <- s_property_seq_com
    ylab <- "EIIP"
  } else if (property == 'volume'){
    y <- s_property_seq_com
    ylab <- "Volume"
  } else if (property == 'avg.hyd'){
    y <- s_property_seq_com
    ylab <- "Averaged Hydrophobicity Index"
  } else if (property == 'pi-helix'){
    y <- s_property_seq_com
    ylab <- "pi-Helix Propensity"
  } else if (property == 'helix'){
    y <- s_property_seq_com
    ylab <- "Alpha-Helix Propensity"
  } else if (property == 'beta-sheet'){
    y <- s_property_seq_com
    ylab <- "Beta-Sheet Propensity"
  } else if (property == 'argos'){
    y <- s_property_seq_com
    ylab <- "Hydrophobicity index"
  } else if (property == 'mutability'){
    y <- s_property_seq_com
    ylab <- "Relative Mutability"
  } else if (property == 'hplc'){
    y <- s_property_seq_com
    ylab <- "Retention time in HPLC"
  } else if (property == "polarizability"){
    y <- s_property_seq_com
    ylab <- "Polarizability"
  } else if (property == 'own'){
    y <- s_property_seq_com
    ylab <- "Own Index"
  }

  ## ---------------------------------------------------------------- ##
  ## ------------------------ Ploting data -------------------------- ##
  ## ---------------------------------------------------------------- ##
  oldpar <- par(no.readonly = TRUE)
  on.exit(oldpar)
  par(mar = c(2, 4.1, 2, 2.1))
  ## -------- When two plots are shown (chain vs complex) ----------- ##
  if (sum(s_property_seq_com != s_property_seq_mon) != 0){ # Two plots
    # layout(matrix(c(1,2), 2, 1, byrow = TRUE))
    par(mfrow = c(2,1))
    xlab = ""
  } else {
    xlab = "Residue Number"
  }

  ## -------------- Plot 1: Property in the complex ----------------- ##
  x <- as.numeric(names(seq))
  y <- s_property_seq_com
  names(y) <- names(seq)

  ymin <- min(y, na.rm = TRUE)
  ymax <- max(y, na.rm = TRUE)
  ylim <- c(ymin - 0.1*(ymax - ymin), ymax + 0.1*(ymax - ymin))
  i <- which(!is.na(y))
  plot(x[i], na.omit(y), type = 'l', xlab = xlab,
       ylim = ylim, ylab = "")
  mtext(text = ylab, side = 2, line = 2)

  if (exists.ptm){ # Add point corresponding to the PTM sites
    ptmScan <- ptmScan[which(ptmScan$n %in% x), ] # Only modified sites present in the PDB
    if (nrow(ptmScan) >= 1) { # PTM sites may not be found in the PDB structure
      for (i in 1:nrow(ptmScan)){
        t <- ptmScan$n[i]
        points(t, y[which(names(y) == t)], pch = 19,
               cex = 0.55*ptmScan$multi[i], col = ptmScan$col[i])
      }
    }
    if (exists.pdb){
      points(coils$respdb, rep(ymin - 0.07*(ymax - ymin), nrow(coils)),
             pch = 15, col = 'pink', cex = 0.5)
      points(helices$respdb, rep(ymin -0.07*(ymax - ymin), nrow(helices)),
             pch = 15, col = 'cyan', cex = 0.5)
      points(strands$respdb, rep(ymin -0.07*(ymax - ymin), nrow(strands)),
             pch = 15, col = 'magenta', cex = 0.5)
    }
  }

  ## -------------- Plot 2: Property in the monomer ----------------- ##
  if (sum(s_property_seq_mon != s_property_seq_com) != 0){
    # Plot 2: Changes in the property |complex - chain|
    names(dy) <- names(y)
    dymin <- min(y, na.rm = TRUE)
    dymax <- max(y, na.rm = TRUE)
    dylim <- c(ymin - 0.1*(dymax - dymin), ymax + 0.1*(dymax - dymin))
    i <- which(!is.na(dy))
    plot(x[i], dy, type = 'l', xlab = 'Residue Number', ylab = "", ylim = dylim)
    mtext(text = dylab, side = 2, line = 2)

    if (exists.ptm){
      if (nrow(ptmScan) >= 1) {
        for (i in 1:nrow(ptmScan)){
          t <- ptmScan$n[i]
          points(t, dy[which(names(dy) == t)], pch = 19,
                 cex = 0.55*ptmScan$multi[i], col = ptmScan$col[i])
        }
      }
    }
    points(coils$respdb, rep(dymin - 0.07*(dymax - dymin), nrow(coils)),
           pch = 15, col = 'pink', cex = 0.5)
    points(helices$respdb, rep(dymin -0.07*(dymax - dymin), nrow(helices)),
           pch = 15, col = 'cyan', cex = 0.5)
    points(strands$respdb, rep(dymin -0.07*(dymax - dymin), nrow(strands)),
           pch = 15, col = 'magenta', cex = 0.5)
  }

  output <- "Work done."
  attr(output, "uniprot") <- up_id
  if (exists.pdb){
    attr(output, "pdb") <- pdb_id
    attr(output, "chain") <- pdb_chain
  } else {
    attr(output, "pdb") <- "No pdb found for this protein"
  }

  if (borrar & file.exists("./split_chain")){
    unlink("./split_chain", recursive = TRUE)
  }

  return(output)
}


## ------------------------------------------------------------------------------- ##
#                               find.aaindex(word)                                    #
## ------------------------------------------------------------------------------- ##
#' Find the Amino Acid Indexes
#' @description Finds amino acid indexes.
#' @usage find.aaindex(word)
#' @param word a character string for the key-word of interest.
#' @return The number and ID of the indexes matching the requested word
#' @author Juan Carlos Aledo
#' @examples find.aaindex('mutability')
#' @examples find.aaindex('Kyte-Doolittle')
#' @seealso ptm.plot()
#' @references https://www.genome.jp/aaindex/AAindex/list_of_indices
#' @export

find.aaindex <- function(word){
  aseveration <- function(x){
    is <- length(grep(word, x$D, ignore.case = TRUE)) != 0
    return(is)
  }
  ind <- which(sapply(bio3d::aa.index, aseveration))
  return(ind)
}
