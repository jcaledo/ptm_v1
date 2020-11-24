
# Q571F8 no existe PDB, no está en MetOSite
up_id <- 'O94925' # sí existe PDB, no está en MetOSite
up_id <- 'P01009' # sí existe PDB, sí está en MetOSite
up_id <- 'P0A392' # no existe PDB, sí está en MetOSite
up_id <- 'P0DP23'
## ------------------------------------------------------------------------------- ##
#   plot.serie <- function(up_id, property, ptm , sdata = FALSE)         #
## ------------------------------------------------------------------------------- ##
library(bio3d)
library(ptm)
rm(list=ls())
# get.pdb
# aa321
# PROPERTY: sasa, acc, dpx.min, dpx.avg, dpx.max, hetatm.min, hetatm.avg, 
# hetatm.max, ncontacts, ncontacts.intra, ncontacts.inter, entropy7.aa,
# entropy7.codon, entropy100.aa, entropy100.codon, aa01, eiip, volume, 
# polarizability, av.hyd, pi.hel, a.hel, b.sheet, B.factor, own

plot.serie <- function(up_id, property, ptm , sdata = FALSE, ...){

  ## ---------------- Is there any PDB for this protein? --------------- ##
  pdbs <- uniprot2pdb(up_id)
  if (nrow(pdbs) == 0){
    exists.pdb <- FALSE
  } else {
    exists.pdb <- TRUE
  }
  if (!exists.pdb & property %in% c('sasa', 'acc', 'dpx.min', 'dpx.avg',
                      'dpx.max', 'hetatm.min', 'hetatm.avg', 'hetatm.max',
                      'ncontacts', 'ncontacts.intra', 'ncontacts.inter')){
    stop('This property cannot be computed because no PDB file could be found')
  }
  
  ## ---------------- Is this protein found in MetOSite? --------------- ##
  metosite <- meto.scan(up_id)
  if (length(metosite$Metosites$met_pos) == 0){
    exists.meto <- FALSE
  } else {
    exists.meto <- TRUE
  }
  ## ------------------------------------------------------------------- ##
  ## ------------ Choose the PDB and chain when they exist ------------- ##
  ## ------------------------------------------------------------------- ##
  if (exists.meto & exists.pdb){
    
    # The chain showing the minimal distance to the metosite sequence
    metopdb <- meto.pdb(up_id, keepfile = TRUE)
    pdb_id <- metopdb[[4]]
    d <- metopdb[[1]]
    diag(d) <- NA
    d <- d[which(rownames(d) == 'metosite'),]
    chain <- names(d)[which(d == min(d, na.rm = TRUE))[1]]
    
  } else if (!exists.meto & exists.pdb){
    
    pdbs$length <-  pdbs$SP_END - pdbs$SP_BEG
    maximo <- which(pdbs$length == max(na.omit(pdbs$length)))
    pdb_id <- toupper(as.character(pdbs$PDB[maximo[1]]))
    chain <- as.character(pdbs$CHAIN[maximo[1]])
  }
  
  ## ------------------------------------------------------------------- ##
  ## --------------- Sequences (renumerated when necessary) ------------ ##
  ## ------------------------------------------------------------------- ##
  if (exists.meto & exists.pdb){
    df <- renum.meto(up_id)
    a <- renum.pdb(pdb_id, chain, up_id)
    df <- cbind(df[,1:4], a$pdb, df$meto_pos, a$pdb_pos)
    names(df)[3:7] <- c('aa_uni', 'aa_meto','aa_pdb', 'meto_pos', 'pdb_pos')
    df[is.na(df)] <- '-'
  } else if (exists.meto & !exists.pdb){
    df <- renum.meto(up_id)[,-6]
    names(df)[3:5] <- c('aa_uni', 'aa_meto', 'meto_pos')
  } else if (!exists.meto & exists.pdb){
    df <- renum.pdb(pdb_id, chain, up_id)
    df <- df[,-6]
    names(df)[3:5] <-  c('aa_uni', 'aa_pdb', 'pdb_pos')
  } else {
    seq <- get.seq(up_id, as.string = FALSE)[1,]
    df <- data.frame(uni_pos = 1:length(seq), uniprot = seq)
  }
  
  ## ------------------------------------------------------------------- ##
  ## --------------------- Computing sse and sasa ---------------------- ##
  ## ------------------------------------------------------------------- ##
  if (exists.pdb){
    df$chain <- chain
    df$ss <- NA
    df$sasa_complex <- df$sasa_chain <- NA
    df$acc_complex <- df$acc_chain <- NA
    
    asa <- acc.dssp(pdb_id, met = FALSE)
    asa <- asa[which(asa$chain == chain),]
  
  
    ## -------------- Matching residue numeration ------------------- ##
    k <- 1
    for (i in 1:nrow(asa)){
      for (j in k:nrow(df)){
        if (df$aa_pdb[j] == asa$aa[i]){
          df$chain[j] <- asa$chain[i]
          df$ss[j] <- asa$ss[i]
          df$sasa_chain[j] <- asa$sasa_chain[i]
          df$sasa_complex[j] <- asa$sasa_complex[i]
          df$acc_chain[j] <- asa$acc_chain[i]
          df$acc_complex[j] <- asa$acc_complex[i]
          k <- j + 1
          break
        }
      }
    }
  }
 ## ---------------------------------------------------------------- ##
 ## ----------- Computing other properties when needed ------------- ##
 ## ---------------------------------------------------------------- ##
  
 ## ------------------------- DPX ---------------------------------- ##
 if (property %in% c('dpx.min', 'dpx.avg', 'dpx.max') & exists.pdb){

    th <- strsplit(property, split = "\\.")[[1]][2]
    property <- strsplit(property, split = "\\.")[[1]][1]

    bio3d::get.pdb(pdb_id)
    path = paste("./", pdb_id, ".pdb", sep = "")
    dpx <- res.dpx(path, met = FALSE)
    dpx <- dpx[which(dpx$chain == chain),]
    df$min_dpx_complex <- df$min_dpx_chain <- NA
    df$avg_dpx_complex <- df$avg_dpx_chain <- NA
    df$max_dpx_complex <- df$max_dpx_chain <- NA
    k <- 1
    for (i in 1:nrow(dpx)){
      for (j in k:nrow(df)){
        if (df$aa_pdb[j] == dpx$resid[i]){
          df$min_dpx_complex[j] <- dpx$min_dpx_complex[i]
          df$min_dpx_chain[j] <- dpx$min_dpx_chain[i]
          df$avg_dpx_complex[j] <- dpx$avg_dpx_complex[i]
          df$avg_dpx_chain[j] <- dpx$avg_dpx_chain[i]
          df$max_dpx_complex[j] <- dpx$max_dpx_complex[i]
          df$max_dpx_chain[j] <- dpx$max_dpx_chain[i]
          k <- j+1
          break
        }
      }
    }
    
  ## -------------------------- HETATM ------------------------------ ##
 } else if (property %in% c('hetatm.min', 'hetatm.avg', 
                            'hetatm.max') & exists.pdb){

    th <- strsplit(property, split = "\\.")[[1]][2]
    property <- strsplit(property, split = "\\.")[[1]][1]

    bio3d::get.pdb(pdb_id)
    path = paste("./", pdb_id, ".pdb", sep = "")
    het <- search.hetatm(path)
    
    if (is.null(nrow(het))){
      df$het <- 'No het found'
    } else {
      
      for (i in 1:nrow(het)){
        t <- strsplit(het$res[i], split = '-')[[1]]
        het$aa[i] <- bio3d::aa321(t[1])
        het$resno[i] <- as.numeric(t[2])
        het$chain[i] <- t[3]
      }
      het <- het[which(het$chain == chain),]
      
      df$max_het <- df$avg_het <- df$min_het <- NA
      k <- 1
      for (i in 1:nrow(het)){
        for (j in k:nrow(df)){
          if (df$aa_pdb[j] == het$aa[i]){
            df$min_het[j] <- het$min[i]
            df$avg_het[j] <- het$avg[i]
            df$max_het[j] <- het$max[i]
            k <- j + 1
            break
          }
        }
      }
      file.remove(paste(pdb_id, "_atm.pdb", sep = ""))
      file.remove(paste(pdb_id, "_het.pdb", sep = ""))
    }
    
  ## --------------------------- CONTACTS ---------------------------- ##
  } else if (grepl('ncontacts', property)  & exists.pdb){
      cont <- contacts(pdb_id, threshold = 6)
      cont$residue <- aa321(cont$residue)
      
      df$inter_contacts <- df$intra_contacts <- NA
      k <- 1
      for (i in 1:nrow(cont)){
        for (j in k:nrow(df)){
          if (df$aa_pdb[j] == cont$residue[i]){
            df$intra_contacts[j] <- cont$intra_contacts[i]
            df$inter_contacts[j] <- cont$inter_contacts[i]
            k <- j + 1
            break
          }
        }
      }
      
  ## ------------------------- SHANNON ENTROPY------------------------- ##
  } else if (property == 'entropy7'){
    sp <- species.mapping(up_id)
    kegg <- id.mapping(up_id, from = 'uniprot', to = 'kegg', sp)
    ent <- shannon(kegg, 'seven')
    df$entropy7_aa <- df$entropy7_codon <- NA
    k <- 1
    for (i in 1:nrow(ent)){
      for (j in k:nrow(df)){
        if (df$aa_uni[j] == ent$aa[i]){
          df$entropy7_codon[j] <- ent$Hcodon[i]
          df$entropy7_aa[j] <- ent$Haa[i]
          k <- j + 1
          break
        }
      }
    }
  } else if (property == 'entropy100'){
    sp <- species.mapping(up_id)
    kegg <- id.mapping(up_id, from = 'uniprot', to = 'kegg', sp)
    ent <- shannon(kegg, 'one-hundred')
    df$entropy100_aa <- df$entropy100_codon <- NA
    k <- 1
    for (i in 1:nrow(ent)){
      for (j in k:nrow(df)){
        if (df$aa_uni[j] == ent$aa[i]){
          df$entropy100_codon[j] <- ent$Hcodon[i]
          df$entropy100_aa[j] <- ent$Haa[i]
          k <- j + 1
          break
        }
      }
    }
 
  ## ------------- AMINO ACID PRESENCE/ABSENCE ---------------------- ##
  } else if (grepl('01', property)){
    aa <- substring(property, 1, 1)
    presence <- rep(0, length(df$aa_uni))
    presence[which(df$aa_uni == aa)] <- 1
    df$presence <- presence
  
  ## --------------------- AMINO ACID INDEX ------------------------- ##
  } else if (property == 'eiip'){
    eiip <- ptm:::aai$eiip
    names(eiip) <- ptm:::aai$aa
    df$eiip <- sapply(df$aa_uni, 
                      function(x) eiip[which(names(eiip) == x)])
    
  } else if (property == 'volume'){
    volume <- ptm:::aai$volume
    names(volume) <- ptm:::aai$aa
    df$volume <- sapply(df$aa_uni, 
                      function(x) volume[which(names(volume) == x)])
  } else if (property == 'av.hyd'){
    av_hyd <- ptm:::aai$av_hyd
    names(av_hyd) <- ptm:::aai$aa
    df$av.hyd <- sapply(df$aa_uni, 
                        function(x) av.hyd[which(names(av_hyd) == x)])
  } else if (property == 'pi.hel'){
    pi_hel <- ptm:::aai$pi_hel
    names(pi_hel) <- ptm:::aai$aa
    df$pi_hel <- sapply(df$aa_uni, 
                        function(x) pi_hel[which(names(pi_hel) == x)])
  } else if (property == 'a.hel'){
    a_hel <- ptm:::aai$a_hel
    names(a_hel) <- ptm:::aai$aa
    df$a_hel <- sapply(df$aa_uni, 
                        function(x) a_hel[which(names(a_hel) == x)])
  } else if (property == 'b.sheet'){
    b_sheet <- ptm:::aai$b_sheet1
    names(b_sheet) <- ptm:::aai$aa
    df$b_sheet <- sapply(df$aa_uni, 
                        function(x) b_sheet[which(names(b_sheet) == x)])
  } else if (property == 'B.factor'){
    B_factor <- ptm:::aai$B_factor
    names(B_factor) <- ptm:::aai$aa
    df$B_factor <- sapply(df$aa_uni, 
                         function(x) B_factor[which(names(B_factor) == x)])
  }
  ## ----------------- OWN AMINO ACID INDEX ------------------------- ##
  z <- list(...)
  if (!is.null(z$own)){
    df$own <- sapply(df$aa_uni, 
                     function(x) z$own[which(names(z$own) == x)])
  }
  
  ## ---------------------------------------------------------------- ##
  ## -------------- Dependent Variable and Data File ---------------- ##
  ## ---------------------------------------------------------------- ##

  if (property == 'sasa'){
    y <- df$sasa_complex
    ylab <- expression(paste('SASA (', ring(A)^2, ')', sep = ""))
    dy <- df$sasa_chain - df$sasa_complex
    dylab <- expression(paste(Delta, 'SASA (', ring(A)^2, ')', sep = ""))
  } else if (property == 'acc'){
    y <- df$acc_complex
    ylab <- 'Accessibility'
    dy <- df$acc_chain - df$acc_complex
    dylab <- expression(paste(Delta, "Accessibility"))
  } else if (property == 'dpx'){
    if (th == 'min'){
      y <- df$min_dpx_complex
      ylab <- expression(paste('Minimal Depth (', ring(A), ')', sep = ""))
      dy <- df$min_dpx_complex - df$min_dpx_chain
      dylab <- expression(paste(Delta,'Minimal Depth (', ring(A), ')', sep = ""))
    } else if (th == 'max'){
      y <- df$max_dpx_complex
      ylab <- expression(paste('Maximal Depth (', ring(A), ')', sep = ""))
      dy <- df$max_dpx_complex - df$max_dpx_chain
      dylab <- expression(paste(Delta,'Maximal Depth (', ring(A), ')', sep = ""))
    } else {
      y <- df$avg_dpx_complex
      ylab <- expression(paste('Averaged Depth (', ring(A), ')', sep = ""))
      dy <- df$avg_dpx_complex - df$avg_dpx_chain
      dylab <- expression(paste(Delta,'Averaged Depth (', ring(A), ')', sep = ""))
    }
  } else if (property == 'hetatm'){
    if (th == 'min'){
      y <- df$min_het
      ylab <- expression(paste('Minimal Distance to AS (', ring(A), ')', sep = ""))
    } else if (th == 'max'){
      y <- df$max_het
      ylab <- expression(paste('Maximal Distance to AS (', ring(A), ')', sep = ""))
    } else {
      y <- df$avg_het
      ylab <- expression(paste('Averaged Distance to AS (', ring(A), ')', sep = ""))
    }
  } else if (property == 'ncontacts'){
    y <- df$intra_contacts + df$inter_contacts
    ylab <-'Total number of contacts'
  } else if (property == 'ncontacts.intra'){
    y <- df$intra_contacts 
    ylab <-'Number of intramolecular contacts'
  } else if (property == 'ncontacts.inter'){
    y <- df$inter_contacts 
    ylab <-'Number of intermolecular contacts'
  } else if (property == 'entropy7.aa'){
    y <- df$entropy7_aa
    ylab <- "Shannon's entropy-aa7"
  } else if (property == 'entropy7.codon'){
    y <- df$entropy7_codon
    ylab <- "Shannon's entropy-codon7"
  } else if (property == 'entropy100.aa'){
    y <- df$entropy100_aa
    ylab <- "Shannon's entropy-aa100"
  } else if (property == 'entropy100.codon'){
    y <- df$entropy100_codon
    ylab <- "Shannon's entropy-codon100"
  } else if (grepl('01', property)){
    y <- df$presence
    ylab <- paste('Presence/Absence of', aa, sep = " ")
  } else if (property == 'eiip'){
    y <- df$eiip
    ylab <- "EIIP"
  } else if (property == 'volume'){
    y <- df$volume
    ylab <- "Volume"
  } else if (property == 'av.hyd'){
    y <- df$av_hyd
    ylab <- "Averaged Hydrophobicity Index"
  } else if (property == 'pi.hel'){
    y <- df$pi_hel
    ylab <- "pi-Helix Propensity"
  } else if (property == 'a.hel'){
    y <- df$a_helix
    ylab <- "alpha-Helix Propensity"
  } else if (property == 'b.sheet'){
    y <- df$b_sheet
    ylab <- "beta-Sheet Propensity"
  } else if (property == 'B.factor'){
    y <- df$B_factor
    ylab <- "Averaged B-Factor"
  } else if (property == 'own'){
    y <- df$own
    ylab <- "Own Index"
  } 
  
  ## -------------- Save data if required ------------------------ ##
  if (sdata == TRUE){
    name <- paste(up_id, '_plot_ptm.Rda', sep = "")
    save(df, file = paste("./", name, sep = ""))
  }

  ## ---------------------------------------------------------------- ##
  ## -------------------- The PTMs of interest ---------------------- ##
  ## ---------------------------------------------------------------- ##
  ac_at <- me_at <- meto_at <- p_at <- numeric()
  su_at <- ub_at <- reg_at <- numeric()
  
  # --------------------- ACETYLATION ------------------------------- ##
  if ('ac' %in% ptm){
    ac <- ac.scan(up_id)
    ac <- unique(ac$modification)
    ac <- strsplit(ac, split = "-")
    ac_at <- as.numeric(sapply(ac, function(x) substring(x, 2))[1,])
  }
  # ---------------------- METHYLATION  ----------------------------- ##
  if ('me' %in% ptm){ 
    me <- me.scan(up_id)
    me <- unique(me$modification)
    me <- strsplit(me, split = "-")
    me_at <- as.numeric(sapply(me, function(x) substring(x, 2))[1,])
    me_at <- as.numeric(unique(me_at))
  }
  ## -------------------- SULFOXIDATION ----------------------------- ##
  if ('meto' %in% ptm){ 
    meto <- meto.scan(up_id, report = 2)
    meto_at <- as.numeric(meto$Metosites[,1])
    # Data quality check:
    a <- sum(df$aa_meto[which(df$meto_pos %in% meto_at)] == 'M')
    if (a != length(meto_at)){
      stop('Check MetO site positions')
    }
  }
  ## ------------------- PHOSPHORYLATION --------------------------- ##
  if ('p' %in% ptm){
    p <- p.scan(up_id)
    p <- unique(p$modification)
    p <- strsplit(p, split = "-")
    p_at <- as.numeric(sapply(p, function(x) substring(x, 2))[1,])
    p_at <- unique(p_at)
  } 
  ## -------------------- SUMOYLATION ------------------------------ ##
  if ('su' %in% ptm){
    su <- su.scan(up_id)
    su <- unique(su$modification)
    su <- strsplit(su, split = "-")
    su_at <- as.numeric(sapply(su, function(x) substring(x, 2))[1,])
    su_at <- unique(su_at)
  } 
  ## ------------------- UBIQUITINATION ---------------------------- ##
  if ('ub' %in% ptm){
    ub.scan(up_id)
    ub <- unique(ub$modification)
    ub <- strsplit(ub, split = "-")
    ub_at <- as.numeric(sapply(ub, function(x) substring(x, 2))[1,])
    ub_at <- unique(ub_at)
  }
  ## --------------------- REGULATORY ------------------------------ ##
  if ('reg' %in% ptm){
    reg.scan(up_id)
    reg <- unique(reg$modification)
    reg <- strsplit(reg, split = "-")
    reg_at <- as.numeric(sapply(reg, function(x) substring(x, 2))[1,])
    reg_at <- unique(reg_at)
  }
  
  ## ---------------------------------------------------------------- ##
  ## -------------------- SECUNDARY STRUCTURE ----------------------- ##
  ## ---------------------------------------------------------------- ##
  if (exists.pdb){
    df$sse <- NA
    for (i in 1:nrow(df)){
      if (df$ss[i] %in% c('G', 'H', 'I')){
          df$sse[i] <- 1 # helixes arecoded as 1
        } else if (df$ss[i] %in% c('E', 'B')){
          df$sse[i] <- -1 # strands are coded as -1
        } else {
          df$sse[i] <- 0 # coils are code as 0
        }
      }
      helices <- df[which(df$sse == 1),]
      strands <- df[which(df$sse == -1),]
      coils <- df[which(df$sse == 0),]
  }
  
  ## ---------------------------------------------------------------- ##
  ## ----------------------- PLOTTING DATA -------------------------- ##
  ## ---------------------------------------------------------------- ##
  if (property %in% c('sasa', 'acc', 'dpx')){ # Two plots to be shown
    layout(matrix(c(1,2), 2, 1, byrow = TRUE))
  }

  # pdf(paste(up_id, "-", pdb_id, ".pdf", sep =""),
  #     width = 11.69, height = 4.13)

  # Plot 1: Property in the complex
  ymin <- min(y, na.rm = TRUE)
  ymax <- max(y, na.rm = TRUE)
  ylim <- c(ymin - 0.1*(ymax - ymin), ymax + 0.1*(ymax - ymin))
  i <- which(!is.na(y))
  plot(df$uni_pos[i], na.omit(y), type = 'l', xlab = 'Residue Number',
       ylim = ylim, ylab = "")
  mtext(text = ylab, side = 2, line = 2)
  # points(df$resno[which(df$aa_chain == 'M')],
  #        y[which(df$aa_chain == 'M')], pch = 19, cex =1)
  points(df$uni_pos[which(df$meto_pos %in% meto_at)],
         y[which(df$meto_pos %in% meto_at)], 
         pch = 19, cex = 1, col = 'darkorange')
  points(df$uni_pos[ac_at], y[ac_at], pch = 19, cex = 1, col = 'blue')
  points(df$uni_pos[me_at], y[me_at], pch = 19, cex = 1, col = 'darkgreen')
  points(df$uni_pos[p_at], y[p_at], pch = 19, cex = 1, col = 'red')
  points(df$uni_pos[su_at], y[su_at], pch = 19, cex = 1, col = 'aquamarine')
  points(df$uni_pos[ub_at], y[ub_at], pch = 19, cex = 1, col = 'deeppink')
  points(df$uni_pos[reg_at], y[reg_at], pch = 19, cex = 1)
  
  points(coils$uni_pos, rep(ymin - 0.07*(ymax - ymin), nrow(coils)),
         pch = 15, col = 'pink', cex = 0.5)
  points(helices$uni_pos, rep(ymin -0.07*(ymax - ymin), nrow(helices)),
         pch = 15, col = 'cyan', cex = 0.5)
  points(strands$uni_pos, rep(ymin -0.07*(ymax - ymin), nrow(strands)),
         pch = 15, col = 'magenta', cex = 0.5)

  if (property %in% c('sasa', 'acc', 'dpx')){
    # Plot 2: Change in the property |complex - chain|
    dymin <- min(dy, na.rm = TRUE)
    dymax <- max(dy, na.rm = TRUE)
    dylim <- c(dymin - 0.1*(dymax - dymin), dymax + 0.1*(dymax - dymin))
    plot(df$uni_pos[i], dy, type = 'l', xlab = 'Residue Number',
         ylab = "", ylim = dylim)
    mtext(text = dylab, side = 2, line = 2)
    # points(df$resno[which(df$aa_chain == 'M')],
    #        dy[which(df$aa_chain == 'M')], pch = 19, cex =1)
    points(meto_at, dy[meto_at], pch = 19, cex = 1, col = 'darkorange')
    points(df$uni_pos[ac_at], dy[ac_at], pch = 19, cex = 1, col = 'blue')
    points(df$uni_pos[me_at], dy[me_at], pch = 19, cex = 1, col = 'darkgreen')
    points(df$uni_pos[p_at], dy[p_at], pch = 19, cex = 1, col = 'red')
    points(df$uni_pos[su_at], dy[su_at], pch = 19, cex = 1, col = 'aquamarine')
    points(df$uni_pos[ub_at], dy[ub_at], pch = 19, cex = 1, col = 'deeppink')
    points(df$uni_pos[reg_at], y[reg_at], pch = 19, cex = 1)
    
    points(coils$uni_pos, rep(dymin - 0.07*(dymax - dymin), nrow(coils)),
           pch = 15, col = 'pink', cex = 0.5)
    points(helices$uni_pos, rep(dymin -0.07*(dymax - dymin), nrow(helices)),
           pch = 15, col = 'cyan', cex = 0.5)
    points(strands$uni_pos, rep(dymin -0.07*(dymax - dymin), nrow(strands)),
           pch = 15, col = 'magenta', cex = 0.5)
  }
  return(df)
  #   # dev.off()
}


