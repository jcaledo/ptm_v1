## ---------------  ddG.R  ------------------- ##
#                                               #
#   imutant                                     #
#   foldx.mut                                   #
#   foldx.stab                                  #
#   foldx.assembly                              #
#   ddG.profile                                 #
#   ddG.ptm                                     #
#                                               #
## ------------------------------------------- ##


## ------------------------------------------------------------------ ##
#     imutant <- function(pdb, ch, pos, newres = "",                   #
#                             pH = 7, Te = 25, timeout = 60)           #
## ------------------------------------------------------------------ ##
#' Compute Changes in Stability (DDG)
#' @description Computes changes in the stability of a protein after a residue mutation using a machine-learning approach.
#' @usage imutant(protein, ch = "_", pos, newres = "", pH = 7, Te = 25, timeout = 60)
#' @param protein either the 4-letter identifier of a PDB structure, or the amino acid sequence (one letter amino acid code) of a protein.
#' @param ch a letter identifying the chain of interest.
#' @param pos the position, in the primary structure, of the residue to be mutated.
#' @param newres the one letter code of the residue to be incorporated. When a value is not entered for this parameter, then the function will compute DDG for the mutation to any possible amino acid.
#' @param pH a numeric value between 0 and 14.
#' @param Te a numeric value indicating the temperature in degrees Celsius.
#' @param timeout maximum time to wait, in seconds, for a response from the I-Mutant server.
#' @details This function implements the I-Mutant v2.0 tool, which is a fast method based on a support vector machine approach to predict protein stability changes upon single point mutations.
#' @return The function computes and returns a dataframe containing the following variables:
#' \itemize{
#' \item{Position:}  {Position in the primary structure of the mutated residue.}
#' \item{WT:}  {Amino acid found at that position in the wild-type protein.}
#' \item{NW:}  {New amino acid found in the mutated protein.}
#' \item{DDG:} {Change in Gibbs free energy (kcal/mol), defined as DDG = DGmt - DGwt, where DG is the change in Gibbs free energy for the folding of the protein from its unfolded state. Thus, a positive value means a stabilizing effect, and vice versa.}
#' \item{pH:}  {-log H+]}
#' \item{T:}   {Temperature in Celsius degrees.}
#' \item{RSA:}  {Relative Solvent Accessible Area (Only if a PDB file has been provided).}
#' }
#' @author Juan Carlos Aledo
#' @examples \dontrun{imutant(protein = '1u8f', ch = 'O', pos = 46, newres = 'K')}
#' @seealso foldx.mut(), ddG.profile()
#' @references Capriotti et al (2005) Nucl. Ac. Res. 33:W306-W310.
#' @importFrom httr POST
#' @importFrom httr http_error
#' @importFrom httr content
#' @importFrom xml2 read_html
#' @importFrom xml2 xml_text
#' @importFrom bio3d get.pdb
#' @export

imutant <- function(protein, ch = "_", pos, newres = "",
                              pH = 7, Te = 25, timeout = 60){


  aa <- c('A','R','N','D','C','Q','E','G','H','I',
          'L','K','M','F','P','S','T','W','Y','V')


  ## ----------------- Checking paramaters ---------------------- ##
  if (nchar(protein) == 4){
    is.pdb <- TRUE
    pdb <- paste(protein, ch, sep = ":")
    nc <- 7 # number variable of the output dataframe
  } else {
    is.pdb <- FALSE
    pdb <- "Protein Sequence"
    protein <- toupper(protein)
    nc <- 6 # number variable of the output dataframe
  }
  pH <- suppressWarnings(as.numeric(pH))
  if (!is.numeric(pH) | pH < 0 | pH > 14){
    stop("Please, provide a suitable pH value")
  }
  Te <- suppressWarnings(as.numeric(Te))
  if (!is.numeric(Te)){
    stop("Please, provide a suitable temperature")
  }
  if (!is.numeric(pos)){
    stop("The argument 'pos' must be a valid numeric position within the primary structure of the protein")
  } else if (! newres %in% aa){
    newres <- ""
  }


  ## ---------------- HTTP Form for I-Mutant -------------------- ##
  cgi_url <- "https://folding.biofold.org/cgi-bin/i-mutant2.0.cgi"
  query_parameters <- list("proteina" = protein,
                           "chain" = ch,
                           "posizione" = as.character(pos),
                           "newres" = newres,
                           "temp" = as.character(Te),
                           "ph" = as.character(pH),
                           "pred" = "ddg")

  resp <- tryCatch(
    {
      httr::POST(cgi_url, body = query_parameters)
    },
    error = function(cond){
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(resp)){
    message("Sorry, the server did not respond")
    return(NULL)
  } else {
    resp_xml <- xml2::read_html(resp)
  }

  ## ---------- Getting the URL for the results page ------------ ##
  a <- xml2::xml_find_all(resp_xml, xpath = "//a")
  b <- xml2::xml_text(a)
  url <- b[sapply(b, function(x) grepl("https://folding", x))]
  if (length(url) == 0){
    message("Sorry, no url for the response could be retrieved")
    return(NULL)
  }

  ## ------------------ Parsing the result page ----------------- ##
  waiting <- TRUE
  counter <- 0
  times <- ceiling(timeout/10)
  while (waiting & counter < times){
    resource <- httr::GET(url)
    check <- regexpr("Position {2,10}", resource)[[1]]
    if (check != -1){
      waiting <- FALSE
    } else {
      Sys.sleep(10)
      counter <- counter + 1
    }
  }

  if (counter >= times){
    message("Sorry timeout waiting the server's response")
    return(NULL)
  }

  resource_ <- httr::content(resource)
  text <- xml2::xml_text(resource_)
  t <- strsplit(text, "\n")[[1]]
  t <- t[which(t != "")]
  inicio <- which(grepl("Position", t))
  final <- which(grepl("WT: \\s", t))

  ## --------------------- Output dataframe -------------------- ##
  df <- as.data.frame(matrix(rep(NA, (final - inicio -1) * nc), ncol = nc))
  if (is.pdb){
    names(df) <- c("Position", "WT", "NEW", "DDG", "pH", "T", "RSA")
  } else {
    names(df) <- c("Position", "WT", "NEW", "DDG", "pH", "T")
  }
  for (i in seq_len(nrow(df))){
    fila <- strsplit(t[inicio + i], split = "\\s")[[1]]
    fila <- fila[which(fila != "")]
    df[i,] <- fila
  }
  # To fit our criterium:DDG = DGmut- DGwt of the unfolded <-> folded process:
  df$DDG <- -as.numeric(df$DDG)

  attr(df, "method") <- "I-Mutant v2.0"
  attr(df, "protein") <- pdb
  return(df)
}

## ----------------------------------------------------------------------------- ##
#  foldx.mut <- function(pdb, ch, pos, newres = "", pH = 7,                       #
#                   method = "buildmodel", keepfiles = FALSE)                     #
## ----------------------------------------------------------------------------- ##
#' Compute Changes in Stability (DDG)
#' @description Computes changes in the stability of a protein after a residue mutation using a force-field approach.
#' @usage foldx.mut(pdb, ch, pos, newres = "", pH =7, method = "buildmodel", keepfiles = FALSE)
#' @param pdb the 4-letter identifier of a PDB structure or the path to a PDB file.
#' @param ch a letter identifying the chain of interest.
#' @param pos the position, in the primary structure, of the residue to be mutated.
#' @param newres the one letter code of the residue to be incorporated. When a value is not entered for this parameter, then the function will compute DDG for the mutation to any possible amino acid (including phosphoserine, phosphothreonine, phosphotyrosine and hydroxiproline in the case of the 'positionscan' method).
#' @param pH a numeric value between 0 and 14.
#' @param method a character string specifying the approach to be used; it must be one of 'buildmodel', 'positionscan'.
#' @param keepfiles logical, when TRUE the repaired PDB file is saved in the working directory.
#' @details Two computational approaches for prediction of the effect of amino acid changes on protein stability are implemented. FoldX (buildmodel and positionscan methods) uses a force field approach and although it has been proved to be satisfactorily accurate, it is also a time-consuming method. An alternative much faster is I-Mutant, a method base on machine-learning
#' @return The function computes and returns the DDG (kcal/mol) for the requested residue change, defined as DDG = DGmt - DGwt, where DG is the Gibbs free energy for the folding of the protein from its unfolded state. Thus, a positive value means a destabilizing effect, and vice versa.
#' @author Juan Carlos Aledo
#' @examples \dontrun{foldx.mut('1aaq', 'A', 45, 'R')}
#' @seealso imutant(), ddG.profile()
#' @references Schymkowitz et al (2005) Nucl. Ac. Res. 33:W382-W388.
#' @importFrom bio3d get.pdb
#' @importFrom bio3d read.pdb
#' @importFrom bio3d aa321
#' @export

foldx.mut <- function(pdb, ch, pos, newres = "", pH = 7,
                      method = "buildmodel", keepfiles = FALSE){

  aa <- c('A','R','N','D','C','Q','E','G','H','I',
          'L','K','M','F','P','S','T','W','Y','V')

  ## ---------------------------- Checking parameters --------------------------------- ##
  war <- NULL
  if (!method %in% c("buildmodel", "positionscan")){
    method = "buildmodel"
    war <- "Method has been set to 'buildmodel'"
  }
  if (!is.numeric(pos)){
    stop("The argument 'pos' must be a valid numeric position within the primary structure of the protein")
  }

  ## --------------- Build the repaired pdb file if it doesn't exist ----------------- ##
  # Regardeless the FoldX method (either BuilModel or PositionScan), we will
  # start repairing the input PDB (see http://foldxsuite.crg.eu/command/RepairPDB).

  if (nchar(pdb) == 4){
    x <- tryCatch(
      {
        bio3d::get.pdb(pdb)
      },
      error = function(cond){
        return(NULL)
      },
      warning = function(w) conditionMessage(w)
    )
    if (is.null(x) | grepl("cannot", x)){
      message("Sorry, get.pdb failed")
      file.remove(paste(pdb, ".pdb", sep = ""))
      return(NULL)
    }
    path <- "./"
    id <- pdb
  } else {
    t <- strsplit(pdb, split = "/")[[1]]
    path <- paste(t[-length(t)], collapse = "/")
    id <- strsplit(t[length(t)], "\\.")[[1]][1]
    if (grepl("_Repair", id)){
      id <- strsplit(id, "_")[[1]][1]
    }
    oldwd <- getwd()
    on.exit(setwd(oldwd))
    setwd(path)
  }

  repaired = FALSE
  if (!grepl("_Repair.pdb", pdb)){
    repaired = TRUE
    A <- paste("foldx --command=RepairPDB --pdb=", pdb, ".pdb --ionStrength=0.05 ",  sep = "")
    B <- paste("--pH=", pH, " --water=CRYSTAL --vdwDesign=2 --pdbHydrogen=false", sep = "")
    AB <- paste(A, B)

    x <- tryCatch(
      {
        system(AB)
      },
      error = function(cond){
        return(NULL)
      },
      warning = function(w) conditionMessage(w)
    )
    if (is.null(x) | grepl("error", x) | x >0){
      message("Sorry, foldx failed")
      return(NULL)
    }
  }
  ## --- Cleaning unnecessary files generated by FoldX --- ##
  if (repaired){
    if (file.exists("Unrecognized_molecules.txt")){
      file.remove("Unrecognized_molecules.txt")
    }
  }

  if (method == 'buildmodel'){ ## ------------------------------ Method FoldX-BuildModel

    ## -------------- Formatting the mutations -------------- ##
    mypdb <- NULL
    if (file.exists(pdb)){
      mypdb <- bio3d::read.pdb(pdb)$atom
    } else {
      mypdb <- bio3d::read.pdb(paste(id, "_Repair.pdb", sep = ""))$atom
    }
    if (is.null(mypdb)){
      message("Sorry, ddG.mut failed")
      return(NULL)
    }

    wtres <- bio3d::aa321(mypdb$resid[which(mypdb$resno == pos & mypdb$chain == ch)[1]])

    con <- tryCatch(
      {
        file('individual_list.txt', 'w')
      },
      error = function(cond){
        return(NULL)
      }
    )
    if (is.null(con) | grepl("No such file", con)){
      message("Sorry, ddG.mut failed")
      return(NULL)
    }

    if (newres %in% aa){
      writeLines(paste(wtres, ch, pos, newres, ';', sep = ""), con)
    } else {
      for (newres in aa){
        if (newres != wtres){
          writeLines(paste(wtres, ch, pos, newres, ';', sep = ""), con)
        }
      }
      newres = ""
    }
    close(con)
    ## -------- Performing the computation ----------------- ##
    if (file.exists(pdb)){
      A <- paste("foldx --command=BuildModel --pdb=", pdb, sep = "")
    } else {
      A <- paste("foldx --command=BuildModel --pdb=", id, "_Repair.pdb", sep = "")
    }
    B <- paste(" --mutant-file=individual_list.txt --ionStrength=0.05 --pH=", pH,sep = "")
    C <- " --water=CRYSTAL --vdwDesign=2 --pdbHydrogens=false --numberOfRuns=3"
    ABC <- paste(A, B, C)

    x <- tryCatch(
      {
        system(ABC)
      },
      error = function(cond){
        return(NULL)
      },
      warning = function(w) conditionMessage(w)
    )
    if (is.null(x) | grepl("error", x) | x >0){
      message("Sorry, foldx failed")
      return(NULL)
    }


    ## -------------- Parsing the result page -------------- ##
    con <- tryCatch(
      {
        file(paste('Average_', id, "_Repair.fxout", sep = ""), 'r')
      },
      error = function(cond){
        return(NULL)
      }
    )
    if (is.null(con) | grepl("No such file", con)){
      message("Sorry, ddG.mut failed")
      return(NULL)
    }
    lines <- readLines(con)
    close(con)

    ## --------------- Output dataframe -------------------- ##
    if (newres %in% aa){
      df <- as.data.frame(matrix(rep(NA, 8), nrow = 1))
      colnames(df) <- c("Method", "PDB", "Chain", "Position", "WT", 'NEW', 'DDG', 'SD')
      l <- strsplit(lines[10], split = "\t")[[1]]
      df$Method[1] <- method
      df$PDB[1] <- id
      df$Chain[1] = ch
      df$Position[1] <- as.numeric(pos)
      df$WT[1] <- wtres
      df$NEW[1] <- newres
      df$DDG[1] <- as.numeric(l[3]) # We define DDG as DGmt - DGwt
      df$SD[1] <- as.numeric(l[2])

    } else if (newres == "") {
      df <- as.data.frame(matrix(rep(NA, 8*19), nrow = 19))
      colnames(df) <- c("Method", "PDB", "Chain", "Position", "WT", 'NEW', 'DDG', 'SD')
      for (i in 1:19){
        l <- strsplit(lines[9 + i], split = "\t")[[1]]
        df$Method[i] <- method
        df$PDB[i] <- pdb
        df$Chain[i] = ch
        df$Position[i] <- as.numeric(pos)
        df$WT[i] <- wtres
        df$NEW[i] <- setdiff(aa, wtres)[i]
        df$DDG[i] <- as.numeric(l[3]) # We define DDG as DGmt - DGwt
        df$SD[i] <- as.numeric(l[2])
      }
    } else {
      message("Something went wrong when parsing the result file")
      return(NULL)
    }

    ## --- Cleaning unnecessary files generated by FoldX --- ##
    file.remove(paste("Raw_", id, "_Repair.fxout", sep = ""))
    file.remove(paste("PdbList_", id, "_Repair.fxout", sep = ""))
    file.remove(paste("Dif_", id,  "_Repair.fxout", sep = ""))
    file.remove("individual_list.txt")
    file.remove(paste("Average_", id, "_Repair.fxout", sep = ""))
    file.remove("rotabase.txt")
    a <- paste(id, "_Repair_*.pdb", sep = "")
    b <- paste("WT_", id, "_Repair_*.pdb", sep = "")
    system(paste("rm ", a, sep = ""))
    system(paste("rm ", b, sep = ""))
    system("rmdir molecules")

    if (!keepfiles){
      if (file.exists(paste(id, "_Repair.pdb", sep = ""))){
        file.remove(paste(id, "_Repair.pdb", sep = ""))
      }
      if (file.exists(paste(id, ".pdb", sep = ""))){
        file.remove(paste(id, ".pdb", sep = ""))
      }
    }


    ## --------------------- Output --------------- ##
    if (!is.null(war)){
      warning(war)
    }
    closeAllConnections()
    return(df)

  } else if (method == 'positionscan'){ ## -------------------- Method FoldX-PositionScan

    nst_aa <- c('y', 'p','s', 'd', 'z', 'k', 'm', 'l', 'o', 'e', 'f')
    # The following non standard amino acids are also mutable:

    # Phosphotyrosine  ... PTR ... y ------ request as part of 'd'
    # Phosphothreonine.... TPO ... p ------ request as part of 'd'
    # Phosphoserine ...... SEP ... s ------ request as part of 'd'
    # Hydroxiproline ..... HYP ... h ------ request as part of 'd
    # Sulfotyrosine ...... TYS ... z
    # Monomethylated Lys . MLZ ... k
    # Dimethylated Lys ... MLY ... m
    # Trimethylated Lys .. M3L ... l
    # Charged ND1 His .... H1S ... o
    # Charged NE2 His .... H2S ... e
    # Neutral His ........ H3S ... f
    ## -------------- Formating the mutations -------------- ##
    mypdb <- NULL
    if (file.exists(pdb)){
      mypdb <- bio3d::read.pdb(pdb)$atom
    } else {
      mypdb <- bio3d::read.pdb(paste(id, "_Repair.pdb", sep = ""))$atom
    }
    if (is.null(mypdb)){
      message("Sorry, ddG.mut failed")
      return(NULL)
    }

    wtres <- bio3d::aa321(mypdb$resid[which(mypdb$resno == pos & mypdb$chain == ch)[1]])
    if (newres %in% c(aa, nst_aa)){
      positions <- paste(wtres, ch, pos, newres, sep = "")
    } else {
      positions <- paste(wtres, ch, pos, 'd', sep = "") # Includes pTyr, pSer, pThr and HO-Pro
    }
    ## -------- Performing the computation ----------------- ##
    A <- paste("foldx --command=PositionScan --pdb=", id, "_Repair.pdb", sep = "")
    B <- paste(" --positions=", positions, sep = "")
    AB <- paste(A, B)

    x <- tryCatch(
      {
        system(AB)
      },
      error = function(cond){
        return(NULL)
      },
      warning = function(w) conditionMessage(w)
    )
    if (is.null(x) | grepl("error", x) | x >0){
      message("Sorry, foldx failed")
      return(NULL)
    }

    ## -------------- Parsing the result page -------------- ##
    con <- tryCatch(
      {
        file(paste('PS_', id, "_Repair_scanning_output.txt", sep = ""), 'r')
      },
      error = function(cond){
        return(NULL)
      }
    )
    if (is.null(con) | grepl("No such file", con)){
      message("Sorry, ddG.mut failed")
      return(NULL)
    }

    lines <- readLines(con)
    close(con)
    ## --------------- Output dataframe -------------------- ##
    if (newres %in% c(aa, nst_aa)){
      df <- as.data.frame(matrix(rep(NA, 8), nrow = 1))
      colnames(df) <- c("Method", "PDB", "Chain", "Position", "WT", "NEW", "DDG", "SD")
      l <- strsplit(lines[2], split = "\t")[[1]]
      df$Method[1] <- method
      df$PDB[1] <- pdb
      df$Chain[1] = ch
      df$Position[1] <- as.numeric(pos)
      df$WT[1] <- wtres
      df$NEW[1] <- newres
      df$DDG[1] <- as.numeric(l[2]) # We define DDG as DGmt - DGwt
      df$SD[1] <- NA

    } else if (newres == "") {
      df <- as.data.frame(matrix(rep(NA, 8*24), nrow = 24))
      colnames(df) <- c("Method", "PDB", "Chain", "Position", "WT", 'NEW', 'DDG', 'SD')
      for (i in 1:24){
        l <- strsplit(lines[1 + i], split = "\t")[[1]]
        df$Method[i] <- method
        df$PDB[i] <- pdb
        df$Chain[i] = ch
        df$Position[i] <- as.numeric(pos)
        df$WT[i] <- wtres
        mt <- strsplit(l[1], split = "")[[1]]
        df$NEW[i] <- mt[length(mt)]
        df$DDG[i] <- round(as.numeric(l[2]), 2) # We define DDG as DGmt - DGwt
        df$SD[i] <- NA
      }
    } else {
      message("Something went wrong when parsing the result file")
      return(NULL)
    }

    if (!keepfiles){
      if (file.exists(paste(id, "_Repair.pdb", sep = ""))){
        file.remove(paste(id, "_Repair.pdb", sep = ""))
      }
      if (file.exists(paste(id, ".pdb", sep = ""))){
        file.remove(paste(id, ".pdb", sep = ""))
      }
    }
    ## --- Cleaning unnecessary files genetared by FoldX --- ##
    file.remove(paste("energies_", pos, "_", id, "_Repair.txt", sep = ""))
    # file.remove(paste("binding_energies_", pos, "_", id, "_Repair.txt", sep = ""))
    file.remove(paste("PS_", id, "_Repair_scanning_output.txt", sep = ""))
    file.remove(paste("PS_", id, "_Repair.fxout", sep = ""))
    file.remove("rotabase.txt")

    a <- paste("*_", id, "_Repair.pdb", sep = "")
    system(paste("rm ", a, sep = ""))
    system("rmdir molecules")
    ## --------------------- Output --------------- ##
    if (!is.null(war)){
      warning(war)
    }
    closeAllConnections()
    return(df)
  }
}


## ----------------------------------------------------------------------------- ##
#           foldx.stab <- function(pdb, pH = 7, I = 0.05)                         #
## ----------------------------------------------------------------------------- ##
#' Compute Folding Free Energy (DG)
#' @description Computes changes in the Gibbs free energy of the folding process of a protein.
#' @usage foldx.stab(pdb, pH = 7, I = 0.05)
#' @param pdb the 4-letter identifier of a PDB structure or the path to a PDB file.
#' @param pH a numeric value between 0 and 14.
#' @param I a value indicating the molar ionic strength.
#' @details This function implements the FoldX's command 'Stability'
#' @return The function computes and returns the DG (kcal/mol) of the folding process of the requested protein.
#' @author Juan Carlos Aledo
#' @examples \dontrun{foldx.stab('5zok')}
#' @seealso foldx.assembly()
#' @references Schymkowitz et al (2005) Nucl. Ac. Res. 33:W382-W388.
#' @importFrom bio3d get.pdb
#' @importFrom utils read.csv2
#' @export

foldx.stab <- function(pdb, pH = 7, I = 0.05) {

  ## --------- Getting the PDB ----------- ##
  if (nchar(pdb) == 4){
    x <- tryCatch(
      {
        bio3d::get.pdb(pdb)
      },
      error = function(cond){
        return(NULL)
      },
      warning = function(w) conditionMessage(w)
    )
    if (is.null(x) | grepl("cannot", x)){
      message("Sorry, get.pdb failed")
      file.remove(paste(pdb, ".pdb", sep = ""))
      return(NULL)
    }
    path <- "./"
    id <- pdb
  } else {
    t <- strsplit(pdb, split = "/")[[1]]
    path <- paste(t[-length(t)], collapse = "/")
    id <- strsplit(t[length(t)], "\\.")[[1]][1]
    if (grepl("_Repair", id)){
      id <- strsplit(id, "_")[[1]][1]
    }
  }

  ## ----- Repair pdb if neccesary -------- ##
  if (grepl("_Repair.pdb", pdb)){
    repaired = FALSE
  } else {
    repaired = TRUE
    A <- paste("foldx --command=RepairPDB --pdb=", id, ".pdb --ionStrength=", I,  sep = "")
    B <- paste(" --pH=", pH, " --water=CRYSTAL --vdwDesign=2 --pdbHydrogen=false", sep = "")
    AB <- paste(A, B)

    x <- tryCatch(
      {
        system(AB)
      },
      error = function(cond){
        return(NULL)
      },
      warning = function(w) conditionMessage(w)
    )
    if (is.null(x) | grepl("error", x) | x >0){
      message("Sorry, foldx failed")
      return(NULL)
    }

  }

  ## ------- Folding energy computation ------- ##
  A <- paste("foldx --command=Stability --pdb=", id, "_Repair.pdb", sep = "")
  B <- paste(" --pdb-dir=", path, " --output-dir=", path, sep = "")
  AB <- paste(A, B)

  x <- tryCatch(
    {
      system(AB)
    },
    error = function(cond){
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(x) | grepl("error", x) | x >0){
    message("Sorry, foldx failed")
    return(NULL)
  }

  ## ---------- Parsing the results ---------- ##
  DG <- read.csv2(file = paste(path, "/", id, "_Repair_0_ST.fxout", sep = ""), header = FALSE, sep = "\t")
  DG <- as.numeric(as.character(DG[1,2]))

  ## -------- Cleaning unnecesary files ------- ##
  f <- paste(path, "/", id, "_Repair_0_ST.fxout", sep = "")
  if (file.exists(f)){
    file.remove(f)
  }
  f <- paste(path, "/", id, "_Repair.fxout", sep = "")
  if (file.exists(f)){
    file.remove(f)
  }
  f <- 'molecules'
  if (file.exists(f)){
    file.remove(f)
  }
  if (file.exists('rotabase.txt')){
    file.remove('rotabase.txt')
  }

  ## --------- Output ----------- ##
  attr(DG, "pdb") <- pdb
  attr(DG, "units") <- "kcal/mol"
  closeAllConnections()
  return(DG)
}


# ------------------------------------------------------------------------------ ##
#         foldx.assembly <- function(pdb, mol1, mol2, pH = 7, I = 0.05)           #
## ----------------------------------------------------------------------------- ##
#' Compute Assembly Free Energy
#' @description Computes changes in the Gibbs free energy of the assembly process of a protein.
#' @usage foldx.assembly(pdb, mol1, mol2, pH = 7, I = 0.05)
#' @param pdb the 4-letter identifier of a PDB structure or the path to a PDB file.
#' @param mol1 molecule or group of molecules interacting with mol2 (see details)
#' @param mol2 molecule or group of molecules interacting with mol1 (see details)
#' @param pH a numeric value between 0 and 14.
#' @param I a value indicating the molar ionic strength.
#' @details This function implements the FoldX's command 'AnalyseComplex', which allows to determine the interaction energy between two molecules or two groups of molecules. For instance, if in a dimeric protein, formed by chain A and B, we may set: mol1 = 'A', mol2 = 'B'. If we are dealing with a trimer, we may set: mol1 = 'A', mol2: 'AB'.
#' @return The function returns a dataframe with the residues that make up the interface between mol1 and mol2, as well as the change in Gibbs free energy, DG, of the assembly process for the requested subunits.
#' @author Juan Carlos Aledo
#' @examples \dontrun{foldx.assembly(pdb = '1sev', mol1 = 'A', mol2 = 'B')}
#' @seealso foldx.stab()
#' @references Schymkowitz et al (2005) Nucl. Ac. Res. 33:W382-W388.
#' @importFrom bio3d get.pdb
#' @importFrom stats complete.cases
#' @export

foldx.assembly <- function(pdb, mol1, mol2, pH = 7, I = 0.05) {

  ## --------- Getting the PDB ----------- ##
  if (nchar(pdb) == 4){

    x <- tryCatch(
      {
        bio3d::get.pdb(pdb)
      },
      error = function(cond){
        return(NULL)
      },
      warning = function(w) conditionMessage(w)
    )
    if (is.null(x) | grepl("cannot", x)){
      message("Sorry, get.pdb failed")
      file.remove(paste(pdb, ".pdb", sep = ""))
      return(NULL)
    }
    path <- "./"
    id <- pdb
  } else {
    t <- strsplit(pdb, split = "/")[[1]]
    path <- paste(t[-length(t)], collapse = "/")
    id <- strsplit(t[length(t)], "\\.")[[1]][1]
    if (grepl("_Repair", id)){
      id <- strsplit(id, "_")[[1]][1]
    }
  }

  ## ----- Repair pdb if neccesary -------- ##
  if (file.exists(paste(path, "/", id, "_Repair.pdb", sep = ""))){
    repaired = FALSE
  } else {
    repaired = TRUE
    A <- paste("foldx --command=RepairPDB --pdb=", id, ".pdb --ionStrength=", I,  sep = "")
    B <- paste(" --pH=", pH, " --water=CRYSTAL --vdwDesign=2 --pdbHydrogen=false", sep = "")
    AB <- paste(A, B)
    x <- tryCatch(
      {
        system(AB)
      },
      error = function(cond){
        return(NULL)
      },
      warning = function(w) conditionMessage(w)
    )
    if (is.null(x) | grepl("error", x) | x >0){
      message("Sorry, foldx failed")
      return(NULL)
    }
  }

  ## ------ Complex Analysis --------------- ##
  A <- paste("foldx --command=AnalyseComplex --pdb=", id, "_Repair.pdb --analyseComplexChains=", mol1, ",", mol2, sep = "")
  B <- paste(" --pdb-dir=", path, " --output-dir=", path, sep = "")
  AB <- paste(A, B, sep = "")
  x <- tryCatch(
    {
      system(AB)
    },
    error = function(cond){
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(x) | grepl("error", x) | x >0){
    message("Sorry, foldx failed")
    return(NULL)
  }

  ## ---------- Parsing the results ---------- ##
  # a <- read.csv2(file = "./Indiv_energies_ABC_Repair_AC.fxout", skip = 8, sep = "\t", fill=TRUE)
  # b <- read.csv2(file = "./Interaction_ABC_Repair_AC.fxout", skip = 8, sep = "\t", fill=TRUE)
  ## Interface Residues:
  fxout <- paste(path, "/Interface_Residues_", id, "_Repair_AC.fxout", sep = "")
  c <- t(read.csv2(file = fxout, header = FALSE, skip = 10, sep = "\t", fill=TRUE))
  if (nrow(c) > 0){
    IR <- data.frame(
                     id = c[,1],
                     pos = apply(c, 2, function(x) substr(x, 3, nchar(x))),
                     aa = apply(c, 2, function(x) substr(x, 1, 1)),
                     chain = apply(c, 2, function(x) substr(x, 2, 2)),
                     stringsAsFactors = FALSE)
    IR <- IR[complete.cases(IR), ]
    rownames(IR) <- NULL
  } else {
    IR <- "No interface residues identified"
  }
  ## Assembly DG
  DGfile <- paste(path, "/Summary_", id, "_Repair_AC.fxout", sep = "")
  d <- read.csv2(file = DGfile, header = TRUE, skip = 8, sep = "\t")
  DG <- as.numeric(as.character(d$Interaction.Energy))

  ## -------- Cleaning unnecesary files ------- ##
  f <- paste(path, "/", id, "_Repair_0_ST.fxout", sep = "")
  if (file.exists(f)){
    file.remove(f)
  }
  f <- paste(path, "/", id, "_Repair.fxout", sep = "")
  if (file.exists(f)){
    file.remove(f)
  }
  f <- paste(path, "/", "Indiv_energies_", id, "_Repair_AC.fxout", sep = "" )
  if (file.exists(f)){
    file.remove(f)
  }
  f <- paste(path, "/", "Interaction_", id, "_Repair_AC.fxout", sep = "" )
  if (file.exists(f)){
    file.remove(f)
  }
  f <- paste(path, "/", "Interface_Residues_", id, "_Repair_AC.fxout", sep = "" )
  if (file.exists(f)){
    file.remove(f)
  }
  f <- paste(path, "/", "Summary_", id, "_Repair_AC.fxout", sep = "" )
  if (file.exists(f)){
    file.remove(f)
  }
  f <- 'molecules'
  if (file.exists(f)){
    file.remove(f)
  }
  if (file.exists('rotabase.txt')){
    file.remove('rotabase.txt')
  }


  ## --------- Output ----------- ##
  attr(IR, "Gibbs") <- DG
  attr(IR, "units") <- "kcal/mol"
  attr(IR, "pdb") <- pdb
  closeAllConnections()
  return(IR)
}

## ------------------------------------------------------------------ ##
#     ddG.profile <- function(prot, ch, pos, pH = 7, Te = 25)          #
## ------------------------------------------------------------------ ##
#' Contribution of a given position to changes in stability
#' @description Represents the sensitivity of a given position to changes in stability of a protein (DDG).
#' @usage ddG.profile(prot, ch, pos, pH = 7, Te = 25)
#' @param prot either the 4-letter identifier of a PDB structure, or the amino acid sequence (one letter amino acid code) of a protein.
#' @param ch a letter identifying the chain of interest.
#' @param pos the position, in the primary structure, of the residue to be mutated.
#' @param pH a numeric value between 0 and 14.
#' @param Te a numeric value indicating the temperature in degrees Celsius.
#' @details It must be remembered that DDG > 0 implies destabilizing change and DDG > 0 implies a stabilizing change.
#' @return The function returns a dataframe with the DDG values (kcal/mol) for each alternative amino acid, and a barplot grouping the amino acids according to their physicochemical nature.
#' @author Juan Carlos Aledo
#' @examples \dontrun{ddG.profile(prot = '1pga', ch = 'A', pos = 27)}
#' @seealso foldx.mut(), imutant()
#' @importFrom graphics barplot
#' @export

ddG.profile <- function(prot, ch, pos, pH = 7, Te = 25){

  acidic <- c('E', 'D')
  basic <- c('K', 'R', 'H')
  hydrophobic <- c('A', 'F', 'L', 'I', 'M', 'V', 'Y', 'W')
  polar <- c('N', 'Q', 'S', 'T', 'C')
  special <- c('G', 'P')

  ## ----------------- Checking paramaters ---------------------- ##
  war <- NULL
  if (nchar(prot) == 4){
    is.pdb <- TRUE
    pdb <- paste(prot, ch, sep = ":")
    nc <- 7 # number variable of the output dataframe
  } else {
    is.pdb <- FALSE
    pdb <- "Protein Sequence"
    protein <- toupper(prot)
    nc <- 6 # number variable of the output dataframe
  }
  pH <- suppressWarnings(as.numeric(pH))
  if (!is.numeric(pH) | pH < 0 | pH > 14){
    pH = 7
    war <- "pH has been set to 7.0"
  }
  Te <- suppressWarnings(as.numeric(Te))
  if (!is.numeric(Te)){
    Te = 25
    war <- "Temperature has been set to 25\u00baC"
  }
  if (!is.numeric(pos)){
    stop("The argument 'pos' must be a valid numeric position within the primary structure of the protein")
  }

  t <- tryCatch(
    {
      imutant(protein = prot, ch = ch, pos = pos, pH = pH, Te = Te)
    },
    error = function(cond){
      return(NULL)
    }
  )
  if (is.null(t)){
    message("Sorry, imutant failed")
    return(NULL)
  }


  ddG <- t$DDG
  names(ddG) <- t$NEW

  ddG_acidic <- ddG[which(names(ddG) %in% acidic)]
  ddG_basic <- ddG[which(names(ddG) %in% basic)]
  ddG_hyd <- ddG[which(names(ddG) %in% hydrophobic)]
  ddG_polar <- ddG[which(names(ddG) %in% polar)]
  ddG_special <- ddG[which(names(ddG) %in% special)]

  all <- c(ddG_acidic, ddG_basic, ddG_hyd, ddG_polar, ddG_special)
  df <- data.frame(ddG = all, aa = names(all))

  df$character <- NA
  df$character[which(df$aa %in% acidic)] <- 'acidi'
  df$character[which(df$aa %in% basic)] <- 'basic'
  df$character[which(df$aa %in% hydrophobic)] <- 'hydrophobic'
  df$character[which(df$aa %in% polar)] <- 'polar'
  df$character[which(df$aa %in% special)] <- 'special'

  df$color <- NA
  df$color[which(df$aa %in% acidic)] <- 1
  df$color[which(df$aa %in% basic)] <- 2
  df$color[which(df$aa %in% hydrophobic)] <- 3
  df$color[which(df$aa %in% polar)] <- 4
  df$color[which(df$aa %in% special)] <- 5



  mycol <- c("red", "blue", "orange", "green", "purple")
  barplot(height = df$ddG,
          names = df$aa,
          col =  mycol[df$color],
          horiz = TRUE, las = 1,
          xlab = expression(paste(Delta, Delta, "G (kcal/mol)", sep = "")),
          main = paste("Wild type: ", t$WT[1], pos, sep = ""))

  attr(df, "protein") <- prot
  attr(df, "wild_type") <- t$WT[1]
  attr(df, "position") <- pos
  attr(df, "pH") <- pH
  attr(df, "T_Celsius") <- Te

  return(df)
}


## ------------------------------------------------------------------ ##
#     ddG.ptm <- function(pdb, ch, pos, ptm, dir = 'f', pH = 7)        #
## ------------------------------------------------------------------ ##
#' PDB Model and Change in Stability of a Modified Protein
#' @description Builds a PDB model of the modified protein and computes the corresponding change in stability.
#' @usage ddG.ptm(pdb, ch, pos, ptm, dir = 'f', pH = 7)
#' @param pdb the 4-letter identifier of a PDB structure or the path to a PDB file.
#' @param ch a letter identifying the chain of interest.
#' @param pos the position, in the primary structure, of the residue to be modified.
#' @param ptm the post-translational modification to be considered. It should be one among: 'pSer', 'pThr', 'pTyr', 'MetO-Q', 'MetO-T'.
#' @param dir indicates the direction of the PTM reaction: either forward ('f'), or backward ('b').
#' @param pH a numeric value between 0 and 14.
#' @details The current function uses FoldX to build the model of the modified protein. Currently, FoldX does not allow to change Met by MetO, so we use glutamine (Q) or threonine (T) to mimic MetO.
#' @return The function computes and returns the DDG (kcal/mol) for the requested modification, defined as DDG = DGmodified - DGunmodified, where DG is the Gibbs free energy for the folding of the protein from its unfolded state. Thus, a positive value means a destabilizing effect, and vice versa. A PDB model containing the modified target is saved in the current directory.
#' @author Juan Carlos Aledo
#' @examples \dontrun{ddG.ptm('./1u8f_Repair.pdb', 'O', pos = 246, ptm = 'pThr')}
#' @seealso imutant(), foldx.mut(), ddG.profile()
#' @references Schymkowitz et al (2005) Nucl. Ac. Res. 33:W382-W388.
#' @importFrom bio3d get.pdb
#' @importFrom bio3d read.pdb
#' @importFrom bio3d aa321
#' @export

ddG.ptm <- function(pdb, ch, pos, ptm, dir = 'f', pH = 7){

  ## --------- Getting the PDB ----------- ##
  if (nchar(pdb) == 4){
    x <- tryCatch(
      {
        bio3d::get.pdb(pdb)
      },
      error = function(cond){
        return(NULL)
      },
      warning = function(w) conditionMessage(w)
    )
    if (is.null(x) | grepl("cannot", x)){
      message("Sorry, get.pdb failed")
      file.remove(paste(pdb, ".pdb", sep = ""))
      return(NULL)
    }
    path <- "./"
    id <- pdb
  } else {
    t <- strsplit(pdb, split = "/")[[1]]
    path <- paste(t[-length(t)], collapse = "/")
    id <- strsplit(t[length(t)], "\\.")[[1]][1]
    if (grepl("_Repair", id)){
      id <- strsplit(id, "_")[[1]][1]
    }
  }

  oldwd <- getwd()
  on.exit(setwd(oldwd))
  setwd(path)

  ## ------ Checking ptm argument --------- ##
  if (!ptm %in% c("pSer", "pThr", "pTyr", "MetO-Q", "MetO-T")){
    stop("Please, choose a proper PTM")
  }

  ## ----- Repair pdb if neccesary -------- ##
  if (file.exists(paste(id, "_Repair.pdb", sep = ""))){
    repaired = FALSE
  } else {
    repaired = TRUE
    A <- paste("foldx --command=RepairPDB --pdb=", id, ".pdb --ionStrength=0.05 ",  sep = "")
    B <- paste("--pH=", pH, " --water=CRYSTAL --vdwDesign=2 --pdbHydrogen=false", sep = "")
    AB <- paste(A, B)
    x <- tryCatch(
      {
        system(AB)
      },
      error = function(cond){
        return(NULL)
      },
      warning = function(w) conditionMessage(w)
    )
    if (is.null(x) | grepl("error", x) | x >0){
      message("Sorry, foldx failed")
      return(NULL)
    }

  }


  ## -------- Formating the Modification --------- ##
  mypdb <- bio3d::read.pdb(paste(id, "_Repair.pdb", sep = ""))$atom
  wtres <- bio3d::aa321(mypdb$resid[which(mypdb$resno == pos & mypdb$chain == ch)[1]])
  wt <- mypdb$resid[which(mypdb$resno == pos & mypdb$chain == ch)[1]]

  if (wt == 'SEP'){
    wtres <- 's'
  } else if (wt == "TPO"){
    wtres <- 'p'
  } else if (wt == 'PTR'){
    wtres <- 's'
  }

  war = NULL
  con <- file('individual_list.txt', 'w')

  if (dir == 'f'){ # Forward PTM
    if (ptm == "pSer"){
      if (wtres != "S") war = "Wild type residue is not Ser"
      writeLines(paste(wtres, ch, pos, newres = 's', ';', sep = ""), con)
    } else if (ptm == "pThr"){
      if (wtres != "T") war = "Wild type residue is not Thr"
      writeLines(paste(wtres, ch, pos, newres = 'p', ';', sep = ""), con)
    } else if (ptm == "pTyr"){
      if (wtres != "Y") war = "Wild type residue is not Tyr"
      writeLines(paste(wtres, ch, pos, newres = 's', ';', sep = ""), con)
    } else if (ptm == "MetO-Q"){
      if (wtres != "M") war = "Wild type residue is not Met"
      writeLines(paste(wtres, ch, pos, newres = 'Q', ';', sep = ""), con)
    } else if (ptm == "MetO-T"){
      if (wtres != "M") war = "Wild type residue is not Met"
      writeLines(paste(wtres, ch, pos, newres = 'T', ';', sep = ""), con)
    } else {
      stop("A suitable PTM modification must be provided!")
    }
  } else {# Backward PTM
    if (ptm == "pSer"){
      if (wt != "SEP"){
        war = "Wild type residue is not pSer"
      }
      writeLines(paste(wtres, ch, pos, newres = 'S', ';', sep = ""), con)
    } else if (ptm == "pThr"){
      if (wt != "TPO"){
        war = "Wild type residue is not pThr"
      }
      writeLines(paste(wtres, ch, pos, newres = 'T', ';', sep = ""), con)
    } else if (ptm == "pTyr"){
      if (wt != "TPR"){
        war = "Wild type residue is not pTyr"
      }
      writeLines(paste(wtres, ch, pos, newres = 'Y', ';', sep = ""), con)
    } else if (ptm == "MetO-Q"){
      if (wt != "GLN") war = "Wild type residue is not Gln"
      writeLines(paste(wtres, ch, pos, newres = 'M', ';', sep = ""), con)
    } else if (ptm == "MetO-T"){
      if (wt != "THR") war = "Wild type residue is not Thr"
      writeLines(paste(wtres, ch, pos, newres = 'M', ';', sep = ""), con)
    } else {
      stop("A suitable PTM modification must be provided!")
    }
  }

  close(con)

  ## -------- Build the PTM model --------- ##
  A <- paste("foldx --command=BuildModel --pdb=", id, "_Repair.pdb", sep = "")
  B <- " --mutant-file=individual_list.txt --ionStrength=0.05 --pH=7 --water=CRYSTAL"
  C <- " --vdwDesign=2 --pdbHydrogens=false --numberOfRuns=1"
  ABC <- paste(A,B,C, sep = "")
  x <- tryCatch(
    {
      system(ABC)
    },
    error = function(cond){
      return(NULL)
    },
    warning = function(w) conditionMessage(w)
  )
  if (is.null(x) | grepl("error", x) | x >0){
    message("Sorry, foldx failed")
    return(NULL)
  }


  ## -------------- Parsing the result page -------------- ##
  result_file <- paste('Average_', id, "_Repair.fxout", sep = "")
  con <- file(paste('Average_', id, "_Repair.fxout", sep = ""), 'r')
  lines <- readLines(con)
  close(con)
  DDG <- strsplit(lines[10], split = "\t")[[1]][3]

  ## --- Cleaning ancillary files genetared by FoldX --- ##
  f <- paste("Raw_", id, "_Repair.fxout", sep = "")
  if (file.exists(f)){
    file.remove(f)
  }
  f <- paste(id, "_Repair.fxout", sep = "")
  if (file.exists(f)){
    file.remove(f)
  }
  f <- paste("PdbList_", id, "_Repair.fxout", sep = "")
  if (file.exists(f)){
    file.remove(f)
  }
  f <- paste("Dif_", id,  "_Repair.fxout", sep = "")
  if (file.exists(f)){
    file.remove(f)
  }
  f <- "individual_list.txt"
  if (file.exists(f)){
    file.remove(f)
  }
  f <- paste("Average_", id, "_Repair.fxout", sep = "")
  if (file.exists(f)){
    file.remove(f)
  }
  f <- "rotabase.txt"
  if (file.exists(f)){
    file.remove(f)
  }
  f <- paste("WT_", id, "_Repair_*.pdb", sep = "")
  system(paste("rm ", f, sep = ""))
  system("rmdir molecules")


  ## --------- Output -------------- ##
  if (!is.null(war)){
    warning(war)
  }
  system(paste("mv ", id, "_Repair_1.pdb ", id, "_", ptm, pos, ".pdb", sep = ""))
  print(paste("THE FILE ", id, "_", ptm, pos, ".pdb", " HAS BEEN SAVED IN THE CURRENT DIRECTORY", sep = ""))

  output <- DDG
  attr(output, "wild-type") <- wt
  attr(output, "position") <- pos
  attr(output, "PTM") <- ptm
  attr(output, "units") <- "kcal/mol"

  closeAllConnections()
  return(output)

}

