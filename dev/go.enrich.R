## ---------------------------------------------------------------- ##
#  go.enrich <- function(target, background, aspect = 'BP', n = 20)  #
## ---------------------------------------------------------------- ##
#' GO Terms Enrichment Tests
#' @description Carry out GO term enrichment tests
#' @usage go.enrich(target, background, aspect = 'BP', n = 20)
#' @param target either a vector containing the UniProt IDs of the set of interest, or the path to the txt file containing the list of identifiers for the sample of interest.
#' @param background  a dataframe with two columns (Uniprot ID and GO terms) and as many rows as different proteins there are in the background set.
#' @param aspect  character string indicating the aspect or sub-ontology. It must be one of 'BP' (Biological Process), 'MF' (Molecular Function) or 'CC' (Cellular Component). acting as background.
#' @param n maximum number of enriched GO terms reported.
#' @details It is essential that the items in the 'sample' vector correspond to items within the background.
#' @return Returns the results of the enrichement test as a dataframe.
#' @author Juan Carlos Aledo
#' @references Rhee et al. (2008) Nature Reviews Genetics 9:509–515.
#' @seealso search.go(), term.go(), background.go(), get.go(), gorilla(), net.go()
#' @examples \dontrun{go.enrich("../bench/sample.txt", "../bench/background.txt", 'CC', n = 10)}
#' @importFrom topGO readMappings
#' @importFrom topGO runTest
#' @importFrom topGO GenTable
#' @importFrom topGO annFUN.gene2GO
#' @importFrom topGO topGO.Env
#' @export

go.enrich <- function(target, background, aspect = 'BP', n = 20){

  ## ----- The target sample to be analyzed
  if (is.character(target) & length(target) == 1){ # input as path to the txt
    if (gregexpr('txt', target)[[1]] != -1){
      sample <- read.csv(target, header = FALSE)
      sample <- trimws(as.character(sample$V1))
    } else {
      stop("A proper path to a txt file should be provided for the target set")
    }
  } else if (is.character(target) & length(target) > 1){ # input as vector
    sample <- trimws(as.character(target))
  } else if (is.data.frame(target) & nrow(target) > 1){ # input as dataframe
    sample <- trimws(as.character(target))
  } else {
    stop("A proper target set must be provided")
  }

  ## ----- Check the input background set
  if (is.data.frame(background) & ncol(background) == 2){
    bg <- trimws(as.character(background[,1]))
  } else {
    stop("A proper background set must be provided")
  }

  ## ----- Check that the target is included into the background set
  sample_bg <- intersect(sample, bg)
  if (length(sample) != sum(sample_bg == sample)){
    stop("Please, make sure that all the target proteins are contained in the background set")
  }

  ## ----- Getting GO ids for the backgraund set
  bg <- data.frame(up_id = bg, GO_id = background[,2])

  bg_proteins <- bg$up_id
  write.table(bg, file = "file_temp.map", quote = FALSE,
              sep = "\t", row.names = FALSE, col.names = FALSE)
  bg2GO <- topGO::readMappings(file = 'file_temp.map')
  bg_proteins <- names(bg2GO)
  file.remove("file_temp.map")

  ## ------- Compare sample vs bg_proteins
  # It is essential that the items in the 'sample' vector
  # correspond to items within the background ie 'bg_proteins'
  compared_proteins <- factor(as.integer(bg_proteins %in% sample))
  names(compared_proteins) <- bg_proteins

  ## ------- Create topGO object
  # ann <- function (whichOnto, feasibleGenes = NULL, gene2GO)
  # {
  #   ontoGO <- get(paste("GO", whichOnto, "Term", sep = ""))
  #   if (!is.null(feasibleGenes))
  #     gene2GO <- gene2GO[intersect(names(gene2GO), feasibleGenes)]
  #   if (any(is.na(gene2GO)))
  #     gene2GO <- gene2GO[!is.na(gene2GO)]
  #   gene2GO <- gene2GO[sapply(gene2GO, length) > 0]
  #   allGO <- unlist(gene2GO, use.names = FALSE)
  #   geneID <- rep(names(gene2GO), sapply(gene2GO, length))
  #   goodGO <- allGO %in% ls(ontoGO)
  #   return(split(geneID[goodGO], allGO[goodGO]))
  # }

  GOdata <- new("topGOdata", ontology = aspect, allGenes = compared_proteins,
                annot = annFUN.gene2GO, gene2GO = bg2GO)
  ## ------- Run Fisher test
  resultFisher <- topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher")

  ## --- Create table with enrichment result
  output <- topGO::GenTable(GOdata, classicFisher = resultFisher, topNodes = n)
  output <- as.data.frame(output)

  return(output)
}



## ---------------------------------------------- ##
#             Testing  go.enrich                   #
## ---------------------------------------------- ##
test_that(" go.enrich() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- go.enrich(target = c('Q14667', 'Q5JSZ5'),
                 background = background.go(c("Q13015", "Q14667", "P08575", "Q5JSZ5", "P13196")),
                 aspect = 'BP', n = 20)
  b <- go.enrich(target = c('Q14667', 'Q5JSZ5'), bg, aspect = 'BP', n = 20)

})
