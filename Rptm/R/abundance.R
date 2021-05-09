## ------ abundance.R ------------- ##
#                                    #
#     abundance                      #
#                                    #
## -------------------------------- ##

## --------------------------------------------------------------- ##
#                         abundance(id, ...)                        #
## --------------------------------------------------------------- ##
#' Protein Abundance Data
#' @description Returns data regarding the abundance of a given protein.
#' @usage abundance(id, ...)
#' @param id the UniProt identifier of the protein of interest.
#' @param ... either 'jarkat' or 'hela' if required.
#' @details For human proteins, in addition to the abundance in the whole organism (by default), the abundance found in Jurkat or HeLa cells can be requested. The data are obtained from the PaxDb.
#' @return A numeric value for the abundance, expressed a parts per million (ppm), of the requested protein.
#' @author Juan Carlos Aledo
#' @examples abundance(id = 'A0AVT1')
#' @examples \dontrun{abundance(id = 'A0AVT1', 'jurkat')}
#' @examples \dontrun{abundance(id = 'A0AVT1', 'hela')}
#' @references Wang et al. Proteomics 2015, 10.1002/pmic.201400441. (PMID: 25656970)
#' @export

abundance <- function(id, ...){

  tissue <- list(...)

  sp <- species.mapping(id)
  if (is.null(sp)){
    message("Sorry, data abundance for the requested protein could not be found")
    return(NULL)
  }

  if (sp == 'Arabidopsis thaliana'){
    paxdb <- pax.ath
  } else if (sp == 'Bos taurus' ){
    paxdb <- pax.bta
  } else if (sp == 'Dictyostelium discoideum'){
    paxdb <- pax.ddi
  } else if (sp == 'Drosophila melanogaster'){
    paxdb <- pax.dme
  } else if (sp == 'Equus caballus'){
    paxdb <- pax.ecb
  } else if (sp == 'Escherichia coli'){
    paxdb <- pax.eco
  } else if (sp == 'Gallus gallus'){
    paxdb <- pax.gga
  } else if (sp == 'Homo sapiens'){
    if (length(tissue) == 0){
      paxdb <- pax.hsa
    } else if (tissue[[1]] == 'jurkat'){
      paxdb <- pax.jurkat
    } else if (tissue[[1]] == 'hela') {
      paxdb <- pax.hela
    } else {
      stop('A proper tissue should be provided!')
    }
  } else if (sp == 'Mus musculus'){
    paxdb <- pax.mmu
  } else if (sp == 'Mycobacterium tuberculosis'){
    paxdb <- pax.mtu
  } else if (sp == 'Rattus norvegicus'){
    paxdb <- pax.rno
  } else if (sp == 'Saccharomyces cerevisiae'){
    paxdb <- pax.sce
  } else if (sp == 'Sus scrofa'){
    paxdb <- pax.ssc
  } else {
    message(paste('Abundance data for proteins of the species ',
                     sp, " couldn't be found", sep = ""))
    return(NULL)
  }

  paxdb$abundance <- as.numeric(as.character(paxdb$abundance))
  output <- paxdb$abundance[which(paxdb$up_id == id)]
  attr(output, 'units') <- 'ppm'
  attr(output, 'species') <- sp
  attr(output, 'string') <- id
  return(output)
}
