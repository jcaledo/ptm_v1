## ------------- aa.R ------------- ##
#                                    #
#     abundance                      #
#                                    #
## -------------------------------- ##

## --------------------------------------------------------------- ##
#                         abundance(id, ...)                        #
## --------------------------------------------------------------- ##
#' Residue Found at the Requested Position
#' @description Returns the residue found at the requested position
#' @usage abundance(id, ...)
#' @param id the UniProt identifier of the protein of interest.
#' @details For human proteins, in addition to the abundance in the whole organism (by default), the abundance found in Jurkat or HeLa cells can be requested.
#' @return A numeric value for the abundance, expressed a parts per million (ppm), of the requested protein.
#' @author Juan Carlos Aledo
#' @examples abundance(id = 'A0AVT1')
#' abundance(id = 'A0AVT1', 'jurkat')
#' abundance(id = 'A0AVT1', 'hela')
#' @references Wang et al. Proteomics 2015, 10.1002/pmic.201400441. (PMID: 25656970)
#' @seealso
#' @export

abundance <- function(id, ...){

  tissue <- list(...)
  sp <- species.mapping(id)

  if (sp == 'Arabidopsis thaliana'){
    paxdb <- ptm:::pax.ath
  } else if (sp == 'Bos taurus' ){
    paxdb <- ptm:::pax.bta
  } else if (sp == 'Dictyostelium discoideum'){
    paxdb <- ptm:::pax.ddi
  } else if (sp == 'Drosophila melanogaster'){
    paxdb <- ptm:::pax.dme
  } else if (sp == 'Equus caballus'){
    paxdb <- ptm:::pax.ecb
  } else if (sp == 'Escherichia coli'){
    paxdb <- ptm:::pax.eco
  } else if (sp == 'Gallus gallus'){
    paxdb <- ptm:::pax.gga
  } else if (sp == 'Homo sapiens'){
    if (length(tissue) == 0){
      paxdb <- ptm:::pax.hsa
    } else if (tissue[[1]] == 'jurkat'){
      paxdb <- ptm:::pax.jurkat
    } else if (tissue[[1]] == 'hela') {
      paxdb <- ptm:::pax.hela
    } else {
      stop('A proper tissue should be provided!')
    }
  } else if (sp == 'Mus musculus'){
    paxdb <- ptm:::pax.mmu
  } else if (sp == 'Mycobacterium tuberculosis'){
    paxdb <- ptm:::pax.mtu
  } else if (sp == 'Rattus norvegicus'){
    paxdb <- ptm:::pax.rno
  } else if (sp == 'Saccharomyces cerevisiae'){
    paxdb <- ptm:::pax.sce
  } else if (sp == 'Sus scrofa'){
    paxdb <- ptm:::pax.scr
  } else {
    warning <- paste('Abundance data for proteins of the species ',
                     sp, " couldn't be found", sep = "")
    return(warning)
  }

  paxdb$abundance <- as.numeric(as.character(paxdb$abundance))
  output <- paxdb$abundance[which(paxdb$up_id == id)]
  attr(output, 'units') <- 'ppm'
  attr(output, 'species') <- sp
  attr(output, 'string') <- id
  return(output)
}
