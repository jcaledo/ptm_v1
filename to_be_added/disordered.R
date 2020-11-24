## ------- disordered.R ------------ ##
#                                     #
#       disor.scan                    #
#       disor.met                     #
#                                     #
## --------------------------------- ##

## --------------------------------------------------------------- ##
#                            disor.scan(id)                           #
## --------------------------------------------------------------- ##
#' Search for Intrinsically Disordered Region 
#' @description Scans the protein in search of disordered regions
#' @usage disor.scan(id)
#' @param id the UniProt identifier of the protein of interest.
#' @details This function carries out a query into DisProt, which is a database of intrinsically disordered proteins. In this database disordered regions are manually curated from literature. 
#' @return  If the protein being scanned contains disordered region, the function returns a dataframe with as many rows as regions found. For each row, the positions at which the region starts and ends are indicated, as well as the type of region and the sequence of the region. If the protein being scanned is not present in the dabase, the function returns 'status code: 204'. If the protein is found in the database but the experimental evidences are ambiguous, the function will return an empty list.
#' @author Juan Carlos Aledo
#' @examples disor.scan("P04637")
#' @references Hatos et al. Nucleic Acids Res. 2019, 10.1093/nar/gkz975. (PMID: 31713636)
#' @seealso disor.target()
#' @export
 

# https://www.disprot.org/help




disor.scan <- function(id){
  
  call <- paste('https://www.disprot.org/api/', id,  sep = "")
  response <- httr::GET(call)
  
  if (response$status_code == 200){
    output <- httr::content(response, 'text')
    output <- jsonlite::fromJSON(output, flatten = TRUE)
    output <- output$disprot_consensus$structural_state
    if (length(output) > 0){
      prot <- ptm::get.seq(id)
      output$seq <-NA
      for (i in 1:nrow(output)){
        output$seq[i] <- substr(prot, start = output$start[i], stop = output$end[i])
      }
    }
    return(output)
  } else {
    return(paste('status code: ', response$status_code, sep = ""))
  }
}

## --------------------------------------------------------------- ##
#                     disor.target(id, target)                      #
## --------------------------------------------------------------- ##
#' Is The Target Within Intrinsically Disordered Regions 
#' @description Checks whether the target is found within a disordered region
#' @usage disor.target(id)
#' @param id the UniProt identifier of the protein of interest.
#' @param target either a position or a character string indicating the sequence of an oligopeptide of interest.
#' @details This function carries out a query into DisProt, which is a database of intrinsically disordered proteins. In this database disordered regions are manually curated from literature. 
#' @return  
#' @author Juan Carlos Aledo
#' @examples disor.target("P04637", )
#' @references Hatos et al. Nucleic Acids Res. 2019, 10.1093/nar/gkz975. (PMID: 31713636)
#' @seealso disor.scan()
#' @export

disor.target <- function(id, target){
  output <- 0
  t <- disor.scan(id)
  if (t[[1]][1] == "status code: 204"){
   return(FALSE)
  } else if (is.data.frame(t) & is.numeric(target)){
      for (i in 1:nrow(t)){
        if (t$start[i] < target & target < t$end[i]){
          output <- output + 1
        } 
      }
      if (output > 0){
        output <- TRUE
      } else {
        output <- FALSE
      }
      attr(output, 'target') <- target
  } else if (is.data.frame(t) & !is.numeric(target)){
      target <- toupper(target)
      for (i in 1:nrow(t)){
        if (gregexpr(target, t$seq[i])[[1]][1] > 0){
          output <- TRUE
          attr(output, 'target') <- gregexpr(target, t$seq[i])[[1]][1]
        }
      }  
  }
  return(output)
}
   

t <- disor.target("P49913") # en db con ambiguedad
t <- disor.target("P01009") # no en db
t <- disor.target("P04637", "PPVAPAPAAP") # ctr positivo

t <- disor.scan("P49913") # en db con ambiguedad
t <- disor.scan("P01009") # no en db
t <- disor.scan("P04637") # ctr positivo
rm(list = ls())
id <- "P04637"
target <- 90
