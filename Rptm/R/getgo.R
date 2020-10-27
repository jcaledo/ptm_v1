get.go <- function(id, filter = TRUE, format = 'dataframe', silent = FALSE){

  if (!silent){
    print(paste("Getting GO terms for ", id, sep = ""))
  }

  ## ------------------------- Subfunction for complet list ---------------------- ##
  complet.list <- function(id){
    requestURL <- paste("https://www.ebi.ac.uk/QuickGO/services/annotation/",
                        "downloadSearch?includeFields=goName&selectedFields=symbol&geneProductId=",
                        id, sep = "")
    r <- httr::GET(requestURL, httr::accept("text/gpad"))
    httr::stop_for_status(r)
    content <- httr::content(r,as = "text")
    a <- strsplit(content, split = "\n")[[1]]
    lines <- a[10:length(a)]
    output <- as.data.frame(matrix(rep(NA, length(lines)*8), ncol = 8))
    names(output) <- c('gene_product', 'qualifier', 'GO_id', 'evidence', 'evidence_code',
                       'reference', 'assigned_by', 'date')
    for (i in seq_len(length(lines))){
      t <- strsplit(lines[i], split = "\t")[[1]]
      output$gene_product[i] <- t[2]
      output$qualifier[i] <- t[3]
      output$GO_id[i] <- t[4]
      output$evidence[i] <- t[6]
      output$evidence_code[i] <- strsplit(t[12], split = "=")[[1]][2]
      output$reference[i] <- t[5]
      output$assigned_by[i] <- t[10]
      output$date[i] <- t[9]
    }
    return(output)
  }


  ## ------------------------- Subfunction for filtered list --------------------- ##
  filtered.list <- function(id){
    baseURL <- 'https://www.uniprot.org/uniprot/?query='
    requestURL <- paste(baseURL, id, '&format=tab&columns=id%2Cgo', sep = "")
    resp <- .get.url(requestURL)
    cont <- httr::content(resp, 'text')
    a <- strsplit(cont, split = '\n')[[1]] # all the lines
    b <- a[which(grepl(id, a))]
    c <- strsplit(b, split = "\t")[[1]][2] # only term names and GO ids
    d <-  strsplit(c, split = ";")[[1]] # a single line by term

    output <- as.data.frame(matrix(rep(NA,length(d)*2), ncol = 2))
    names(output) <- c('term_name', 'GO_id')
    for (i in 1:length(d)){
      output$term_name[i] <- trimws( strsplit(d[i], split = '\\[')[[1]][1] )
      output$GO_id[i] <- gsub('\\]', '', strsplit(d[i], split = '\\[')[[1]][2])
    }
    ## ---- Removing spurious rows if necessary
    output <- output[which(substr(output$GO_id, 1, 2) == "GO"), ]

    if (sum(is.na(output$GO_id)) == nrow(output)){
      return(paste("Sorry, no GO terms found for the ", id, " entry", sep = ""))
    } else {
      output$obsolete <- output$definition_text <- output$aspect <- NA
      for (i in 1:nrow(output)){
        t <- strsplit(output$GO_id[i], split = ":")[[1]][2]
        url <- 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query=GO%3A'
        call <- paste(url, t, '&limit=1&page=1', sep = "")
        resp <- .get.url(call)
        cont <- httr::content(resp, 'text')
        cont <- jsonlite::fromJSON(cont, flatten = TRUE)$results

        if ("isObsolete" %in% names(cont)){
          output$obsolete[i] <- cont$isObsolete
        }
        if ("definition.text" %in% names(cont)){
          output$definition_text[i] <- cont$definition.text
        }
        if ("aspect" %in% names(cont)){
          output$aspect[i] <- cont$aspect
        }
      }
    }
    return(output)
  }

  ## ------- Building the output dataframe ----------------- ##
  if (filter){
    output <- filtered.list(id)
  } else {
    output <- complet.list(id)
  }

  if (format == 'string'){
    output <- paste(output$GO_id, collapse = ", ")
  }

  return(output)
}

