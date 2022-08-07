## ---------------- go.R --------------- ##
#                                         #
#     search.go                           #
#     term.go                             #
#     get.go                              #
#     bg.go                               #
#     hdfisher.go                         #
#     net.go                              #
#                                         #
## ------------------------------------- ##

## ---------------------------------------------------------------- ##
#           search.go <- function(query)                             #
## ---------------------------------------------------------------- ##
#' Search a Simple User Query
#' @description Searches a simple user query.
#' @usage search.go(query)
#' @param query character string defining the query.
#' @return Returns a dataframe containing the GO IDs found associated to the query, as well as other information related to these terms.
#' @author Juan Carlos Aledo
#' @seealso term.go(), get.go(), bg.go(), hdfisher.go(), gorilla(), net.go()
#' @references Rhee et al. (2008) Nature Reviews Genetics 9:509–515.
#' @examples \dontrun{search.go('oxidative stress')}
#' @importFrom jsonlite fromJSON
#' @export

search.go <- function(query){

  query <- tolower(query)
  query_ <- gsub(" ", '%20', query)
  baseURL <- "https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query="
  requestURL <- paste(baseURL, query_, "&limit=600", sep = "")

  r <- gracefully_fail(requestURL)

  if (is.null(r)){
    message("Sorry, no result could be retrieved")
    return(NULL)
  } else if (jsonlite::fromJSON(r)$numberOfHits == 0) {
    message("Sorry, no hits found")
    return(NULL)
  } else {
    output <- as.data.frame(jsonlite::fromJSON(r))[, 2:6]
    names(output) <- c("GO_id", "obsolete", "term_name", "definition_text", "aspect")
    attr(output, "query") <- query
    return(output)
  }
}


## ---------------------------------------------------------------- ##
#           term.go <- function(go, children = FALSE)                #
## ---------------------------------------------------------------- ##
#' Get Core Information About the GO Term
#' @description Gets core information about the GO term of interest.
#' @usage term.go(go, children = FALSE)
#' @param go GO id.
#' @param children logical, when true GO children terms are returned.
#' @details When the argument children is set to TRUE, the output of this function is a list with two elements: the first one is a dataframe with the core information, and the second one is a dataframe containing the children terms.
#' @return Returns a dataframe containing core information such as term name and definition, reference, aspect, and whether or not the term is obsolete. If children is set to TRUE, the function returns a list.
#' @author Juan Carlos Aledo
#' @seealso search.go(), get.go(), bg.go(), hdfisher.go(), gorilla(), net.go()
#' @references Rhee et al. (2008) Nature Reviews Genetics 9:509–515.
#' @examples \dontrun{term.go('GO:0034599')}
#' @importFrom jsonlite fromJSON
#' @export

term.go <- function(go, children = FALSE){

  baseURL <- "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A"
  id <- strsplit(go, split = ":")[[1]][2]
  requestURL <- paste(baseURL, id, sep = "")

  r <- gracefully_fail(requestURL)
  if (is.null(r)){
    message("Sorry, no result could be retrieved")
    return(NULL)
  }

  t <- as.list(jsonlite::fromJSON(r)[[2]])
  core <- as.data.frame(matrix(rep(NA, 7), ncol = 7))
  names(core) <- c("term_name","GO_id","aspect","definition_text","reference", "synonyms","obsolete")
  if (!is.null(t$name)) core$term_name[1] = t$name[[1]]
  if (!is.null(t$id)) core$GO_id[1] = t$id[[1]]
  if (!is.null(t$aspect)) core$aspect[1] = t$aspect[[1]]
  if (!is.null(t$definition$text)) core$definition_text[1] = t$definition[[1]][[1]]
  if (!is.null(t$definition$xrefs)) core$reference[1] = paste(unlist(t$definition$xrefs[[1]][1]),
                                                              unlist(t$definition$xrefs[[1]][2]), sep = ":")
  if (!is.null(t$isObsolete)) core$obsolete[1] = t$isObsolete[[1]]
  if (!is.null(t$children)){
    children_terms <- as.data.frame(t$children)
  } else {
    children_terms <- "Sorry, no children couldn't be found"
  }

  # r <- httr::GET(requestURL, httr::accept("application/json"))
  # httr::stop_for_status(r)
  # json <- jsonlite::toJSON(httr::content(r))
  # t <- as.list(jsonlite::fromJSON(json)[[2]])

  if (children){
    output <- list(core, children_terms)
  } else {
    output <- core
  }
  return(output)
}

## ---------------------------------------------------------------- ##
#   get.go <- function(id, format = 'dataframe', silent = FALSE)     #
## ---------------------------------------------------------------- ##
#' Get Gene Ontology Annotation
#' @description Gets the gene ontology annotations for a given protein.
#' @usage get.go(id, format = 'dataframe', silent = FALSE)
#' @param id the UniProt identifier of the protein of interest.
#' @param format string indicating the output's format. It should be either 'dataframe' or 'string'. The 'string' format may be convenient when subsequent GO terms enrichment analysis is intended.
#' @param silent logical, if FALSE print details of the reading process.
#' @return Returns a dataframe (by default) with GO IDs linked to the protein of interest, as well as additional information related to these GO ids. A string with the GO ids can be obtained as output if indicated by means of the argument 'format'.
#' @author Juan Carlos Aledo
#' @seealso search.go, term.go(), bg.go(), hdfisher.go(), gorilla(), net.go()
#' @references Rhee et al. (2008) Nature Reviews Genetics 9:509–515.
#' @examples \dontrun{get.go('P01009')}
#' @importFrom httr GET
#' @importFrom httr accept
#' @importFrom httr content
#' @importFrom jsonlite fromJSON
#' @export

get.go <- function(id, format = 'dataframe', silent = FALSE){

  if (!silent){
    print(paste("Getting GO terms for ", id, sep = ""))
  }

  ## ------------------------- Subfunction for complet list ---------------------- ##
  complet.list <- function(id){
    requestURL <- paste("https://www.ebi.ac.uk/QuickGO/services/annotation/",
                        "downloadSearch?includeFields=goName&selectedFields=symbol&geneProductId=",
                        id, sep = "")

    r <- tryCatch(
      {
        httr::GET(requestURL, httr::accept("text/gpad"))
      },
      error = function(cond){
        return(NULL)
      }
    )
    if (is.null(r)){
      message("Sorry, no result could be retrieved")
      return(NULL)
    } else if (r$status_code > 300){
      message("Sorry, no result could be retrieved")
      return(NULL)
    }
    content <- httr::content(r, as = "text")

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


  # ## ------------------------- Subfunction for filtered list --------------------- ##
  # filtered.list <- function(id){
  #   baseURL <- 'https://www.uniprot.org/uniprot/?query='
  #   requestURL <- paste(baseURL, id, '&format=tab&columns=id%2Cgo', sep = "")
  #   cont <- gracefully_fail(requestURL)
  #   if (is.null(cont)){
  #     message("Sorry, no result could be retrieved")
  #     return(NULL)
  #   } else if (cont == ""){
  #     message("Sorry, no result could be retrieved")
  #     return(NULL)
  #   }
  #
  #   a <- strsplit(cont, split = '\n')[[1]] # all the lines
  #   b <- a[which(grepl(id, a))]
  #   if (length(b) == 0){
  #     return(NULL)
  #   }
  #   c <- strsplit(b, split = "\t")[[1]][2] # only term names and GO ids
  #   d <-  strsplit(c, split = ";")[[1]] # a single line by term
  #
  #   output <- as.data.frame(matrix(rep(NA,length(d)*2), ncol = 2))
  #   names(output) <- c('term_name', 'GO_id')
  #
  #   for (i in 1:length(d)){
  #     output$term_name[i] <- trimws( strsplit(d[i], split = '\\[')[[1]][1] )
  #     output$GO_id[i] <- gsub('\\]', '', strsplit(d[i], split = '\\[')[[1]][2])
  #   }
  #   ## ---- Removing spurious rows if necessary
  #   output <- output[which(substr(output$GO_id, 1, 2) == "GO"), ]
  #
  #   if (sum(is.na(output$GO_id)) == nrow(output)){
  #     message(paste("Sorry, no GO terms found for the ", id, " entry", sep = ""))
  #     return(NULL)
  #   } else {
  #     output$obsolete <- output$definition_text <- output$aspect <- NA
  #     for (i in 1:nrow(output)){
  #       t <- strsplit(output$GO_id[i], split = ":")[[1]][2]
  #       url <- 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query=GO%3A'
  #       call <- paste(url, t, '&limit=1&page=1', sep = "")
  #       resp <- .get.url(call)
  #       cont <- httr::content(resp, 'text')
  #       cont <- jsonlite::fromJSON(cont, flatten = TRUE)$results
  #
  #       if ("isObsolete" %in% names(cont)){
  #         output$obsolete[i] <- cont$isObsolete
  #       }
  #       if ("definition.text" %in% names(cont)){
  #         output$definition_text[i] <- cont$definition.text
  #       }
  #       if ("aspect" %in% names(cont)){
  #         output$aspect[i] <- cont$aspect
  #       }
  #     }
  #   }
  #   return(output)
  # }

  ## ------- Building the output dataframe ----------------- ##
  # if (filter){
  #   output <- filtered.list(id)
  # } else {
  #   output <- complet.list(id)
  # }

  output <- complet.list(id)
  if (format == 'string' & !is.atomic(output)){
    output <- paste(output$GO_id, collapse = ", ")
  }

  return(output)
}

## ---------------------------------------------------------------- ##
#                   bg.go <- function(ids)                           #
## ---------------------------------------------------------------- ##
#' Search GO Terms for Background Set
#' @description Searches the GO terms of the protein contained in a given set.
#' @usage bg.go(ids)
#' @param ids either a vector containing the UniProt IDs of the background set or the path to the txt file containing the list of IDs acting as background.
#' @return Returns a dataframe with two columns (Uniprot ID, GO terms) and as many rows as different proteins there are in the input set.
#' @author Juan Carlos Aledo
#' @references Rhee et al. (2008) Nature Reviews Genetics 9:509–515.
#' @seealso search.go(), term.go(), get.go(), go.enrich(), gorilla(), net.go()
#' @examples \dontrun{bg.go(c('P01009', 'P01374', 'Q86UP4'))}
#' @importFrom utils read.csv
#' @export

bg.go <- function(ids){
  ## ----- The background set
  if (is.character(ids) & length(ids) == 1){ # input as path to the txt
    if (gregexpr('txt', ids)[[1]] != -1){
      bg <- utils::read.csv(ids, header = FALSE)
      bg <- trimws(as.character(bg$V1))
    } else {
      stop("A proper path to a txt file should be provided for the background set")
    }
  } else if (is.character(ids) & length(ids) > 1){ # input as vector
    bg <- as.character(ids)
  } else if (is.data.frame(ids) & nrow(ids) > 1){ # input as dataframe
    bg <- trimws(as.character(ids))
  } else {
    stop("A proper background set must be provided")
  }

  ## ----- Getting GO ids for the backgraund set
  bg <- data.frame(up_id = bg, GO_id = rep(NA, length(bg)))

  for (i in 1:nrow(bg)){
    bg$GO_id[i] <- tryCatch(
      {
        get.go(trimws(bg$up_id[i]), format = 'string')
      },
      error = function(cond){
        return(NA)
      }
    )

  }
  return(bg)
}


## ---------------------------------------------------------------- ##
#   hdfisher.go <- function(target, background, query,               #
#                                  analysis = 'enrichment')          #
## ---------------------------------------------------------------- ##
#' Hypothesis-Driven Fisher Test
#' @description Carries out an enrichment Fisher's test using a hypothesis driven approach.
#' @usage hdfisher.go(target, background, query, analysis = 'enrichment')
#' @param target either a vector containing the UniProt IDs of the target set or the path to the txt file containing the list of IDs.
#' @param background  a dataframe with two columns (Uniprot ID and GO terms) and as many rows as different proteins there are in the background set.
#' @param query character string defining the query.
#' @param analysis a character string indicating whether the desired analysis is the enrichment ('enrichment') or depletion ('depletion').
#' @return Returns a list that contains the contingency table and the p-Value.
#' @author Juan Carlos Aledo
#' @references Rhee et al. (2008) Nature Reviews Genetics 9:509–515.
#' @seealso search.go(), term.go(), get.go(), bg.go(), go.enrich(), gorilla(), net.go()
#' @examples \dontrun{hdfisher.go(c('Q14667', 'Q5JSZ5'), bg.go(c('Q14667', 'Q5JSZ5', 'P13196')), 'ion')}
#' @importFrom utils read.csv
#' @importFrom stats fisher.test
#' @export

hdfisher.go <- function(target, background, query, analysis = 'enrichment'){
  ## ----- The target sample to be analyzed
  if (is.character(target) & length(target) == 1){ # input as path to the txt
    if (gregexpr('txt', target)[[1]] != -1){
      sample <- utils::read.csv(target, header = FALSE)
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
  target_c <- setdiff(bg, sample) # target complement

  rquery <- tryCatch(
    {
      search.go(query)
    },
    error = function(cond){
      return(NULL)
    }
  )
  if (is.null(rquery)){
    message("Sorry, no hits were found for the current query")
    return(NULL)
  } else {
    rquery <- unlist(rquery$GO_id)

    # -- Contingency Table

    df_target <- background[which(background$up_id %in% sample), ]
    df_target_c <- background[which(background$up_id %in% target_c), ]
    df_target$query <- df_target_c$query <- NA

    for (i in 1:nrow(df_target)){
      t <- strsplit(df_target$GO_id[i], split = ',')[[1]]
      t <- trimws(t)
      if (length(intersect(t, rquery)) > 0){
        df_target$query[i] <- TRUE
      } else {
        df_target$query[i] <- FALSE
      }
    }
    # a: number of proteins from the target set with terms present into the query
    a <- sum(df_target$query)
    # c: number of protein from the target set which terms are absent from the query
    c <- nrow(df_target) - a


    for (i in 1:nrow(df_target_c)){
      t <- strsplit(df_target_c$GO_id[i], split = ',')[[1]]
      t <- trimws(t)
      if (length(intersect(t, rquery)) > 0){
        df_target_c$query[i] <- TRUE
      } else {
        df_target_c$query[i] <- FALSE
      }
    }
    # b: number of proteins from the target complement set with terms present into the query
    b <- sum(df_target_c$query)
    # d: number of protein from the target complement set which terms are absent from the query
    d <- nrow(df_target_c) - b

    # -- Fisher's test
    if (analysis == 'enrichment'){
      alternative <- 'greater'
    } else if (analysis == 'depletion'){
      alternative <- 'less'
    } else {
      stop("a suitable analysis, either 'enrichment' or 'deplation' must be indicated")
    }
    ct <- matrix(c(a,b,c,d), nrow = 2, byrow = TRUE)
    colnames(ct) <- c('target', 'target-complement')
    rownames(ct) <- c('query', 'non-query')

    ft <- stats::fisher.test(ct, alternative = alternative)
    output <- list(contingency_table = ct, pv = ft$p.value)
  }

  attr(output, 'query') <- query
  attr(output, 'analysis') <- analysis
  return(output)
}


## ---------------------------------------------------------------- ##
#     net.go <- function(data, threshold = 0.2, silent = FALSE)      #                             #
## ---------------------------------------------------------------- ##
#' Gene Ontology Network
#' @description Explores the relationship among proteins from a given set.
#' @usage net.go(data, threshold = 0.2, silent = FALSE)
#' @param data either a vector containing the UniProt IDs (vertices) or the path to the txt or rda file containing them.
#' @param threshold threshold value of the Jaccard index above which two proteins are considered to be linked.
#' @param silent logical, if FALSE print details of the running process.
#' @details This function first searches the GO terms for each vertex and then computes the Jaccard index for each protein pair, based on their GO terms. Afterwards, an adjacency matrix is computed, where two proteins are linked if their Jaccard index is greater than the selected threshold.
#' @return Returns a list containing (i) the dataframe corresponding to the computed Jaccard matrix, (ii) the adjacency matrix, (iii) a vector containing the vertices, and (iv) a matrix describing the edges of the network.
#' @author Pablo Aledo & Juan Carlos Aledo
#' @seealso search.go(), term.go(), get.go(), bg.go(), gorilla()
#' @references Aledo & Aledo (2020) Antioxidants 9(10), 987.
#' @references Rhee et al. (2008) Nature Reviews Genetics 9:509–515.
#' @examples \dontrun{net.go(path2data = "./GOvivo.txt")}
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph get.edgelist
#' @export

net.go <- function(data, threshold = 0.2, silent = FALSE){

  ## ------- Assessing whether data are in txt or rda format -------- ##
  format <- format_ <- ""
  if (is.vector(data) & length(data) > 1){ # --- Vertices are directly provided as a vector
    vertices <- data
    # id <- vertices[,1]
  } else if (is.character(data)){
    format <- strsplit(data, split = "\\.")[[1]]
    format <- format[length(format)]
    format_ <- tolower(format)
    if (! format_ %in% c('txt', 'rda')){
      stop("Please, make sure that data are in either 'txt' or 'rda' format")
    }
  }
  if (format_ == 'txt'){ # ----- Vertices are provided as txt file
    con <- file(data, 'r')
    vertices <- readLines(con)
    close(con)
    # vertices.df <- data.frame(vertex = id)
  } else if (format_ == 'rda'){ # ----- Vertices are provided as rda file
    load(data)
    l <- ls()[! ls() %in% c('format', 'format_', 'silent', 'threshold')]
    if (length(l) > 1){
      vertices <- get(l[which(l != "data")])
    } else {
      vertices <- get(l)
    }
  }
  id <- vertices
  ## ----------------- Computing f(id) = GO_subset ------------------ ##
  fid <- lapply(id, function(x) get.go(id = x, format = "string"))
  valid_id <- which(fid != "NULL") # removing proteing without GO annotations
  fid <- fid[valid_id]
  id <- id[valid_id]
  ## ---------- Computing Jaccard index in the id x id set ---------- ##
  jaccard <- matrix(rep(NA, length(id)^2), ncol =length(id))
  colnames(jaccard) <- rownames(jaccard) <- id

  for (i in 1:(length(id) - 1)){
    if (!silent){
      print(paste(i, "  .....  ", id[i], sep = ""))
    }
    for (j in (i+1):length(id)){
      A <- unique(lapply(strsplit(fid[[i]], split = ","), function(x) trimws(x))[[1]])
      B <- unique(lapply(strsplit(fid[[j]], split = ","), function(x) trimws(x))[[1]])
      AuB <- length(union(A,B))
      AB <- length(intersect(A,B))
      jaccard[i,j] <- round(AB/AuB, 3)
    }
  }

  ## -------------------- From Jaccard to Adjacency ------------------ ##
  A <- as.matrix(jaccard)
  A[A >= threshold] <- 1
  A[A < threshold] <- 0
  diag(A) <- 0
  A[is.na(A)] <- 0
  A <- A + t(A)

  ## ------------------------- Network ------------------------------- ##
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
  edges.df <- igraph::get.edgelist(g, names=TRUE)

  ## -------------------------- Output ------------------------------- ##
  output <- list(jaccard, A, vertices, trimws(edges.df))
  attr(output, 'Jaccard threshold') <- threshold
  return(output)
}


