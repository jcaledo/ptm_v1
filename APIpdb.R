# https://data.rcsb.org/#data-api
# http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items
# https://www.rcsb.org/pages/download/http

# former endpoint:
#      "https://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList="
# current endpoints:
#  https://www.rcsb.org/fasta/entry/{entry}/download
#  https://data.rcsb.org/rest/v1/core/entry/{entry_id}
#  https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_id}


pdb.seq <- function(pdb){

  ## ----------------- Check PDB argument ------------ #
  if (nchar(pdb) != 4){
    stop( "Please, provide a proper PDB ID")
  } else {
    call <- paste('https://www.rcsb.org/fasta/entry/',
                  pdb, '/download', sep = "")
  }

  ## -------- Client <-> Server Communication -------- #
  if (!is.null(call)){
    resp <- try(httr::GET(call), silent = FALSE)
    if (inherits(resp, "try-error")) {
      text <- "LOST CONNECTION"
    } else {
      text <- httr::content(resp, "text", encoding = "utf-8")
    }
    retry = 0
    while ((grepl("^LOST CONNECTION", text) || httr::http_error(resp) ||
            grepl("^ERROR", text) || grepl("Nothing has been found",
                                           text)) && retry < 3) {
      retry = retry + 1
      Sys.sleep(5)
      resp <- try(httr::GET(url), silent = TRUE)
      if (inherits(resp, "try-error")) {
        text <- "LOST CONNECTION"
      }
      else {
        text <- httr::content(resp, "text", encoding = "utf-8")
      }
    }
  }
  ## ----------- Parsing the response ------------- #
  t <- strsplit(text, split = ">")[[1]][-1]
  seq <- data.frame(entry = rep(NA, length(t)),
                    entity = rep(NA, length(t)),
                    chain = rep(NA, length(t)),
                    name = rep(NA, length(t)),
                    species = rep(NA, length(t)),
                    sequence = rep(NA, length(t)))

  for (i in 1:nrow(seq)){
    tt <- strsplit(t[i], split = "\n")[[1]]
    z <- strsplit(tt[1], "\\|")[[1]]
    seq$entry[i] <- strsplit(z[1], split = "_")[[1]][1]
    seq$entity[i] <- strsplit(z[1], split = "_")[[1]][2]
    chains <- gsub('Chains', "", z[2])
    seq$chain[i] <- trimws(chains)
    seq$name[i] <- z[3]
    seq$species[i] <- z[4]
    seq$sequence[i] <- tt[2]
  }
  ## ------------------ Output ------------------- #
  return(seq)
}



