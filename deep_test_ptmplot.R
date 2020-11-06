rm(list = ls())
library(ptm)
setwd("/Users/JCA/Dropbox/Investigacion/Proteomes")
load("./human_proteome.Rda")
h <- h[1:100, ]
h$test <- h$time <- NA
for (i in 44:100){
  print(i)
  a <- Sys.time()
  h$test[i] <- ptm.plot(h$prot_id[i],
                        property = 'acc',
                        ptm = 'reg')
  b <- Sys.time()
  h$time[i] <- b-a
  closeAllConnections()
}


## ----
# testado con i de 1 a 209
# Error server didn't responde para i = 14 (5xtc), 90 (5a2q), 125 (3pja),
# 174 (3p8c) y 190 (4k1a).
# Error en bio3d::read.pdb() para i = 146 (5t2c) y 209
# up_id <- variable-i
# pdb <- ""
# property <- 'eiip'
# ptm <- 'meto'
# dssp <- 'compute'
# window <- 1
# sdata <- TRUE

## ----
# testado con i de 1 a 67
# Error para i = 14 (5xtc)
# Error list.home(target) status code 400 en i = 27, 35, 67 (id.mapping: "Sorry, no KEGG...")
# Error bio3d::seqaln() en 50 (Bus error 10)
# up_id <- variable-i
# pdb <- ""
# property <- 'entropy7.aa'
# ptm <- 'ac'
# dssp <- 'compute'
# window <- 1
# sdata <- TRUE


## ----
# testado con i de 1 a 49
# Solo se queja cuando no se encuentra pdb {11,15,19,21,23,24,28,31,32,34,35,36,43,49}
# up_id <- variable-i
# pdb <- ""
# property <- 'acc'
# ptm <- 'reg'
# dssp <- 'compute'
# window <- 1
# sdata <- TRUE

rm(list = ls())
up_id <- "Q6A555"
kegg <- id.mapping("Q6A555", 'uniprot', 'kegg')
# custom.aln(target = kegg, species = "vertebrates")
a <- shannon(target = kegg,
             species = "vertebrates",
             base = 21, alphabet = 21)


