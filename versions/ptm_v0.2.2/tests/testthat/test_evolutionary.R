library(ptm)
context("Evo traces")

## ---------------------------------------------- ##
#             Testing msa                          #
## ---------------------------------------------- ##
test_that("msa() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- msa(sequences = sapply(c("P19446", "P40925", "P40926"), ptm::get.seq),
           ids = c("wmelon", "cyt", "mit"))

  expect_is(a, "list")
  expect_equal(length(a), 4)
  expect_is(a[[1]], 'character')
  expect_equal(length(a[[1]]), 3)
  expect_equal(names(a[[1]]), c("P19446", "P40925", "P40926"))
  expect_equivalent(nchar(a[[1]][1]), 356)
  expect_is(a[[2]], 'character')
  expect_equal(length(a[[2]]), 3)
  expect_equal(a[[2]], c("wmelon", "cyt", "mit" ))
  expect_equal(length(a[[3]]), 3)
  expect_equivalent(nchar(a[[3]][1]), 383)
  expect_is(a[[4]], 'matrix')
  expect_equal(dim(a[[4]]), c(3, 383))
})

## ---------------------------------------------- ##
#             Testing custom.aln                   #
## ---------------------------------------------- ##
test_that("custom.aln() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- custom.aln('hsa:4069', species = c('pps', 'pon', 'mcc', 'ssc'))
  b <- custom.aln('cge:100773737', species = 'vertebrates', molecule = 'dna' )
  c <- custom.aln('eco:b2114', species = 'plants')
  d <- custom.aln('eco:b2114', species = 'one-hundred')
  e <- custom.aln('eco:b2114', species = 'two-hundred')

  expect_is(a, 'fasta')
  expect_equal(length(a), 3)
  expect_equal(a[[1]], c("ref", "pps", "pon", "mcc", "ssc"))

  expect_is(b, 'fasta')
  expect_equal(length(b), 3)
  expect_equal(b[[1]], c("ref", "hsa", "ptr", "rno", "ggo", "gga", "bta", "xtr", "dre"))

  expect_is(c, 'fasta')
  expect_equal(length(c), 3)
  expect_equal(c[[1]], c("ref", "cre", "boe", "aly", "gmx", "sly", "osa", "ara"))

  expect_is(d, 'fasta')
  expect_equal(length(d), 3)
  expect_equal(dim(d[[2]]), c(101, 680))

  expect_is(e, 'fasta')
  expect_equal(length(e), 3)
  expect_equal(dim(e[[2]]), c(201, 727))

})

## ---------------------------------------------- ##
#             Testing list.hom                     #
## ---------------------------------------------- ##
test_that("list.hom() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- list.hom(target = 'hsa:4069', homology = 'o')
  b <- list.hom(target = 'hsa:4069', homology = 'p') # "Lysozyme"
  c <- list.hom(target = 'hsa:2744', homology = 'o')
  d <- list.hom(target = 'hsa:2744', homology = 'p') # Glnase
  e <- list.hom(target = 'hsa:2597', homology = 'p') # GAPDH

  expect_is(a, 'data.frame')
  expect_gte(nrow(a), 307)
  expect_equal(ncol(a), 7)
  expect_true(a$name[1] == 'lysozyme')

  expect_is(b, 'data.frame')
  expect_gte(nrow(b), 8)
  expect_equal(ncol(b), 7)
  expect_true('lysozyme like 2' %in% b$name)

  expect_is(c, 'data.frame')
  expect_gte(nrow(c), 4174)
  expect_equal(ncol(c), 7)
  expect_true(c$name[1] == 'glutaminase kidney isoform, mitochondrial')

  expect_is(d, 'data.frame')
  expect_gte(nrow(d), 300)
  expect_equal(ncol(d), 7)
  expect_true('glutaminase 2' %in% d$name)

  expect_is(e, 'data.frame')
  expect_gte(nrow(e), 1)
  expect_equal(ncol(d), 7)
})


## ---------------------------------------------- ##
#             Testing parse.hssp                   #
## ---------------------------------------------- ##
test_that("parse.hssp() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- parse.hssp(file = "./pdb/1u8f.hssp", keepfiles = TRUE)

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 333)
  expect_equal(ncol(a), 28)

  expect_true(file.exists("1u8f_insertions.Rda"))
  expect_true(file.exists("1u8f_aln.Rda"))
  expect_true(file.exists("1u8f_profile.Rda"))
  expect_true(file.exists("1u8f_seq_list.Rda"))
})
file.remove("1u8f_insertions.Rda", "1u8f_aln.Rda", "1u8f_profile.Rda", "1u8f_seq_list.Rda")



## ---------------------------------------------- ##
#             Testing get.hssp                     #
## ---------------------------------------------- ##
test_that("get.hssp() works properly", {

  skip_on_cran()
  skip_on_travis()

  b <- get.hssp(pdb = "5xge")
  if (class(b) != "character"){ # skip if no local HSSP repo is present

      expect_is(b, 'data.frame')
      expect_equal(nrow(b), 549)
      expect_equal(ncol(b), 28)

      expect_true(file.exists("5xge_insertions.Rda"))
      expect_true(file.exists("5xge_aln.Rda"))
      expect_true(file.exists("5xge_profile.Rda"))
      expect_true(file.exists("5xge_seq_list.Rda"))

      file.remove("5xge_insertions.Rda", "5xge_aln.Rda", "5xge_profile.Rda", "5xge_seq_list.Rda")
  }
})


## ---------------------------------------------- ##
#             Testing shannon                      #
## ---------------------------------------------- ##
test_that("shannon() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- shannon('hsa:4069', species = 'vertebrates', base = 2, alphabet = 21)
  b <- shannon('hsa:4069', species = 'vertebrates', base = 2, alphabet = 4)
  c <- shannon('hsa:4069', species = 'vertebrates', base = 21, alphabet = 21)
  d <- shannon('hsa:4069', species = 'vertebrates', base = 4, alphabet = 4)


  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 148)
  expect_equal(ncol(a), 5)

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 148)
  expect_equal(ncol(b), 5)

  expect_is(c, 'data.frame')
  expect_equal(nrow(c), 148)
  expect_equal(ncol(c), 5)

  expect_is(d, 'data.frame')
  expect_equal(nrow(d), 148)
  expect_equal(ncol(d), 5)

  expect_gte(a$Haa[100], b$Haa[100])
})

## ---------------------------------------------- ##
#             Testing site.type                   #
## ---------------------------------------------- ##
test_that("site.type() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- site.type(target = 'hsa:4069', species = 'vertebrates') # lysozime

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 148)
  expect_equal(ncol(a), 5)
  expect_gte(a$H21[50], a$H4[50])

})
