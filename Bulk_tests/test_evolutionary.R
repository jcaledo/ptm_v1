library(ptm)
context("Evo traces")

## ---------------------------------------------- ##
#             Testing msa                          #
## ---------------------------------------------- ##
test_that("msa() works properly", {

  skip_on_cran()
  skip_on_travis()

  sequences <- c("MQPIPDVNQRIARISAHLHPPKSQMEESSALRRANCRAKGGAPGFKV",
                 "MSEPIRVLVTGAAGQIAYSLLYSIGNGSVFGKDQPIILVLLDITPMM",
                 "MSEPIRVLVTGAAGQIAYSLLYSIGNGSVFGKDQPIILVLLDITPMM")
  names(sequences) <- c("P19446", "P40925", "P40926")

  a <- msa(sequences = sequences, ids = c("wmelon", "cyt", "mit"))

  b <- msa(sequences = sequences,
           ids = c("wmelon", "cyt", "mit"), inhouse = TRUE)

  expect_is(a, "list")
  expect_equal(length(a), 4)
  expect_is(a[[1]], 'character')
  expect_equal(length(a[[1]]), 3)
  expect_equal(names(a[[1]]), c("P19446", "P40925", "P40926"))
  expect_equivalent(nchar(a[[1]][1]), 47)
  expect_is(a[[2]], 'character')
  expect_equal(length(a[[2]]), 3)
  expect_equal(a[[2]], c("wmelon", "cyt", "mit" ))
  expect_equal(length(a[[3]]), 3)
  expect_equivalent(nchar(a[[3]][1]), 57)
  expect_is(a[[4]], 'matrix')
  expect_equal(dim(a[[4]]), c(3, 57))

  expect_is(b, "fasta")
  expect_equal(length(b), 4)
  expect_is(b[[1]], 'character')
  expect_equal(length(b[[1]]), 3)
  expect_equal(dim(b$ali), c(3,57))
  expect_equal(b$id, c("wmelon", "cyt", "mit" ))
})

