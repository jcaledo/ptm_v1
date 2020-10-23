library(ptm)
context("Environment Functions Tests")

## ---------------------------------------------- ##
#               Testing env.extract                #
## ---------------------------------------------- ##
test_that("env.extract() works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- env.extract('P01009', db = 'uniprot', c = 351+24, r = 10, ctr = 'random')
  b <- env.extract('P01009', db = 'uniprot', c = 351+24, r = 10, ctr = 'closest')
  c <- env.extract('P01009', db = 'uniprot', c = 7, r = 14)
  d <- env.extract('P01009', db = 'uniprot', c = 415, r = 20)

  expect_is(a, 'list')
  expect_equal(length(a), 2)
  expect_equal(nchar(a$Positive), 21)
  expect_equal(nchar(a$Positive), nchar(a$Control))

  expect_is(b, 'list')
  expect_equal(length(b), 2)
  expect_equal(nchar(b$Positive), 21)
  expect_false(b$Positive == b$Control)
  expect_true(a$Positive == b$Positive)

  expect_is(c, 'list')
  expect_equal(length(c), 2)
  expect_true('X' %in% strsplit(c$Positive, split = "")[[1]])
  expect_equal(nchar(c$Control), 0)

  expect_is(d, 'list')
  expect_equal(length(d), 2)
  expect_true('X' %in% strsplit(d$Positive, split = "")[[1]])
  expect_equal(nchar(d$Control), 0)
})


test_that("env.extract() works properly",{

  seq1 <- "ARSTVWXXWVTSRAYPILNMSSQQTTWWYYRTGFLIVSTHKRED"
  seq2 <- "ARSTVWXXWVTSRAYPILNMSSQQTTWWYYMRTGFLIVSTHKRED"
  a <- env.extract(seq1,  c = 20, r = 5, ctr = 'random')
  b <- env.extract(seq2, c = 20, r = 5, ctr = 'random')
  c <- env.extract(seq2,  c = 20, r = 5, ctr = 'closest')
  d <- env.extract(seq2,  c = 20, r = 5, ctr = 'random', exclude = 31)

  expect_is(a, 'list')
  expect_equal(length(a), 2)
  expect_equal(nchar(a$Positive), 11)
  expect_gt(nchar(a$Positive), nchar(a$Control))

  expect_is(b, 'list')
  expect_equal(length(b), 2)
  expect_equal(nchar(b$Positive), 11)
  expect_false(b$Positive == b$Control)
  expect_true(a$Positive == b$Positive)

  expect_is(c, 'list')
  expect_equal(length(c), 2)
  expect_true('W' %in% strsplit(c$Control, split = "")[[1]])
  expect_equal(nchar(c$Positive), nchar(c$Control))

  expect_is(d, 'list')
  expect_equal(length(d), 2)
  expect_true('S' %in% strsplit(d$Positive, split = "")[[1]])
  expect_equal(nchar(d$Control), 0)
})


test_that("env.extract() works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- env.extract('P01009', db = 'metosite', c = 351+24, r = 10, ctr = 'random')

  expect_is(a, 'list')
  expect_equal(length(a), 2)
  expect_equal(nchar(a$Positive), 21)
  expect_equal(nchar(a$Positive), nchar(a$Control))
})


## ---------------------------------------------- ##
#               Testing env.matrices               #
## ---------------------------------------------- ##
test_that("env.matrices() works properly",{

  data(hmeto)
  hmeto <- hmeto[1:100, ]
  a <- env.matrices(hmeto$positive)

  expect_is(a, 'list')
  expect_is(a[[1]], 'data.frame')
  expect_equal(nrow(a[[1]]), 100)
  expect_equal(ncol(a[[1]]), 21)
  expect_equal(as.character(a[[1]]$'0'), rep('m', 100))
  expect_is(a[[2]], 'data.frame')
  expect_equal(nrow(a[[2]]), 21)
  expect_equal(ncol(a[[2]]), 21)
  expect_equal(as.vector(a[[2]]$'0'), c(rep(0, 10), 100, rep(0, 10)))
})


## ---------------------------------------------- ##
#               Testing env.Ztest                  #
## ---------------------------------------------- ##
test_that("env.Ztest() works properly",{

  data(hmeto)
  hmeto <- hmeto[1:100, ]
  pos <- env.matrices(hmeto$positive)
  ctr <- env.matrices(hmeto$control)
  t <- env.Ztest(pos = pos[[2]], ctr = ctr[[2]])
  uno <- t[[1]]
  dos <- t[[2]]
  tres <- t[[3]]

  expect_is(uno, 'matrix')
  expect_equal(dim(uno), c(21,21))
  expect_equivalent(uno[, 11], rep(0, 21))
  expect_is(dos, 'data.frame')
  expect_equal(nrow(dos), 31)
  expect_equal(ncol(dos), 4)
  expect_equal(dos$pValue < 0.05,  rep(TRUE, nrow(dos)))
  expect_equal(dos$Z > 1.5, rep(TRUE, nrow(dos)))
  expect_is(tres, 'data.frame')
  expect_equal(nrow(tres), 31)
  expect_equal(ncol(tres), 4)
  expect_equal(tres$pValue < 0.05,  rep(TRUE, nrow(tres)))
  expect_equal(tres$Z < -1.5, rep(TRUE, nrow(tres)))
})


## ---------------------------------------------- ##
#               Testing env.plot                   #
## ---------------------------------------------- ##
# test_that("env.plot() works properly",{
#
#   skip_on_cran()
#   skip_on_travis()
#
#   data(hmeto)
#   hmeto <- hmeto[1:100, ]
#   pos <- env.matrices(hmeto$positive)
#   ctr <- env.matrices(hmeto$control)
#   t <- env.Ztest(pos = pos[[2]], ctr = ctr[[2]])
#   z <- t[[1]]
#
#   png(filename = "myplot.png")
#   myplot <- env.plot(z, aa = 'G', pValue = 0.05)
#   dev.off()
#
#   if (requireNamespace("visualTest", quietly = TRUE)){
#     es <- visualTest::isSimilar(file = "myplot.png",
#                                 fingerprint = visualTest::getFingerprint(file = 'Zplot.png'),
#                                 threshold = 0.1)
#     expect_true(es)
#   }
#
# })
