library(ptm)
context("GO-related")

## ---------------------------------------------- ##
#              Testing search.go                   #
## ---------------------------------------------- ##
test_that("search.go() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- search.go(query = 'oxidative stress')
  b <- search.go(query = 'methionine')

  expect_is(a, 'data.frame')
  expect_gte(nrow(a), 300)
  expect_equal(ncol(a), 5)
  expect_true('GO:0070994' %in% a$GO_id)
})

## ---------------------------------------------- ##
#              Testing term.go                     #
## ---------------------------------------------- ##
test_that("term.go() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- term.go('GO:0034599')
  b <- term.go('GO:0005886', children = TRUE)

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 1)
  expect_equal(ncol(a), 7)
  expect_equal(a$term_name, "cellular response to oxidative stress")

  expect_is(b, 'list')
  expect_is(b[[1]], 'data.frame')
  expect_equal(nrow(b[[1]]), 1)
  expect_equal(ncol(b[[1]]), 7)
  expect_equal(b[[1]]$term_name, "plasma membrane")
  expect_is(b[[2]], 'data.frame')
  expect_equal(ncol(b[[2]]), 2)
})


## ---------------------------------------------- ##
#              Testing get.go                     #
## ---------------------------------------------- ##
test_that("get.go() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- get.go(id = "P01009")
  b <- get.go(id = 'P01009', filter = FALSE)
  c <- get.go(id = 'P04406', format = 'string')
  d <- get.go(id = 'P04406', filter = FALSE, format = 'string')
  e <- get.go(id = 'P00367') # requires removing spurious rows
  f <- get.go(id = 'Q14687') # no GO terms found

  expect_is(a, 'data.frame')
  expect_gte(nrow(a), 25)
  expect_equal(ncol(a), 5)

  expect_is(b, 'data.frame')
  expect_gte(nrow(b), 76)
  expect_equal(ncol(b), 8)

  expect_is(c, 'character')
  expect_gte(nchar(c), 460)
  expect_is('d', "character")
  expect_gte(nchar(d), 1100)

  expect_is(e, 'data.frame')
  expect_gte(nrow(e), 15)
  expect_equal(ncol(e), 5)

  expect_is(f, 'character')
  expect_true(f == "Sorry, no GO terms found for the Q14687 entry")
})

## ---------------------------------------------- ##
#                Testing bg.go                     #
## ---------------------------------------------- ##
test_that("bg.go() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- bg.go(ids = "./go/id_set.txt")
  b <-  bg.go(ids = c("Q13015", "Q14667", "P08575", "Q5JSZ5", "P13196", "H7C4H7"))

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 6)
  expect_equal(ncol(a), 2)

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 6)
  expect_equal(ncol(b), 2)

  expect_equal(a, b)

})

## ---------------------------------------------- ##
#           Testing  hdfisher.go                  #
## ---------------------------------------------- ##
test_that(" hdfisher.go() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- hdfisher.go(target = c('Q14667', 'Q5JSZ5'),
                   background = bg.go(c("Q13015", "Q14667", "P08575", "Q5JSZ5", "P13196")),
                   query = 'extracellular')

  expect_is(a, 'list')
  expect_is(a[[1]], 'matrix')
  expect_is(a[[2]], 'numeric')
  expect_true(attributes(a)$query == 'extracellular')
})


## ---------------------------------------------- ##
#              Testing net.go                      #
## ---------------------------------------------- ##
test_that("net.go() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- net.go(data = "./go/id_set.txt", threshold = 0.5)
  b <- net.go(data = "./go/id_set.Rda", threshold = 0.5)

  expect_is(a, 'list')
  expect_is(a[[1]], 'matrix')
  expect_equal(dim(a[[1]]), c(6,6))
  expect_true(!isSymmetric(a[[1]]))
  expect_is(a[[2]], 'matrix')
  expect_equal(dim(a[[2]]), c(6,6))
  expect_true(isSymmetric(a[[2]]))
  expect_is(a[[3]], 'character')
  expect_equal(length(a[[3]]), 6)
  expect_is(a[[4]], 'matrix')
  expect_equal(ncol(a[[4]]), 2)

  expect_is(b, 'list')
  expect_is(b[[1]], 'matrix')
  expect_equal(dim(b[[1]]), c(6,6))
  expect_true(!isSymmetric(b[[1]]))
  expect_is(b[[2]], 'matrix')
  expect_equal(dim(b[[2]]), c(6,6))
  expect_true(isSymmetric(b[[2]]))
  expect_is(b[[3]], 'character')
  expect_equal(length(b[[3]]), 6)
  expect_is(b[[4]], 'matrix')
  expect_equal(ncol(b[[4]]), 2)
})

## ---------------------------------------------- ##
#              Testing gorilla                     #
## ---------------------------------------------- ##
test_that("gorila() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- gorilla(target = './go/GOvivo.txt', spe = 'Homo sapiens')
  b <- gorilla(target = './go/GOvivo.txt', db = 'all', spe = 'Homo sapiens')

  c <- gorilla(target = './go/GOvivo.txt',
               background = './go/GObackground.txt',
               mode = 'hg', spe = 'Homo sapiens')
  d <- gorilla(target = './go/GOvivo.txt',
               background = './go/GObackground.txt',
               db = 'func', mode = 'hg', spe = 'Homo sapiens')
  e <- gorilla(target = './go/GOvivo.txt',
               background = './go/GObackground.txt',
               db = 'comp', mode = 'hg', spe = 'Homo sapiens')
  f <- gorilla(target = './go/GOvivo.txt',
               background = './go/GObackground.txt',
               db = 'all', mode = 'hg', spe = 'Homo sapiens')

  expect_is(a, 'data.frame')
  expect_gt(nrow(a), 100)
  expect_gt(ncol(a), 8)

  expect_is(b, 'list')
  expect_equal(length(b), 3)
  expect_is(b[[1]], 'data.frame')
  expect_gt(nrow(b[[1]]), 10)
  expect_equal(ncol(b[[1]]), 10)
  expect_is(b[[2]], 'data.frame')
  expect_gt(nrow(b[[2]]), 10)
  expect_equal(ncol(b[[2]]), 10)
  expect_is(b[[3]], 'data.frame')
  expect_gt(nrow(b[[3]]), 10)
  expect_equal(ncol(b[[3]]), 10)

  expect_is(c, 'data.frame')
  expect_gt(nrow(c), 100)
  expect_gt(ncol(c), 8)

  expect_is(d, 'data.frame')
  expect_gt(nrow(d), 50)
  expect_gt(ncol(d), 8)

  expect_is(e, 'data.frame')
  expect_gt(nrow(e), 100)
  expect_gt(ncol(e), 8)

  expect_is(f, 'list')
  expect_equal(length(f), 3)
  expect_is(f[[1]], 'data.frame')
  expect_gt(nrow(f[[1]]), 10)
  expect_equal(ncol(f[[1]]), 10)
  expect_is(f[[2]], 'data.frame')
  expect_gt(nrow(f[[2]]), 10)
  expect_equal(ncol(f[[2]]), 10)
  expect_is(f[[3]], 'data.frame')
  expect_gt(nrow(f[[3]]), 10)
  expect_equal(ncol(f[[3]]), 10)

})

