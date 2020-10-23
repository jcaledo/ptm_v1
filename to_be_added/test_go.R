library(ptm)
context("GO-related")

## ---------------------------------------------- ##
#              Testing search.go                   #
## ---------------------------------------------- ##
test_that("search.go() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- search.go(query = 'oxidative stress')
  b <- search.go(query = 'apoptosis')


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
})

## ---------------------------------------------- ##
#              Testing gorilla                     #
## ---------------------------------------------- ##
test_that("gorilla() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- gorilla(target = './go/GOvivo.txt', mode = 'mhg', db = 'proc', pvalue = 0.001, species = 'Homo sapiens')
  b <- gorilla(target = './go/GOvivo.txt', mode = 'mhg', db = 'all', pvalue = 0.001)
  c <- gorilla(target = './go/GOvivo.txt', background = './go/GObackground.txt', mode = 'hg', db = 'all')

  expect_is(a, 'data.frame')
  expect_gt(nrow(a), 100)
  expect_equal(ncol(a), 10)

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

  expect_is(c, 'list')
  expect_equal(length(c), 3)
  expect_is(c[[1]], 'data.frame')
  expect_gt(nrow(c[[1]]), 10)
  expect_equal(ncol(c[[1]]), 10)
  expect_is(c[[2]], 'data.frame')
  expect_gt(nrow(c[[2]]), 10)
  expect_equal(ncol(c[[2]]), 10)
  expect_is(c[[3]], 'data.frame')
  expect_gt(nrow(c[[3]]), 10)
  expect_equal(ncol(c[[3]]), 10)

})

## ---------------------------------------------- ##
#              Testing net.go                      #
## ---------------------------------------------- ##
test_that("net.go() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- net.go(data = "./go/sample_vivo.txt", threshold = 0.5)
  b <- net.go(data = "./go/sample_vv.Rda", threshold = 0.5)

  expect_is(a, 'list')
  expect_is(a[[1]], 'matrix')
  expect_equal(dim(a[[1]]), c(100,100))
  expect_true(!isSymmetric(a[[1]]))
  expect_is(a[[2]], 'matrix')
  expect_equal(dim(a[[2]]), c(100,100))
  expect_true(isSymmetric(a[[2]]))
  expect_is(a[[3]], 'data.frame')
  expect_equal(nrow(a[[3]]), 100)
  expect_is(a[[4]], 'matrix')
  expect_equal(ncol(a[[4]]), 2)

  expect_is(b, 'list')
  expect_is(b[[1]], 'matrix')
  expect_equal(dim(b[[1]]), c(150,150))
  expect_true(!isSymmetric(b[[1]]))
  expect_is(b[[2]], 'matrix')
  expect_equal(dim(b[[2]]), c(150,150))
  expect_true(isSymmetric(b[[2]]))
  expect_is(b[[3]], 'data.frame')
  expect_equal(nrow(b[[3]]), 150)
  expect_is(b[[4]], 'matrix')
  expect_equal(ncol(b[[4]]), 2)
})

