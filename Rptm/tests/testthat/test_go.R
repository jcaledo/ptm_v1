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
  c <- search.go(query = 'xxxxxx')

  if (!is.null(a)){
    expect_is(a, 'data.frame')
    expect_gte(nrow(a), 300)
    expect_equal(ncol(a), 5)
    expect_true('GO:0070994' %in% a$GO_id)
  }

  if (!is.null(b)){
    expect_is(b, 'data.frame')
    expect_gte(nrow(b), 100)
    expect_equal(ncol(b), 5)
    expect_true('GO:0015821' %in% b$GO_id)
  }

  expect_is(c, 'NULL')

})

## ---------------------------------------------- ##
#              Testing term.go                     #
## ---------------------------------------------- ##
test_that("term.go() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- term.go('GO:0034599')
  b <- term.go('GO:0005886', children = TRUE)
  c <- term.go('xxxxxxx')

  if (!is.null(a)){
    expect_is(a, 'data.frame')
    expect_equal(nrow(a), 1)
    expect_equal(ncol(a), 7)
    expect_equal(a$term_name, "cellular response to oxidative stress")
  }

  if (!is.null(b)){
    expect_is(b, 'list')
    expect_is(b[[1]], 'data.frame')
    expect_equal(nrow(b[[1]]), 1)
    expect_equal(ncol(b[[1]]), 7)
    expect_equal(b[[1]]$term_name, "plasma membrane")
    expect_is(b[[2]], 'data.frame')
    expect_equal(ncol(b[[2]]), 2)
  }

  expect_is(c, 'NULL')
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
  g <- get.go(id = "P010091") # dummy protein ID
  h <- get.go("AAA58698", filter = FALSE) # currently P06748
  i <- get.go("AAA58698")

  if (!is.null(a)){
    expect_is(a, 'data.frame')
    expect_gte(nrow(a), 20)
    expect_equal(ncol(a), 5)
  }

  if (!is.null(b)){
    expect_is(b, 'data.frame')
    expect_gte(nrow(b), 70)
    expect_equal(ncol(b), 8)
  }

  if (!is.null(c)){
    expect_is(c, 'character')
    expect_gte(nchar(c), 400)
  }

  if (!is.null(d)){
    expect_is(d, "character")
    expect_gte(nchar(d), 1100)
  }

  if (!is.null(e)){
    expect_is(e, 'data.frame')
    expect_gte(nrow(e), 15)
    expect_equal(ncol(e), 5)
  }

  expect_is(f, 'NULL')
  expect_is(g, 'NULL')
  expect_is(h, 'NULL')
  expect_is(i, 'NULL')

})

## ---------------------------------------------- ##
#                Testing bg.go                     #
## ---------------------------------------------- ##
test_that("bg.go() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- bg.go(ids = "./go/id_set.txt")
  b <-  bg.go(ids = c("Q13015", "Q14667", "P08575", "Q5JSZ5", "P13196", "H7C4H7"))

  if (!is.null(a)){
    expect_is(a, 'data.frame')
    expect_equal(nrow(a), 6)
    expect_equal(ncol(a), 2)
  }

  if (!is.null(b)){
    expect_is(b, 'data.frame')
    expect_equal(nrow(b), 6)
    expect_equal(ncol(b), 2)
    expect_equal(a, b)
  }

})

## ---------------------------------------------- ##
#           Testing  hdfisher.go                  #
## ---------------------------------------------- ##
test_that(" hdfisher.go() works properly", {

  skip_on_cran()
  skip_on_travis()

  backg <- bg.go(c("Q13015", "Q14667", "P08575", "Q5JSZ5", "P13196"))

  if (!is.null(backg)){
    a <- hdfisher.go(target = c('Q14667', 'Q5JSZ5'),
                     background = backg,  query = 'extracellular')
  } else {
    a <- NULL
  }

  if (!is.null(backg)){
    b <- hdfisher.go(target = c('Q14667', 'Q5JSZ5'),
                     background = backg,
                     query = 'xxxxxxx')
  } else {
    b <- NULL
  }

  if (!is.null(a)){
    expect_is(a, 'list')
    expect_is(a[[1]], 'matrix')
    expect_is(a[[2]], 'numeric')
    expect_true(attributes(a)$query == 'extracellular')
  }

  expect_is(b, 'NULL')
})


## ---------------------------------------------- ##
#              Testing net.go                      #
## ---------------------------------------------- ##
test_that("net.go() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- net.go(data = "./go/id_set.txt", threshold = 0.1)
  b <- net.go(data = "./go/id_set.Rda", threshold = 0.1)
  c <- net.go(data = "./go/id_set_dummy.txt", threshold = 0.1)

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

  expect_is(c, 'list')
  expect_is(c[[1]], 'matrix')
  expect_equal(dim(c[[1]]), c(5,5))
  expect_true(!isSymmetric(c[[1]]))
  expect_is(c[[2]], 'matrix')
  expect_equal(dim(c[[2]]), c(5,5))
  expect_true(isSymmetric(c[[2]]))
  expect_is(c[[3]], 'character')
  expect_equal(length(c[[3]]), 6)
  expect_is(c[[4]], 'matrix')
  expect_equal(ncol(c[[4]]), 2)
})

