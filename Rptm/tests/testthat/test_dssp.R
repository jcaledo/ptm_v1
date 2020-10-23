library(ptm)
context("DSSP manipulations")


## ---------------------------------------------- ##
#             Testing parse.dssp                   #
## ---------------------------------------------- ##
test_that("parse.dssp() works properly", {

  skip_on_cran()
  skip_on_travis()

  compute.dssp('3cwm')
  a <- parse.dssp('3cwm.dssp')

  expect_equal(nrow(a), 369)
  expect_equal(ncol(a), 8)
  expect_equal(sum(c('B', 'C', 'E', 'G', 'H', 'S', 'T') %in% a$ss), 7)
})



## ---------------------------------------------- ##
#            Testing compute.dssp                  #
## ---------------------------------------------- ##

test_that("compute.dssp() works properly", {

  skip_on_cran()
  skip_on_travis()

  compute.dssp(pdb = "./pdb/1U8F.pdb", destfile = "./pdb/")
  expect_true(file.exists('./pdb/1U8F.dssp'))
  b <- parse.dssp("./pdb/1U8F.dssp")

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 1332)
  expect_equal(ncol(b), 8)
})

## ---------------------------------------------- ##
#               Testing mkdssp                     #
## ---------------------------------------------- ##

test_that("mkdssp() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- mkdssp('2xqq', method = 'ptm')
  b <- mkdssp('2xqq', method = 'bio3d')
  c <- mkdssp('./pdb/1u8f.pdb', method = 'ptm')
  d <- mkdssp('./pdb/1u8f.pdb', method = 'bio3d')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 371)
  expect_equal(ncol(a), 8)
  expect_true('ss' %in% names(a))
  expect_gt(sum(a$sasa), 10000)

  expect_is(b, 'sse')
  expect_equal(length(b), 9)
  expect_gte(length(b$sse), 371)
  expect_gt(sum(b$acc), 10000)

  expect_is(c, 'data.frame')
  expect_equal(nrow(c), 1332)
  expect_equal(ncol(c), 8)
  expect_true('ss' %in% names(c))
  expect_gt(sum(c$sasa), 40000)

  expect_is(d, 'sse')
  expect_equal(length(d), 9)
  expect_gte(length(d$sse), 1332)
  expect_gt(sum(d$acc), 40000)
})
