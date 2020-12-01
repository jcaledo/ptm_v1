library(ptm)
context("Residue-Focused Functions Tests")

## ----------------------------------------------- ##
#         Testing the function aa.at                #
## ----------------------------------------------- ##
test_that("aa.at() using named arguments with order", {

  t <- aa.at(at=28 , target='P01009' , uniprot=TRUE)
  w <- aa.at(at=500, target='P01009' , uniprot=TRUE)
  x <- aa.at(at=-5, target='P01009' , uniprot=TRUE)

  expect_equal(t, "Q")
  expect_equal(length(t), 1)
  expect_equal(nchar(t), 1)
  expect_equal(length(w), 1)
  expect_gt(nchar(w), 1)
  expect_gt(nchar(x), 1)
})

test_that("aa.at() using named arguments without order", {

  t <- aa.at(target='P01009' , uniprot=TRUE, at=28)
  w <- aa.at(at=500, uniprot=TRUE, target='P01009')
  x <- aa.at(uniprot=TRUE, at=-5, target='P01009')

  expect_equal(t, "Q")
  expect_equal(length(t), 1)
  expect_equal(nchar(t), 1)
  expect_equal(length(w), 1)
  expect_gt(nchar(w), 1)
  expect_gt(nchar(x), 1)
})

test_that("aa.at() using uniprot = TRUE",{
  expect_equal(aa.at(28, 'P01009'), "Q")
  expect_equal(length(aa.at(80, 'P00109')), 1)
  expect_equal(nchar(aa.at(100, 'P01009')), 1)
  expect_equal(length(aa.at(500, 'P01009')), 1)
  expect_gt(nchar(aa.at(500, 'P01009')), 1)
  expect_gt(nchar(aa.at(-5, 'P01009')), 1)
})

test_that("aa.at using uniprot = FALSE",{
  expect_equal(aa.at(6, "MARSSRAM", FALSE), 'R')
  expect_gt(nchar(aa.at(10, "MARSSRAM", FALSE)), 1)
  expect_gt(nchar(aa.at(-6, "MARSSRAM", FALSE)), 1)
})


## ---------------------------------------------- ##
#            Testing the function is.at            #
## ---------------------------------------------- ##
test_that("is.at() works properly", {
  expect_true(is.at(at = 5, target = 'MARTMTRAM', uniprot = FALSE))
  expect_true(is.at(28, 'P01009', 'Q'))
  expect_false(is.at(80, 'P00004', 'R'))
})


## ---------------------------------------------- ##
#        Testing the function renum.pdb            #
## ---------------------------------------------- ##
test_that("the function renum.pdb() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- renum.pdb(pdb = '121P', chain = 'A', uniprot = 'P01112')
  b <- renum.pdb(pdb = "./pdb/1u8f.pdb", chain = 'O', 'P04406')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 189)
  expect_equal(ncol(a), 6)
  expect_equal(a$uni_pos[10], a$pdb_pos[10])
  expect_false(a$uniprot[180] == a$pdb[180])

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 335)
  expect_equal(ncol(b), 6)
  expect_equal(b$uni_pos[10], b$pdb_renum[10])
})


## ---------------------------------------------- ##
#        Testing the function renum.meto           #
## ---------------------------------------------- ##
test_that("the function renum.meto() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- renum.meto(uniprot = 'P01009')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 418)
  expect_equal(ncol(a), 6)
  expect_equal(a$uni_pos[180], a$meto_pos[180] + 24)
  expect_false(a$uniprot[10] == a$meto[10])

})


## ---------------------------------------------- ##
#            Testing the function renum            #
## ---------------------------------------------- ##
test_that('renum() works properly', {

  skip_on_cran()
  skip_on_travis()

  a <- renum(up_id = 'P01009', pos = 60,
             from = 'uniprot', to = 'metosite')
  b <- renum(up_id = 'P01009', pos = 36,
             from = 'metosite', to = 'uniprot')
  c <- renum(up_id = 'P01009', pos = 60,
             from = 'uniprot' , to = 'pdb',
             pdb = '1ATU', chain = 'A')
  d <- renum(up_id = 'P01009', pos = 16,
             from = 'pdb' , to = 'uniprot',
             pdb = '1ATU', chain = 'A')

  expect_is(a, 'numeric')
  expect_equal(a, 36)
  expect_is(b, 'numeric')
  expect_equal(b, 60)
  expect_is(c, 'numeric')
  expect_equal(c, 16)
  expect_is(d, 'numeric')
  expect_equal(d, 60)
})

