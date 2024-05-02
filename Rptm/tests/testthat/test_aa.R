library(ptm)
context("Residue-Focused Functions Tests")

## ----------------------------------------------- ##
#         Testing the function aa.at                #
## ----------------------------------------------- ##
test_that("aa.at() using named arguments with order", {

  skip_on_cran()
  skip_on_travis()

  a <- aa.at(at=28 , target='P01009' , uniprot=TRUE)
  b <- aa.at(at=500, target='P01009' , uniprot=TRUE)

  if (!is.null(a)){
    expect_is(a, 'character')
    expect_equal(a, "Q")
  }
  expect_is(b, 'NULL')
})

test_that("aa.at() using named arguments without order", {

  skip_on_cran()
  skip_on_travis()

  a <- aa.at(target='P01009' , uniprot=TRUE, at=28)
  b <- aa.at(at=500, uniprot=TRUE, target='P01009')

  if (!is.null(a)){
    expect_is(a, 'character')
    expect_equal(a, "Q")
    expect_equal(length(a), 1)
    expect_equal(nchar(a), 1)
  }
  expect_is(b, 'NULL')

})

test_that("aa.at using uniprot = FALSE",{

  expect_equal(aa.at(6, "MARSSRAM", FALSE), 'R')
  expect_is(aa.at(10, "MARSSRAM", FALSE), 'NULL')

})


## ---------------------------------------------- ##
#            Testing the function is.at            #
## ---------------------------------------------- ##
test_that("is.at() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- is.at(at = 5, target = 'MARTMTRAM', uniprot = FALSE)
  b <- is.at(28, 'P01009', aa = 'Q')
  c <- is.at(at =5, target = 'MARTMTRAM', aa = "P", uniprot = FALSE)
  d <- is.at(28, 'P01009', aa = 'P')
  e <- is.at(28, 'P010091')

  expect_is(a, 'logical')
  expect_true(a)

  if (!is.null(b)){
    expect_is(b, 'logical')
    expect_true(b)
  }

  if (!is.null(c)){
    expect_is(c, 'logical')
    expect_false(c)
  }

  if (!is.null(d)){
    expect_is(d, 'logical')
    expect_false(d)
  }

  expect_is(e, 'NULL')

})

## ---------------------------------------------- ##
#        Testing the function aa.comp              #
## ---------------------------------------------- ##
test_that("the function aa.comp() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- suppressWarnings(aa.comp("ACEEAGRKDNW", uniprot = FALSE))
  b <- aa.comp("P01009")
  c <- aa.comp('P010091')

  expect_is(a, 'list')
  expect_equal(dim(a[[1]]), c(20, 6))

  if (!is.null(b)){
    expect_is(b, 'list')
    expect_equal(dim(b[[1]]), c(20, 6))
    expect_equal(sum(b[[1]]$observed), 418)
    expect_equal(length(b[[2]]), 9)
    expect_lt(b[[2]]$p.value, 0.5)
  }

  expect_is(c, 'NULL')
})


## ---------------------------------------------- ##
#        Testing the function renum.pdb            #
## ---------------------------------------------- ##
test_that("the function renum.pdb() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- renum.pdb(pdb = '121P', chain = 'A', uniprot = 'P01112')
  b <- renum.pdb(pdb = "./pdb/1u8f.pdb", chain = 'O', 'P04406')
  c <- renum.pdb(pdb = 'xxxx', chain = 'X', uniprot = 'P010091')

  if (!is.null(a)){
    expect_is(a, 'data.frame')
    expect_equal(nrow(a), 189)
    expect_equal(ncol(a), 6)
    expect_equal(a$uni_pos[10], a$pdb_pos[10])
    expect_false(a$uniprot[180] == a$pdb[180])
  }

  if (!is.null(b)){
    expect_is(b, 'data.frame')
    expect_equal(nrow(b), 335)
    expect_equal(ncol(b), 6)
    expect_equal(b$uni_pos[10], b$pdb_renum[10])
  }

  expect_is(c, 'NULL')
})


## ---------------------------------------------- ##
#        Testing the function renum.meto           #
## ---------------------------------------------- ##
test_that("the function renum.meto() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- renum.meto(uniprot = 'P01009')
  b <- renum.meto(uniprot = 'P010091')

  if (!is.null(a)){
    expect_is(a, 'data.frame')
    expect_equal(nrow(a), 418)
    expect_equal(ncol(a), 6)
    expect_equal(a$uni_pos[180], a$meto_pos[180] + 24)
    expect_false(a$uniprot[10] == a$meto[10])
  }
  expect_is(b, 'NULL')
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
  e <- renum(up_id = 'P010091', pos = 60,
             from = 'uniprot', to = 'metosite')

  expect_is(a, 'numeric')
  expect_equal(a, 36)
  expect_is(b, 'numeric')
  expect_equal(b, 60)
  if (!is.null(c)){
    expect_is(c, 'numeric')
    expect_equal(c, 16)
  }
  if (!is.null(d)){
    expect_is(d, 'numeric')
    expect_equal(d, 60)
  }
  expect_is(e, 'NULL')
})

