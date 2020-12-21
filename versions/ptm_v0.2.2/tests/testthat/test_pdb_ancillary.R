library(ptm)
context("PDB Ancillary")

## ---------------------------------------------- ##
#               Testing  pdb.seq                   #
## ---------------------------------------------- ##
test_that('pdb.quaternary works properly', {

  skip_on_cran()
  skip_on_travis()

  a <- pdb.seq('1bpl')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 2)
  expect_equal(ncol(a), 6)
})

## ---------------------------------------------- ##
#            Testing  pdb.quaternary               #
## ---------------------------------------------- ##
test_that('pdb.quaternary works properly', {

  skip_on_cran()
  skip_on_travis()

  a <- pdb.quaternary("./pdb/1u8f.pdb")
  b <- pdb.quaternary("2hhb", keepfiles = TRUE)

  expect_is(a, 'list')
  expect_equal(length(a), 4)
  expect_true(is.matrix(a[[1]]))
  expect_equal(a[[1]][1], 0)
  expect_equal(nchar(a[[2]][1]), 335)
  expect_equal(a[[3]], c("O", "P", "Q", "R"))
  expect_equal(a[[4]], '1u8f')

  expect_is(b, 'list')
  expect_equal(length(b), 4)
  expect_true(is.matrix(b[[1]]))
  expect_equal(b[[1]][1], 0)
  expect_equal(nchar(b[[2]][1]), 141)
  expect_equal(b[[3]], c("A", "B", "C", "D"))
  expect_equal(b[[4]], "2hhb")
  expect_true(file.exists("./2hhb.fa"))
})
file.remove("2hhb.fa")
unlink('./split_chain', recursive = TRUE)

## ---------------------------------------------- ##
#            Testing pdb.chain                     #
## ---------------------------------------------- ##
test_that('pdb.chain works properly',{

  skip_on_cran()
  skip_on_travis()

  a <- pdb.chain("2gls")
  b <- pdb.chain("./pdb/1u8f.pdb", keepfiles = TRUE)


  expect_equal(length(a), 12)
  expect_equal(a, c("A","B","C","D","E","F","G","H","I","J","K","L"))
  expect_equal(length(b), 4)
  expect_true(file.exists('./split_chain/1u8f_P.pdb'))

})
unlink("./split_chain", recursive = TRUE)


## ---------------------------------------------- ##
#            Testing pdb2uniprot                   #
## ---------------------------------------------- ##


test_that('pdb2uniprot works properly', {

  skip_on_cran()
  skip_on_travis()

  a <- pdb2uniprot('2OCC', 'M')

  expect_is(a, 'character')
  expect_match(a, 'P10175')
})

## ---------------------------------------------- ##
#            Testing uniprot2pdb                   #
## ---------------------------------------------- ##


test_that('uniprot2pdb works properly', {

  skip_on_cran()
  skip_on_travis()

  b <- uniprot2pdb('P10175')

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 100)
  expect_equal(ncol(b), 7)
})

## ---------------------------------------------- ##
#             Testing pdb.res                      #
## ---------------------------------------------- ##
test_that('pdb.res() works properly', {

  skip_on_cran()
  skip_on_travis()

  t <- pdb.res(at = 361, up = 'P48163', pdb = '2aw5', chain = 'A')
  f <- pdb.res(at = 15, up = 'P01009', pdb = '3cwm', chain = 'A')

  expect_true(t)
  expect_equal(attributes(t)$pdb_pos, 341)
  expect_false(f)
})

## ---------------------------------------------- ##
#             Testing pdb.pep                      #
## ---------------------------------------------- ##
test_that('pdb.pep() works properly', {

  skip_on_cran()
  skip_on_travis()

  t <- pdb.pep(pep = 'IVKG' , pdb = '2aw5')
  f <- pdb.pep(pep = 'IVKGRASLTQEQ' , pdb = '2aw5')

  expect_true(t)
  expect_equal(attributes(t)$at[1], 325)
  expect_equal(attributes(t)$chain, c("A", "B", "C"))
  expect_false(f)
})

## ---------------------------------------------- ##
#            Testing  pdb.select                   #
## ---------------------------------------------- ##
test_that('pdb.select works properly', {

  skip_on_cran()
  skip_on_travis()

  pdb.select('P04406') -> a
  pdb.select('P0DP23') -> b # CaM, contain Ca2+ in the structure
  pdb.select('P23246') -> c # Only 42.6 % coverage
  suppressWarnings(pdb.select('G3SB67')) -> d # Glutaminase from gorilla (no pdb found)

  expect_is(a, 'list')
  expect_equal(length(a), 2)
  expect_equal(a[[1]], '1u8f')
  expect_equal(a[[2]], 'O')
  expect_gt(attributes(a)$coverage, 0.9)

  expect_is(b, 'list')
  expect_equal(length(b), 2)
  expect_equal(b[[1]], '5jqa')
  expect_equal(b[[2]], 'A')
  expect_gt(attributes(a)$coverage, 0.9)

  expect_is(c, 'list')
  expect_equal(length(c), 2)
  expect_equal(c[[1]], '4wij')
  expect_equal(c[[2]], 'B')
  expect_lt(attributes(c)$coverage, 0.43)

  expect_is(d, 'character')
  expect_equal(d, 'NO PDB FOUND')
})
