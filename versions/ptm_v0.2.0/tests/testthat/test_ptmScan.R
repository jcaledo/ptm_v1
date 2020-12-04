library(ptm)
context("PTM Scan Tests")

## ---------------------------------------------- ##
#                 Testing p.scan                   #
## ---------------------------------------------- ##
test_that('p.scan() works properly', {

  a <- p.scan('P01009')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 29)
  expect_equal(ncol(a), 4)
  expect_equal(a$up_id[1], 'P01009')
  expect_equal(a$modification[1], 'T35-p')
  expect_equal(a$database[1], 'PSP')
  expect_equal(a$database[22], 'dbPAF')
})

## ---------------------------------------------- ##
#                 Testing ac.scan                  #
## ---------------------------------------------- ##
test_that('ac.scan() works properly', {

  a <- ac.scan('P01009')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 13)
  expect_equal(ncol(a), 4)
  expect_equal(a$up_id[1], 'P01009')
  expect_equal(a$modification[1], 'K159-ac')
  expect_equal(a$database[1], 'PSP')
})

## ---------------------------------------------- ##
#                 Testing me.scan                  #
## ---------------------------------------------- ##
test_that('me.scan() works properly', {

  a <- me.scan('B0I1T2')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 6)
  expect_equal(ncol(a), 4)
  expect_equal(a$modification[1], 'K140-m1')
  expect_equal(a$database[1], 'PSP')

})

## ---------------------------------------------- ##
#                 Testing ub.scan                  #
## ---------------------------------------------- ##
test_that('ub.scan() works properly', {

  a <- ub.scan('B0I1T2')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 36)
  expect_equal(ncol(a), 4)
  expect_equal(a$modification[1], 'K87-ub')
  expect_equal(a$database[1], 'PSP')
})

## ---------------------------------------------- ##
#                 Testing su.scan                  #
## ---------------------------------------------- ##
test_that('su.scan() works properly',{

  a <- su.scan('A6NHL2')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 3)
  expect_equal(ncol(a), 4)
  expect_equal(a$modification[1],'K103-sm')
  expect_equal(a$database[1], 'PSP')
})

## ---------------------------------------------- ##
#                 Testing sni.scan                 #
## ---------------------------------------------- ##
test_that('sni.scan() works properly',{

  a <- sni.scan('P27348')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 3)
  expect_equal(ncol(a), 4)
  expect_equal(a$modification[1],'C25-Sni')
  expect_equal(a$database[1], 'dbPTM')
})

## ---------------------------------------------- ##
#                 Testing ni.scan                  #
## ---------------------------------------------- ##
test_that('ni.scan() works properly',{

  a <- ni.scan('P97427')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 1)
  expect_equal(ncol(a), 4)
  expect_equal(a$modification[1],'Y316-ni')
  expect_equal(a$database[1], 'dbPTM')
})

## ---------------------------------------------- ##
#                 Testing ptm.scan                 #
## ---------------------------------------------- ##
test_that('ptm.scan() works properly', {

  a <- suppressWarnings(ptm.scan('P01009'))
  b <- suppressWarnings(ptm.scan('O14757'))
  c <- suppressWarnings(ptm.scan('G3SB67')) # without described ptm
  d <- suppressWarnings(ptm.scan('Q15796')) # without meto
  # Avoid warnings when a modification is not present in the protein

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 27)
  expect_equal(ncol(a), 15)

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 54)
  expect_equal(ncol(b), 15)

  expect_is(c, 'character')
  expect_true(grepl("Sorry", c))

  expect_is(d, 'data.frame')
  expect_equal(nrow(d), 37)
  expect_equal(ncol(d), 15)
})

## ---------------------------------------------- ##
#                 Testing reg.scan                 #
## ---------------------------------------------- ##
test_that('reg.scan() works properly', {

  a <- reg.scan('O14757')
  b <- reg.scan('P01009')
  c <- suppressWarnings(reg.scan('G3SB67'))  # without described ptm
  d <- reg.scan('Q15796') # without meto

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 6)
  expect_equal(ncol(a), 4)
  expect_equal(a$modification[6], "S296-p")

  expect_is(a, 'data.frame')
  expect_equal(nrow(b), 3)
  expect_equal(ncol(b), 4)
  expect_equal(b$modification[3], "M358-ox")

  expect_is(c, 'character')
  expect_true(grepl("Sorry", c))

  expect_is(d, 'data.frame')
  expect_equal(nrow(d), 12)
  expect_equal(ncol(d), 4)
  expect_true(d$database[1] == 'PSP')

})

## ---------------------------------------------- ##
#                 Testing dis.scan                 #
## ---------------------------------------------- ##
test_that('dis.scan() works properly', {

  a <- dis.scan('O14757')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 2)
  expect_equal(ncol(a), 5)
  expect_equal(a$modification[1], "S280-p")
  expect_match(a$disease[1], 'breast cancer')
})
