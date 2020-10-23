library(ptm)
context("ACC and SASA Computations")

## ---------------------------------------------- ##
#             Testing acc.dssp                     #
## ---------------------------------------------- ##
test_that("acc.dssp() works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- acc.dssp('./pdb/1U8F.pdb')
  b <- acc.dssp('./pdb/1U8F.pdb', aa = 'M')
  c <- acc.dssp('./pdb/1U8F.pdb', dssp = 'mkdssp')
  d <- acc.dssp('./pdb/1U8F.pdb', dssp = 'mkdssp', aa = 'S')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 1332)
  expect_equal(ncol(a), 11)
  expect_true('P' %in% unique(a$chain))

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 36)
  expect_equal(ncol(b), 11)
  expect_true('Q' %in% unique(b$chain))

  expect_is(c, 'data.frame')
  expect_equal(nrow(c), 1332)
  expect_equal(ncol(c), 11)
  expect_true('P' %in% unique(c$chain))

  expect_is(d, 'data.frame')
  expect_equal(nrow(d), 84)
  expect_equal(ncol(d), 11)
  expect_true('O' %in% unique(d$chain))

  a <- acc.dssp('6lu7')
  b <- acc.dssp('6lu7', aa = 'M')
  c <- acc.dssp('6lu7', dssp = 'mkdssp')
  d <- acc.dssp('6lu7', dssp = 'mkdssp', aa = 'K')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 309)
  expect_equal(ncol(a), 11)
  expect_true('A' %in% unique(a$chain))

  expect_is(a, 'data.frame')
  expect_equal(nrow(b), 10)
  expect_equal(ncol(b), 11)
  expect_true('A' %in% unique(a$chain))

  expect_is(c, 'data.frame')
  expect_equal(nrow(c), 309)
  expect_equal(ncol(c), 11)
  expect_true('A' %in% unique(a$chain))

  expect_is(d, 'data.frame')
  expect_equal(nrow(d), 11)
  expect_equal(ncol(d), 11)
  expect_true('A' %in% unique(a$chain))

  a <- acc.dssp('1h9d') # pdb contain non-protein chains

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 504)
  expect_equal(ncol(a), 11)

})

## ---------------------------------------------- ##
#             Testing get.area                     #
## ---------------------------------------------- ##
test_that("get.area() works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- get.area('3cwm')
  b <- get.area("./pdb/1U8F.pdb")
  c <- get.area("./pdb/1U8F.pdb", keepfiles = TRUE)

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 2929)
  expect_equal(ncol(a), 5)
  expect_true('areaenergy' %in% names(a))

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 10092)
  expect_equal(ncol(b), 5)
  expect_true('resno' %in% names(b))

  expect_is(c, 'data.frame')
  expect_equal(nrow(c), 10092)
  expect_equal(ncol(c), 5)
  expect_true('resid' %in% names(c))
  expect_true(file.exists('./pdb/1U8F_getarea.Rda'))
  expect_true(file.exists('./pdb/1U8F_getarea.txt'))
})
file.remove('./pdb/1U8F_getarea.Rda')
file.remove('./pdb/1U8F_getarea.txt')


## ---------------------------------------------- ##
#                 Testing dpx                      #
## ---------------------------------------------- ##
test_that("dpx() works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- dpx('3cwm')
  b <- dpx('./pdb/1U8F.pdb')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 2929)
  expect_equal(ncol(a), 9)
  expect_gt(sum(a$dpx), 1)

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 10092)
  expect_equal(ncol(b), 9)
  expect_gt(sum(b$dpx), 1)

})


## ---------------------------------------------- ##
#                 Testing atom.dpx                 #
## ---------------------------------------------- ##
test_that("atom.dpx() works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- atom.dpx('1hxp')
  b <- atom.dpx('./pdb/1U8F.pdb')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 5354)
  expect_equal(ncol(a), 9)
  expect_gt(sum(a$delta_dpx), 1)

  expect_is(b, 'data.frame')
  expect_equal(nrow(b),10092)
  expect_equal(ncol(b), 9)
  expect_gt(sum(b$delta_dpx), 1)
})

## ---------------------------------------------- ##
#                 Testing res.dpx                  #
## ---------------------------------------------- ##
test_that("res.dpx() works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- res.dpx('1hxp')
  b <- res.dpx('./pdb/1U8F.pdb')
  c <- res.dpx('1cll', aa = 'M')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 669)
  expect_equal(ncol(a), 10)
  expect_gt(sum(a$max_dpx_complex), sum(a$max_dpx_chain))

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 1332)
  expect_equal(ncol(b), 10)
  expect_gt(sum(b$max_dpx_complex), sum(b$max_dpx_chain))

  expect_is(c, 'data.frame')
  expect_equal(nrow(c),9)
  expect_equal(ncol(c), 10)
  expect_gte(c$max_dpx_chain[1], c$min_dpx_chain[1])
})

## ---------------------------------------------- ##
#             Testing stru.part                    #
## ---------------------------------------------- ##
test_that("stru.part() works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- stru.part(pdb = '2xqq')
  b <- stru.part(pdb = '2xqq', cutoff = 0.05)
  c <- stru.part('./pdb/1u8f.pdb')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 371)
  expect_equal(ncol(a), 9)
  expect_equal(as.numeric(table(a$str)['surface']), 96)
  expect_equal(as.numeric(table(a$str)['interior']), 104)
  expect_equal(as.numeric(table(a$str)['rim']), 68)
  expect_equal(as.numeric(table(a$str)['support']), 34)
  expect_equal(as.numeric(table(a$str)['core']), 69)

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 371)
  expect_equal(ncol(b), 9)
  expect_gt(as.numeric(table(b$str)['surface']), as.numeric(table(a$str)['surface']))
  expect_lt(as.numeric(table(b$str)['interior']), as.numeric(table(a$str)['interior']))
  expect_gt(as.numeric(table(b$str)['rim']), as.numeric(table(a$str)['rim']))
  expect_lt(as.numeric(table(b$str)['support']), as.numeric(table(a$str)['support']))

  expect_is(c, 'data.frame')
  expect_equal(nrow(c), 1332)
  expect_equal(ncol(c), 9)
  expect_equal(as.numeric(table(c$str)['surface']), 358)
  expect_equal(as.numeric(table(c$str)['interior']), 615)
  expect_equal(as.numeric(table(c$str)['rim']), 92)
  expect_equal(as.numeric(table(c$str)['support']), 114)
  expect_equal(as.numeric(table(c$str)['core']), 153)
})
