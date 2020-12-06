library(ptm)
context("Sequence Recovering Functions Tests")

## ---------------------------------------------- ##
#               Testing get.seq                    #
## ---------------------------------------------- ##
test_that("get.seq() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- get.seq('P01009')
  t <- get.seq('P01009', as.string = FALSE)

  expect_is(a, 'character')
  expect_equal(nchar(a), 418)
  expect_is(t, 'list')
  expect_equal(length(t[[1]]), 418)
  expect_equal(t[[1]][28], 'Q')
})

test_that("get.seq() works properly with MetOSite",{

  skip_on_cran()
  skip_on_travis()

  a <- get.seq('P01009', db = 'metosite')
  t <- get.seq('P01009', db = 'metosite', as.string = FALSE)

  expect_is(a, 'character')
  expect_equal(nchar(a), 394)
  expect_is(t, 'list')
  expect_equal(length(t[[1]]), 394)
  expect_equal(t[[1]][28], 'P')
})

test_that("get.seq() works properly with PDB",{

  skip_on_cran()
  skip_on_travis()

  a <- get.seq('1u8f', db = 'pdb')
  b <- get.seq('1u8f:O', db = 'pdb')
  c <- get.seq('1u8f:P', db = 'pdb')
  d <- get.seq('1u8f:R', db = 'pdb')
  e <- get.seq('1u8f:Q', db = 'pdb', as.string = FALSE)

  expect_is(a, "character")
  expect_equal(nchar(a), 335)
  expect_is(b, "character")
  expect_equal(nchar(b), 335)
  expect_is(c, "character")
  expect_equal(nchar(c), 335)
  expect_is(d, "character")
  expect_equal(nchar(d), 335)
  expect_is(e, "list")
  expect_is(e[[1]], 'character')
  expect_equal(length(e[[1]]), 335)
})


test_that('get.seq() works properly with KEGG', {

  skip_on_cran()
  skip_on_travis()

  a <- get.seq('hsa:5265', 'kegg-aa')
  b <- get.seq('hsa:5265', db = 'kegg-aa', as.string = FALSE)
  c <- get.seq('hsa:5265', 'kegg-nt')
  d <- get.seq('hsa:5265', db = 'kegg-nt', as.string = FALSE)
  e <- get.seq('eih:ECOK1_2344', 'kegg-aa')

  expect_is(a, 'character')
  expect_equivalent(nchar(a), 418)
  expect_is(b, 'list')
  expect_is(b[[1]], 'character')
  expect_equal(length(b[[1]]), 418)
  expect_is(c, 'character')
  expect_equivalent(nchar(c), 1257)
  expect_is(d, 'list')
  expect_is(d[[1]], 'character')
  expect_equal(length(d[[1]]), 1257)
  expect_is(e, 'character')
})


## ---------------------------------------------- ##
#             Testing prot2codon                   #
## ---------------------------------------------- ##
test_that("prot2codon() works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- prot2codon("P01009")
  b <- prot2codon('1ATU', 'A')
  c <- prot2codon("./pdb/1U8F.pdb", chain = 'O')
  d <- prot2codon('1A00', chain = 'A')
  e <- prot2codon('1F8A', 'B')

  expect_is(a, 'data.frame')
  expect_equal(sum(table(a$check)), nrow(a))
  expect_equal(ncol(a), 6)

  expect_is(b, 'data.frame')
  expect_equal(sum(table(b$check)), nrow(b))
  expect_equal(ncol(b), 6)

  expect_is(c, 'data.frame')
  expect_equal(sum(table(c$check)), nrow(c))
  expect_equal(ncol(c), 6)

  expect_is(d, 'data.frame')
  expect_equal(nrow(d), 141)
  expect_equal(ncol(d), 6)
  expect_equal(sum(table(d$check)), nrow(d))

  # expect_is(e, 'data.frame')
  # expect_equal(nrow(e), 167)
  # expect_equal(ncol(e), 6)
  # expect_equal(sum(table(e$check)), nrow(e) - 4)

})


## ---------------------------------------------- ##
#               Testing id.mapping                 #
## ---------------------------------------------- ##
test_that('id.mapping() works properly when pdb -> uniprot', {

  skip_on_cran()
  skip_on_travis()

  a <- id.mapping('3cwm', 'pdb','uniprot')
  b <- id.mapping('2occ', 'pdb', 'uniprot')

  expect_equal(a, 'P01009')
  expect_is(b, 'character')
  expect_equal(length(b), 13)
  expect_true('P00415' %in% b)
})


test_that('id.mapping() works properly when up -> pdb', {

  skip_on_cran()
  skip_on_travis()

  a <- id.mapping('P01009', 'uniprot', 'pdb')
  b <- id.mapping('P00415', 'uniprot', 'pdb')

  expect_is(a, 'character')
  expect_true('3CWM' %in% a)
  expect_is(b, 'character')
  expect_true('2OCC' %in% b)
})


test_that('id.mapping() works properly when up -> kegg', {

  skip_on_cran()
  skip_on_travis()

  up_kegg <- id.mapping('P01009', 'uniprot', 'kegg')

  expect_equivalent(up_kegg, 'hsa:5265')
})


test_that('id.mapping() works properly when kegg -> up', {

  skip_on_cran()
  skip_on_travis()

  kegg_up <- id.mapping('hsa:5265', 'kegg', 'uniprot')

  expect_is(kegg_up, 'character')
  expect_true("P01009" %in% kegg_up)

})


test_that('id.mapping() works properly when pdb -> kegg', {

  skip_on_cran()
  skip_on_travis()

  a <- id.mapping('3cwm', 'pdb','kegg')
  b <- id.mapping('2occ', 'pdb', 'kegg')

  expect_is(a, 'list')
  expect_is(a[[1]], 'character')
  expect_equivalent(a[[1]], "hsa:5265")
  expect_is(b, 'list')
  expect_is(b[[3]], 'character')
  expect_true("bta:281090" %in% b[[3]])
})


test_that('id.mapping() works properly when kegg <- pdb', {

  skip_on_cran()
  skip_on_travis()

  a <- id.mapping('hsa:5265', 'kegg','pdb')
  b <- id.mapping('bta:281090', 'kegg','pdb')

  expect_is(a, 'character')
  expect_true('3CWM' %in% a)
  expect_is(b, 'character')
  expect_true('2OCC' %in% b)
})


## ---------------------------------------------- ##
#            Testing id.features                   #
## ---------------------------------------------- ##
test_that("id.features() works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- id.features('P01009')
  b <- id.features('P01009', features = 'ec,keywords,database(PDB)')

  expect_is(a, 'list')
  expect_equal(length(a), 5)
  expect_equal(a$Status, "reviewed")
  expect_is(b, 'list')
  expect_equal(length(b), 8)
  expect_equal(a$Organism, b$Organism)
})

## ---------------------------------------------- ##
#           Testing species.mapping                #
## ---------------------------------------------- ##
test_that("species.mapping works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- species.mapping('P01009')
  b <- species.mapping('2OCC', db = 'pdb')

  expect_is(a, 'character')
  expect_equal(a, "Homo sapiens")
  expect_is(b, 'character')
  expect_equal(b, "Bos taurus")
})

## ---------------------------------------------- ##
#           Testing species.kegg                   #
## ---------------------------------------------- ##
test_that("species.kegg works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- species.kegg(organism = 'rat', from = 'vulgar')
  b <- species.kegg(organism = 'Rat', from = 'vulgar')
  c <- species.kegg(organism = 'yeast', from = 'vulgar')
  d <- species.kegg('chimpanzee', from = 'vulgar')
  e <- species.kegg(organism = 'Homo sapiens', from = 'vulgar')
  f <- species.kegg(organism = 'Homo sapiens')
  g <- species.kegg(organism = 'Escherichia coli (strain K12)')
  h <- species.kegg('poppycok')
  i <- species.kegg('cfa', from = '3-letter')
  j <- species.kegg('www', from = '3-letter')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 4)
  expect_equal(ncol(a), 2)
  expect_equal(a$organism[1], 'rno')

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 4)
  expect_equal(ncol(b), 2)
  expect_equal(b$organism[1], 'rno')

  expect_is(c, 'data.frame')
  expect_equal(nrow(c), 2)
  expect_equal(ncol(c), 2)
  expect_equal(c$organism[1], 'sce')

  expect_is(d, 'data.frame')
  expect_equal(nrow(d), 1)
  expect_equal(ncol(d), 2)
  expect_equal(d$organism[1], 'ptr')

  expect_is(e, 'character')
  expect_is(f, 'data.frame')
  expect_is(g, 'data.frame')
  expect_equal(nrow(g), 65)
  expect_is(h, 'character')
  expect_is(i, 'data.frame')
  expect_equal(nrow(i), 1)
  expect_is(j, 'character')
})
