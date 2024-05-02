library(ptm)
context("Sequence Recovering Functions Tests")

## ---------------------------------------------- ##
#               Testing get.seq                    #
## ---------------------------------------------- ##
test_that("get.seq() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- get.seq('P01009')
  b <- get.seq('P01009', as.string = FALSE)
  c <- get.seq('P010091')

  if(!is.null(a)){
    expect_is(a, 'character')
    expect_equal(nchar(a), 418)
  }

  if (!is.null(b)){
    expect_is(b, 'list')
    expect_equal(length(b[[1]]), 418)
    expect_equal(b[[1]][28], 'Q')
  }

  expect_is(c, 'NULL')
})

test_that("get.seq() works properly with MetOSite",{

  skip_on_cran()
  skip_on_travis()

  a <- get.seq('P01009', db = 'metosite')
  t <- get.seq('P01009', db = 'metosite', as.string = FALSE)

  if (!is.null(a)){
    expect_is(a, 'character')
    expect_equal(nchar(a), 394)
  }
  if (!is.null(t)){
    expect_is(t, 'list')
    expect_equal(length(t[[1]]), 394)
    expect_equal(t[[1]][28], 'P')
  }

})

test_that("get.seq() works properly with PDB",{

  skip_on_cran()
  skip_on_travis()

  a <- get.seq('1u8f', db = 'pdb')
  b <- get.seq('1u8f:O', db = 'pdb')
  c <- get.seq('1u8f:P', db = 'pdb')
  d <- get.seq('1u8f:R', db = 'pdb')
  e <- get.seq('1u8f:Q', db = 'pdb', as.string = FALSE)

  if (!is.null(a)){
    expect_is(a, "character")
    expect_equal(nchar(a), 335)
  }

  if (!is.null(b)){
    expect_is(b, "character")
    expect_equal(nchar(b), 335)
  }

  if (!is.null(c)){
    expect_is(c, "character")
    expect_equal(nchar(c), 335)
  }

  if (!is.null(d)){
    expect_is(d, "character")
    expect_equal(nchar(d), 335)
  }

  if (!is.null(e)){
    expect_is(e, "list")
    expect_is(e[[1]], 'character')
    expect_equal(length(e[[1]]), 335)
  }

})


test_that('get.seq() works properly with KEGG', {

  skip_on_cran()
  skip_on_travis()

  a <- get.seq('hsa:5265', 'kegg-aa')
  b <- get.seq('hsa:5265', db = 'kegg-aa', as.string = FALSE)
  c <- get.seq('hsa:5265', 'kegg-nt')
  d <- get.seq('hsa:5265', db = 'kegg-nt', as.string = FALSE)
  e <- get.seq('eih:ECOK1_2344', 'kegg-aa')

  if (!is.null(a)){
    expect_is(a, 'character')
    expect_equivalent(nchar(a), 418)
  }

  if (!is.null(b)){
    expect_is(b, 'list')
    expect_is(b[[1]], 'character')
    expect_equal(length(b[[1]]), 418)
  }

  if (!is.null(c)){
    expect_is(c, 'character')
    expect_equivalent(nchar(c), 1257)
  }

  if (!is.null(d)){
    expect_is(d, 'list')
    expect_is(d[[1]], 'character')
    expect_equal(length(d[[1]]), 1257)
  }

  if (!is.null(e)){
    expect_is(e, 'character')
  }
})


## ---------------------------------------------- ##
#               Testing  pdb.seq                   #
## ---------------------------------------------- ##
test_that('pdb.quaternary works properly', {

  skip_on_cran()
  skip_on_travis()

  a <- pdb.seq('1bpl')
  if (!is.null(a)){
    expect_is(a, 'data.frame')
    expect_equal(nrow(a), 2)
    expect_equal(ncol(a), 6)
  }
})

## ---------------------------------------------- ##
#            Testing  pdb.quaternary               #
## ---------------------------------------------- ##
test_that('pdb.quaternary works properly', {

  skip_on_cran()
  skip_on_travis()

  a <- pdb.quaternary("./pdb/1u8f.pdb")
  b <- pdb.quaternary("2hhb", keepfiles = TRUE)

  if (!is.null(a)){
    expect_is(a, 'list')
    expect_equal(length(a), 4)
    expect_true(is.matrix(a[[1]]))
    expect_equal(a[[1]][1], 0)
    expect_equal(nchar(a[[2]][1]), 335)
    expect_equal(a[[3]], c("O", "P", "Q", "R"))
    expect_equal(a[[4]], '1u8f')
  }

  if (!is.null(b)){
    expect_is(b, 'list')
    expect_equal(length(b), 4)
    expect_true(is.matrix(b[[1]]))
    expect_equal(b[[1]][1], 0)
    expect_equal(nchar(b[[2]][1]), 141)
    expect_equal(b[[3]], c("A", "B", "C", "D"))
    expect_equal(b[[4]], "2hhb")
    expect_true(file.exists("./2hhb.fa"))
  }
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

  if (!is.null(a)){
    expect_equal(length(a), 12)
    expect_equal(a, c("A","B","C","D","E","F","G","H","I","J","K","L"))
  }
  if (!is.null(b)){
    expect_equal(length(b), 4)
    expect_true(file.exists('./split_chain/1u8f_P.pdb'))
  }
  unlink("./split_chain", recursive = TRUE)
})

## ---------------------------------------------- ##
#            Testing pdb2uniprot                   #
## ---------------------------------------------- ##

test_that('pdb2uniprot works properly', {

  skip_on_cran()
  skip_on_travis()

  a <- pdb2uniprot('2OCC', 'M')
  x <- pdb2uniprot('xxxx', 'X')

  if (!is.null(a)){
    expect_is(a, 'character')
    expect_match(a, 'P10175')
  }

  expect_is(x, 'NULL')
})


## ---------------------------------------------- ##
#            Testing uniprot2pdb                   #
## ---------------------------------------------- ##

test_that('uniprot2pdb works properly',{

  skip_on_cran()
  skip_on_travis()

  b <- uniprot2pdb('P10175')
  x <- uniprot2pdb('P010091')

  if (!is.null(b)){
    expect_is(b, 'data.frame')
    expect_equal(nrow(b), 100)
    expect_equal(ncol(b), 7)
  }

  expect_is(x, 'NULL')
})


## ---------------------------------------------- ##
#             Testing prot2codon                   #
## ---------------------------------------------- ##
test_that("prot2codon() works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- prot2codon("P01009")
  b <- prot2codon('1ATU', 'A')
  c <- prot2codon("./pdb/1u8f.pdb", chain = 'O')
  d <- prot2codon('1A00', chain = 'A')
  # e <- prot2codon('1F8A', 'B')
  e <- prot2codon('1cll', 'A')

  if (!is.null(a)){
    expect_is(a, 'data.frame')
    expect_equal(sum(table(a$check)), nrow(a))
    expect_equal(ncol(a), 6)
  }

  if (!is.null(b)){
    expect_is(b, 'data.frame')
    expect_equal(sum(table(b$check)), nrow(b))
    expect_equal(ncol(b), 6)
  }

  if (!is.null(c)){
    expect_is(c, 'data.frame')
    expect_equal(sum(table(c$check)), nrow(c))
    expect_equal(ncol(c), 6)
  }

  if (!is.null(d)){
    expect_is(d, 'data.frame')
    expect_equal(nrow(d), 141)
    expect_equal(ncol(d), 6)
    expect_equal(sum(table(d$check)), nrow(d))
  }

  if (!is.null(e)){
    expect_is(e, 'data.frame')
    expect_equal(nrow(e), 148)
    expect_equal(ncol(e), 6)
    expect_equal(sum(table(e$check)), nrow(e))

    # expect_is(e, 'data.frame')
    # expect_equal(nrow(e), 167)
    # expect_equal(ncol(e), 6)
    # expect_equal(sum(table(e$check)), nrow(e) - 4)
  }

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

  if (!is.null(a)){
    expect_is(a, 'list')
    expect_equal(length(a), 2)
    expect_equal(a[[1]], '1u8f')
    expect_equal(a[[2]], 'O')
    expect_gt(attributes(a)$coverage, 0.9)
  }

  if (!is.null(b)){
    expect_is(b, 'list')
    expect_equal(length(b), 2)
    expect_equal(b[[1]], '5jqa')
    expect_equal(b[[2]], 'A')
    expect_gt(attributes(b)$coverage, 0.9)
  }

  if (!is.null(c)){
    expect_is(c, 'list')
    expect_equal(length(c), 2)
    expect_equal(c[[1]], '4wij')
    expect_equal(c[[2]], 'B')
    expect_lt(attributes(c)$coverage, 0.43)
  }

  expect_is(d, 'NULL')
})


## ---------------------------------------------- ##
#               Testing id.mapping                 #
## ---------------------------------------------- ##
test_that('id.mapping() works properly when pdb -> uniprot', {

  skip_on_cran()
  skip_on_travis()

  a <- id.mapping('3cwm', 'pdb','uniprot')
  b <- id.mapping('2occ', 'pdb', 'uniprot')

  if (!is.null(a)){
    expect_is(a, 'character')
    expect_true(a == 'P01009')
  }

  if (!is.null(b)){
    expect_is(b, 'character')
    expect_equal(length(b), 13)
    expect_true('P00415' %in% b)
  }

})


test_that('id.mapping() works properly when up -> pdb', {

  skip_on_cran()
  skip_on_travis()

  a <- id.mapping('P01009', 'uniprot', 'pdb')
  b <- id.mapping('P00415', 'uniprot', 'pdb')

  if (!is.null(a)){
    expect_is(a, 'character')
    expect_true('3CWM' %in% a)
  }

  if (!is.null(b)){
    expect_is(b, 'character')
    expect_true('2OCC' %in% b)
  }

})


test_that('id.mapping() works properly when up -> kegg', {

  skip_on_cran()
  skip_on_travis()

  up_kegg <- id.mapping('P01009', 'uniprot', 'kegg')

  if (!is.null(up_kegg)){
    expect_is(up_kegg, 'character')
    expect_equivalent(up_kegg, 'hsa:5265')
  }

})


test_that('id.mapping() works properly when kegg -> up', {

  skip_on_cran()
  skip_on_travis()

  kegg_up <- id.mapping('hsa:5265', 'kegg', 'uniprot')

  if (!is.null(kegg_up)){
    expect_is(kegg_up, 'character')
    expect_true("P01009" %in% kegg_up)
  }

})


test_that('id.mapping() works properly when pdb -> kegg', {

  skip_on_cran()
  skip_on_travis()

  a <- id.mapping('3cwm', 'pdb','kegg')
  # b <- id.mapping('2occ', 'pdb', 'kegg')

  if (!is.null(a)){
    expect_is(a, 'list')
    expect_is(a[[1]], 'character')
    expect_equivalent(a[[1]], "hsa:5265")
  }

  # if (!is.null(b)){
  #   expect_is(b, 'list')
  #   expect_is(b[[3]], 'character')
  #   expect_true("bta:281090" %in% b[[3]])
  # }

})


test_that('id.mapping() works properly when kegg <- pdb', {

  skip_on_cran()
  skip_on_travis()

  a <- id.mapping('hsa:5265', 'kegg','pdb')
  b <- id.mapping('bta:281090', 'kegg','pdb')

  if (!is.null(a)){
    expect_is(a, 'character')
    expect_true('3CWM' %in% a)
  }

  if (!is.null(b)){
    expect_is(b, 'character')
    expect_true('2OCC' %in% b)
  }

})


## ---------------------------------------------- ##
#            Testing id.features                   #
## ---------------------------------------------- ##
test_that("id.features() works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- id.features('P01009')
  b <- id.features('P01009', features = 'ec,keyword,xref_pdb')

  if (!is.null(a)){
    expect_is(a, 'data.frame')
    expect_equal(length(a), 5)
    expect_equal(a$Reviewed, "reviewed")
  }

  if (!is.null(b)){
    expect_is(b, 'data.frame')
    expect_equal(length(b), 8)
  }

  if (!is.null(a) & !is.null(b)){
    expect_equal(a$Organism, b$Organism)
  }

})


## ---------------------------------------------- ##
#           Testing species.mapping                #
## ---------------------------------------------- ##
test_that("species.mapping works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- species.mapping('P01009')
  b <- species.mapping('2OCC', db = 'pdb')
  c <- species.mapping('P010091')

  if (!is.null(a)){
    expect_is(a, 'character')
    expect_equal(a, "Homo sapiens")
  }

  if (!is.null(b)){
    expect_is(b, 'character')
    expect_equal(b, "Bos taurus")
  }

  expect_is(c, "NULL")
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

  if (!is.null(a)){
    expect_is(a, 'data.frame')
    expect_equal(nrow(a), 4)
    expect_equal(ncol(a), 2)
    expect_equal(a$organism[1], 'rno')
  }

  if (!is.null(b)){
    expect_is(b, 'data.frame')
    expect_equal(nrow(b), 4)
    expect_equal(ncol(b), 2)
    expect_equal(b$organism[1], 'rno')
  }

  if (!is.null(c)){
    expect_is(c, 'data.frame')
    expect_equal(nrow(c), 2)
    expect_equal(ncol(c), 2)
    expect_equal(c$organism[1], 'sce')
  }

  if (!is.null(d)){
    expect_is(d, 'data.frame')
    expect_equal(nrow(d), 1)
    expect_equal(ncol(d), 2)
    expect_equal(d$organism[1], 'ptr')
  }

  expect_is(e, 'NULL')

  if (!is.null(f)){
    expect_is(f, 'data.frame')
  }

  if (!is.null(g)){
    expect_is(g, 'data.frame')
    expect_equal(nrow(g), 65)
  }

  expect_is(h, 'NULL')

  if (!is.null(i)){
    expect_is(i, 'data.frame')
    expect_equal(nrow(i), 1)
  }

  expect_is(j, 'NULL')
})
