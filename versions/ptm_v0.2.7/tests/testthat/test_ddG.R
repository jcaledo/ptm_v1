library(ptm)
context("Compute Changes in Stability (DDG)")

## ----------------------------------------------- ##
#         Testing the function imutant              #
## ----------------------------------------------- ##
test_that("imutant() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- imutant(protein = '1b4i', pos = 8, newres = 'E')
  b <- imutant(protein = "MISSVCVSSYRGPKSGNKPPSKTCLKEEMA", pos = 8, newres = 'E')
  c <- imutant(protein = 'xxxx', pos = 8, newres = 'E')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 1)
  expect_equal(ncol(a), 7)
  expect_equal(a$WT, "S")
  expect_equal(a$Position, "8")
  expect_equal(a$DDG, 0.72)

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 1)
  expect_equal(ncol(b), 6)
  expect_equal(b$WT, "S")
  expect_equal(b$Position, "8")
  expect_gt(b$DDG, 0.7)

  expect_is(c, 'NULL')

  # a <- imutant(protein = '1pga', pos = 31, newres = 'D')
  # b <- imutant(protein = get.seq("P06654"), pos = 31 + 226)
  #
  # expect_is(a, 'data.frame')
  # expect_equal(nrow(a), 1)
  # expect_equal(ncol(a), 7)
  # expect_equal(a$WT, "K")
  # expect_equal(a$Position, "31")
  # expect_equal(a$DDG, 1.18)
  #
  # expect_is(b, 'data.frame')
  # expect_equal(nrow(b), 19)
  # expect_equal(ncol(b), 6)
  # expect_equal(b$WT[1], "K")
  # expect_equal(b$Position[1], "257") # 31 + 226
  # expect_equal(b$DDG[19], 1.56) # K -> D

})


## ----------------------------------------------- ##
#         Testing the function foldx.mut            #
## ----------------------------------------------- ##
test_that("foldx.mut() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- foldx.mut(pdb = './pdb/1b4i_Repair.pdb', ch = 'A', pos = 8, method = 'buildmodel', keepfiles = TRUE)
  b <- foldx.mut(pdb = "./pdb/1b4i_Repair.pdb", ch = 'A', pos = 8, method = 'positionscan', keepfiles = TRUE)
  c <- foldx.mut(pdb = "xxxx", ch = 'A', pos = 8, method = 'buildmodel')

  if (!is.null(a)){ # if FoldX is installed
    expect_is(a, 'data.frame')
    expect_equal(nrow(a), 19)
    expect_equal(ncol(a), 8)
    expect_equal(a$WT[1], "S")
    expect_equal(a$Position[1], 8)
    expect_lt(a$DDG[7], 0) # S -> E
  }

  if (!is.null(b)){
    expect_is(b, 'data.frame')
    expect_equal(nrow(b), 24)
    expect_equal(ncol(b), 8)
    expect_equal(b$WT[1], "S")
    expect_equal(b$Position[1], 8)
    expect_gt(b$DDG[23], 0) # S -> pS
  }

  expect_is(c, 'NULL')
})


## ----------------------------------------------- ##
#         Testing the function foldx.stab           #
## ----------------------------------------------- ##
test_that("foldx.stab() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- foldx.stab(pdb = "./pdb/1b4i_Repair.pdb")
  b <- foldx.stab(pdb = "xxxx")

  if (!is.null(a)){ # If FoldX is installed
    expect_is(a, 'numeric')
    expect_lt(a, 55.9)
  }
  expect_is(b, 'NULL')
})


## ----------------------------------------------- ##
#       Testing the function foldx.assembly         #
## ----------------------------------------------- ##
test_that("foldx.assembly() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- foldx.assembly(pdb = './pdb/2dfd_Repair.pdb', mol1 = 'A', mol2 = 'B')
  b <- foldx.assembly(pdb = "xxxx", mol1 = 'A', mol2 = 'B')

  if (!is.null(a)){ # If FoldX is installed
    expect_is(a, 'data.frame')
    expect_equal(nrow(a), 90)
    expect_equal(ncol(a), 4)
    expect_equal(a$id[13], "LA50")
  }
  expect_is(b, 'NULL')

})

## ----------------------------------------------- ##
#         Testing the function ddG.profile          #
## ----------------------------------------------- ##
test_that("ddG.profile() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- ddG.profile(prot = "1b4i", ch = "A", pos = 8)
  b <- ddG.profile(prot = "xxxx", ch = "A", pos = 8)

  if (!is.null(a)){
    expect_is(a, 'data.frame')
    expect_equal(nrow(a), 19)
    expect_equal(ncol(a), 4)
  }
  expect_is(b, 'NULL')
})
if (file.exists("Rplots.pdf")){
  file.remove("Rplots.pdf")
}


## ----------------------------------------------- ##
#         Testing the function ddG.ptm              #
## ----------------------------------------------- ##
test_that("ddG.ptm() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- ddG.ptm(pdb = "./pdb/1b4i_Repair.pdb", ch = 'A', pos = 8, ptm = 'pSer', dir = 'f', pH = 7)
  b <- ddG.ptm(pdb = "./pdb/1b4i_Repair.pdb", ch = 'A', pos = 15, ptm = 'pSer', dir = 'b', pH = 7)
  c <- ddG.ptm(pdb = "./pdb/1b4i_Repair.pdb", ch = 'A', pos = 23, ptm = 'pThr', dir = 'f', pH = 7)
  d <- ddG.ptm(pdb = "./pdb/1b4i_Repair.pdb", ch = 'A', pos = 29, ptm = 'MetO-Q', dir = 'f', pH = 7)
  e <- ddG.ptm(pdb = "./pdb/1b4i_Repair.pdb", ch = 'A', pos = 29, ptm = 'MetO-T', dir = 'f', pH = 7)
  f <- ddG.ptm(pdb = "xxxx", ch = 'A', pos = 29, ptm = 'MetO-T', dir = 'f', pH = 7)

  if (!is.null(a)){
    expect_is(a, 'character')
    expect_gt(as.numeric(a), 0)
    expect_equal(attributes(a)$units, "kcal/mol")
    expect_equal(attributes(a)$PTM, "pSer")
    expect_equal(attributes(a)$position, 8)
    expect_equal(attributes(a)$`wild-type`, "SER")
  }

  if (!is.null(b)){
    expect_is(b, 'character')
    expect_gt(as.numeric(b), 0)
    expect_equal(attributes(b)$units, "kcal/mol")
    expect_equal(attributes(b)$PTM, "pSer")
    expect_equal(attributes(b)$position, 15)
    expect_equal(attributes(b)$`wild-type`, "SEP")
  }

  if (!is.null(c)){
    expect_is(c, 'character')
    expect_gt(as.numeric(c), 0)
    expect_equal(attributes(c)$units, "kcal/mol")
    expect_equal(attributes(c)$PTM, "pThr")
    expect_equal(attributes(c)$position, 23)
    expect_equal(attributes(c)$`wild-type`, "THR")
  }

  if (!is.null(d)){
    expect_is(d, 'character')
    expect_lt(as.numeric(d), 0)
    expect_equal(attributes(d)$units, "kcal/mol")
    expect_equal(attributes(d)$PTM, "MetO-Q")
    expect_equal(attributes(d)$position, 29)
    expect_equal(attributes(d)$`wild-type`, "MET")
  }

  if (!is.null(e)){
    expect_is(e, 'character')
    expect_lt(as.numeric(e), 0)
    expect_equal(attributes(e)$units, "kcal/mol")
    expect_equal(attributes(e)$PTM, "MetO-T")
    expect_equal(attributes(e)$position, 29)
    expect_equal(attributes(e)$`wild-type`, "MET")
  }

  # expect_warning(ddG.ptm(pdb = "./pdb/1b4i_Repair.pdb", ch = 'A', pos = 21, ptm = 'pThr', dir = 'f'),
  #                "Wild type residue is not Thr")

  expect_is(f, 'NULL')


  if (file.exists('./pdb/1b4i_MetO-Q29.pdb')){
    file.remove('./pdb/1b4i_MetO-Q29.pdb')
  }
  if (file.exists('./pdb/1b4i_MetO-T29.pdb')){
    file.remove('./pdb/1b4i_MetO-T29.pdb')
  }
  if (file.exists('./pdb/1b4i_pSer8.pdb')){
    file.remove('./pdb/1b4i_pSer8.pdb')
  }
  if (file.exists('./pdb/1b4i_pSer15.pdb')){
    file.remove('./pdb/1b4i_pSer15.pdb')
  }
  if (file.exists('./pdb/1b4i_pThr23.pdb')){
    file.remove('./pdb/1b4i_pThr23.pdb')
  }
})
