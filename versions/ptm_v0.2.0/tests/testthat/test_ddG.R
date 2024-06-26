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

  a <- foldx.mut(pdb = './pdb/1b4i_Repair.pdb', ch = 'A', pos = 8, method = 'buildmodel')
  b <- foldx.mut(pdb = "./pdb/1b4i_Repair.pdb", ch = 'A', pos = 8, method = 'positionscan')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 19)
  expect_equal(ncol(a), 8)
  expect_equal(a$WT[1], "S")
  expect_equal(a$Position[1], 8)
  expect_lt(a$DDG[7], 0) # S -> E

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 24)
  expect_equal(ncol(b), 8)
  expect_equal(b$WT[1], "S")
  expect_equal(b$Position[1], 8)
  expect_gt(b$DDG[23], 0) # S -> pS

  # a <- foldx(pdb = '1pga', ch = 'A', pos = 31, method = 'buildmodel')
  # b <- foldx(pdb = '1pga_Repair.pdb', ch = 'A', pos = 31, method = 'positionscan', wfile = "./pdb")
  #
  # expect_is(a, 'data.frame')
  # expect_equal(nrow(a), 19)
  # expect_equal(ncol(a), 8)
  # expect_equal(a$WT[1], "K")
  # expect_equal(a$Position[1], 31)
  # expect_lt(a$DDG[11], 0) # K -> L
  #
  # expect_is(b, 'data.frame')
  # expect_equal(nrow(b), 24)
  # expect_equal(ncol(b), 8)
  # expect_equal(b$WT[1], "K")
  # expect_lt(b$DDG[3], 0) # K -> L

  # a <- foldx(pdb = './pdb/1u8f_Repair.pdb', ch = 'O', pos = 46, method = "buildmodel", wfile = "./pdb")
  # b <- foldx(pdb = './pdb/1u8f_Repair.pdb', ch = 'O', pos = 46, method = "positionscan", wfile = "./pdb")

  # expect_is(a, 'data.frame')
  # expect_equal(nrow(a), 19)
  # expect_equal(ncol(a), 8)
  # expect_equal(a$WT[1], "M")
  # expect_equal(a$Position[1], 46)
  # expect_lt(a$DDG[11], 0) # M -> L

  # expect_is(b, 'data.frame')
  # expect_equal(nrow(b), 24)
  # expect_equal(ncol(b), 8)
  # expect_equal(b$WT[1], "M")
  # expect_lt(b$DDG[3], 0) # M -> L

})


## ----------------------------------------------- ##
#         Testing the function foldx.stab           #
## ----------------------------------------------- ##
test_that("foldx.stab() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- foldx.stab(pdb = "./pdb/2dfd_Repair.pdb")

  expect_is(a, 'numeric')
  expect_lt(a, -105)
})


## ----------------------------------------------- ##
#       Testing the function foldx.assembly         #
## ----------------------------------------------- ##
test_that("foldx.assembly() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- foldx.assembly(pdb = './pdb/2dfd_Repair.pdb', mol1 = 'A', mol2 = 'B')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 90)
  expect_equal(ncol(a), 4)
  expect_equal(a$id[13], "LA50")

})

## ----------------------------------------------- ##
#         Testing the function ddG.profile          #
## ----------------------------------------------- ##
test_that("ddG.profile() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- ddG.profile(prot = "1b4i", ch = "A", pos = 8)

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 19)
  expect_equal(ncol(a), 4)

  # a <- ddG.profile(prot = "./pdb/1pga_Reapir.pdb", ch = "A")
  #
  # expect_is(a, 'data.frame')
  # expect_equal(nrow(a), 19)
  # expect_equal(ncol(a), 4)
})


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


  # a <- ddG.ptm(pdb = "1fin", ch = "A", pos = 160, ptm = "pThr", dir = 'f', pH = 7)
  # b <- ddG.ptm(pdb = "./pdb/1fin_Repair.pdb", ch = "A", pos = 160, ptm = "pThr", dir = 'f', pH = 7)
  # c <- ddG.ptm(pdb = "./pdb/1fin_Repair.pdb", ch = "A", pos = 160, ptm = "pSer", dir = 'f', pH = 7)

  expect_is(a, 'character')
  expect_gt(as.numeric(a), 0)
  expect_equal(attributes(a)$units, "kcal/mol")
  expect_equal(attributes(a)$PTM, "pSer")
  expect_equal(attributes(a)$position, 8)
  expect_equal(attributes(a)$`wild-type`, "SER")

  expect_is(b, 'character')
  expect_gt(as.numeric(b), 0)
  expect_equal(attributes(b)$units, "kcal/mol")
  expect_equal(attributes(b)$PTM, "pSer")
  expect_equal(attributes(b)$position, 15)
  expect_equal(attributes(b)$`wild-type`, "SEP")

  expect_is(c, 'character')
  expect_gt(as.numeric(c), 0)
  expect_equal(attributes(c)$units, "kcal/mol")
  expect_equal(attributes(c)$PTM, "pThr")
  expect_equal(attributes(c)$position, 23)
  expect_equal(attributes(c)$`wild-type`, "THR")

  expect_is(d, 'character')
  expect_lt(as.numeric(d), 0)
  expect_equal(attributes(d)$units, "kcal/mol")
  expect_equal(attributes(d)$PTM, "MetO-Q")
  expect_equal(attributes(d)$position, 29)
  expect_equal(attributes(d)$`wild-type`, "MET")

  expect_is(e, 'character')
  expect_lt(as.numeric(e), 0)
  expect_equal(attributes(e)$units, "kcal/mol")
  expect_equal(attributes(e)$PTM, "MetO-T")
  expect_equal(attributes(e)$position, 29)
  expect_equal(attributes(e)$`wild-type`, "MET")

  expect_warning(ddG.ptm(pdb = "./pdb/1b4i_Repair.pdb", ch = 'A', pos = 21, ptm = 'pThr', dir = 'f'),
                 "Wild type residue is not Thr")
})
if (file.exists('1b4i_MetO-Q29.pdb')){
  file.remove('1b4i_MetO-Q29.pdb')
}
if (file.exists('1b4i_MetO-T29.pdb')){
  file.remove('1b4i_MetO-T29.pdb')
}
if (file.exists('1b4i_pSer8.pdb')){
  file.remove('1b4i_pSer8.pdb')
}
if (file.exists('1b4i_pSer15.pdb')){
  file.remove('1b4i_pSer15.pdb')
}
if (file.exists('1b4i_pThr21.pdb')){
  file.remove('1b4i_pThr23.pdb')
}
