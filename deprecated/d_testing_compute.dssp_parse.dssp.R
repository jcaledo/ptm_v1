## ---------------------------------------------- ##
#             Testing parse.dssp                   #
## ---------------------------------------------- ##
# test_that("parse.dssp() works properly", {
#
#   skip_on_cran()
#   skip_on_travis()
#
#   compute.dssp('3cwm')
#   a <- parse.dssp('3cwm.dssp')
#
#   expect_equal(nrow(a), 369)
#   expect_equal(ncol(a), 8)
#   expect_equal(sum(c('B', 'C', 'E', 'G', 'H', 'S', 'T') %in% a$ss), 7)
# })



## ---------------------------------------------- ##
#            Testing compute.dssp                  #
## ---------------------------------------------- ##
# test_that("compute.dssp() works properly", {
#
#   skip_on_cran()
#   skip_on_travis()
#
#   compute.dssp(pdb = "./pdb/1u8f.pdb", destfile = "./pdb/")
#   if (file.exists('./pdb/1u8f.dssp')){
#     b <- parse.dssp("./pdb/1u8f.dssp")
#     if (!is.null(b)){
#       expect_is(b, 'data.frame')
#       expect_equal(nrow(b), 1332)
#       expect_equal(ncol(b), 8)
#     }
#   }
# })

