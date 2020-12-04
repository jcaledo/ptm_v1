library(ptm)
context("Distances")

## ---------------------------------------------- ##
#            Testing pairwise.dist                 #
## ---------------------------------------------- ##
test_that('pairwise.dist() works properly',{

  a <- matrix(runif(30*3), ncol=3)
  b <- matrix(runif(10*3), ncol=3)
  M <- pairwise.dist(a, b, square = FALSE)
  a <- matrix(c(1,1,1,1,0,0,1,0,1), ncol = 3, byrow = TRUE)
  b <- matrix(c(0,0,0,1,1,1), ncol = 3, byrow = TRUE)
  D <- pairwise.dist(a, b, square = FALSE)

  expect_equal(dim(M), c(30,10))
  expect_equal(dim(D), c(3,2))
  expect_equivalent(D, matrix(c(sqrt(3), 0, 1, sqrt(2), sqrt(2), 1), ncol = 2, byrow = TRUE))
})

## ---------------------------------------------- ##
#            Testing res.dist                      #
## ---------------------------------------------- ##
test_that('res.dist() works properly', {

  skip_on_cran()
  skip_on_travis()

  a <- res.dist(pdb = '1q8k', 51, 'A', 55, 'A', backbone = TRUE, hatoms = TRUE)
  b <- res.dist(pdb = '1q8k', 51, 'A', 55, 'A', backbone = FALSE, hatoms = FALSE)
  c <- res.dist(pdb = '2occ', 51, 'A', 51, 'N')
  d <- res.dist(pdb = '1q8k', 43, 'A', 65, 'A') # two G residues)
  e <- res.dist(pdb = '1q8k', 94, 'A', 187, 'A') # two A residues


  expect_equal(length(a), 3)
  expect_equal(length(a[[1]]), 4)
  expect_equal(length(a[[2]]), 4)
  expect_equal(length(a[[3]]), 4)
  expect_equal(a[[1]][[1]], 2.96)
  expect_true(is.character(a[[1]][[2]]))
  expect_true(is.character(a[[1]][[3]]))
  expect_equal(a[[2]][[1]], 9.28)
  expect_true(is.character(a[[2]][[2]]))
  expect_true(is.character(a[[2]][[3]]))
  expect_equal(a[[3]][[1]], 5.96)
  expect_true(is.character(a[[3]][[2]]))
  expect_true(is.character(a[[3]][[3]]))
  expect_lt(a[[1]][[1]], b[[1]][[1]])
  expect_gt(a[[2]][[1]], b[[2]][[1]])
  expect_gt(a[[3]][[1]], b[[3]][[1]])

  expect_equal(length(c), 3)
  expect_equal(length(c[[1]]), 4)
  expect_equal(length(c[[2]]), 4)
  expect_equal(length(c[[3]]), 4)
  expect_equal(c[[1]][[1]], 89.11)
  expect_true(is.character(c[[1]][[2]]))
  expect_true(is.character(c[[1]][[3]]))
  expect_equal(c[[2]][[1]], 92.68)
  expect_true(is.character(c[[2]][[2]]))
  expect_true(is.character(d[[2]][[3]]))
  expect_equal(c[[3]][[1]], 90.5)
  expect_true(is.character(c[[3]][[2]]))
  expect_true(is.character(c[[3]][[3]]))

  expect_equal(d[[1]][[1]], 15.23)
  expect_equal(d[[2]][[1]], 15.23)
  expect_equal(d[[3]][[1]], 15.23)
  expect_equal(e[[1]][[1]], e[[2]][[1]])
  expect_equal(e[[1]][[1]], e[[3]][[1]])
})


## ---------------------------------------------- ##
#              Testing dist2closest                #
## ---------------------------------------------- ##
test_that('dist2closest() works properly', {

  skip_on_cran()
  skip_on_travis()

  a <- dist2closest(pdb = '1a1x', res = 103, chain = 'A', aa = 'E' )
  b <- dist2closest('1GZX', 101, chain = 'A')

  expect_is(a, 'numeric')
  expect_equal(a[1], 6.72)
  expect_equal(attr(a, "closest residue"), "CG-GLU-37-A <--> CE-LYS-103-A")

  expect_is(b, 'numeric')
  expect_equal(b[1], 3.66)
  expect_equal(attr(b, "closest residue"), "SD-MET-32-A <--> CB-LEU-101-A")
})


## ---------------------------------------------- ##
#                 Testing ball                     #
## ---------------------------------------------- ##
test_that('ball() works properly', {

  skip_on_cran()
  skip_on_travis()

  M <- ball(pdb = '6e7f', chain = 'A', res = 181, r = 6)
  Mb <- ball(pdb = '6e7f', chain = 'A', res = 181, r = 6, backbone = TRUE)

  expect_is(M, 'data.frame')
  expect_equal(dim(M), c(23, 8))
  expect_gt(min(M$distance), 3.5)
  expect_lt(max(M$distance), 6)
  expect_is(Mb, 'data.frame')
  expect_equal(dim(Mb), c(41, 8))
  expect_gt(min(Mb$distance), 3.5)
  expect_lt(max(Mb$distance), 6)

  C <- ball(pdb = '6e7f', chain = 'A', res = 98, r = 6)
  Cb <- ball(pdb = '6e7f', chain = 'A', res = 98, r = 6, backbone = TRUE)

  expect_is(C, 'data.frame')
  expect_equal(dim(C), c(16, 8))
  expect_gt(min(C$distance), 3.5)
  expect_lt(max(C$distance), 6)
  expect_is(Cb, 'data.frame')
  expect_equal(dim(Cb), c(43, 8))
  expect_gt(min(Cb$distance), 3.5)
  expect_lt(max(Cb$distance), 6)

  Thr <- ball(pdb = '6e7f', chain = 'A', res = 72, r = 6)
  Thrb <- ball(pdb = '6e7f', chain = 'A', res = 72, r = 6, backbone = TRUE)

  expect_is(Thr, 'data.frame')
  expect_equal(dim(Thr), c(12, 8))
  expect_gt(min(Thr$distance), 3.5)
  expect_lt(max(Thr$distance), 6)
  expect_is(Thrb, 'data.frame')
  expect_equal(dim(Thrb), c(28, 8))
  expect_gt(min(Thrb$distance), 3.5)
  expect_lt(max(Thrb$distance), 6)

  S <- ball(pdb = '6e7f', chain = 'A', res = 118, r = 6)
  Sb <- ball(pdb = '6e7f', chain = 'A', res = 118, r = 6, backbone = TRUE)

  expect_is(S, 'data.frame')
  expect_equal(dim(S), c(9, 8))
  expect_gt(min(S$distance), 3)
  expect_lt(max(S$distance), 6)
  expect_is(Sb, 'data.frame')
  expect_equal(dim(Sb), c(35, 8))
  expect_gt(min(Sb$distance), 2)
  expect_lt(max(Sb$distance), 6)

  D <- ball(pdb = '6e7f', chain = 'A', res = 100, r = 6)
  Db <- ball(pdb = '6e7f', chain = 'A', res = 100, r = 6, backbone = TRUE)

  expect_is(D, 'data.frame')
  expect_equal(dim(D), c(12, 8))
  expect_gt(min(D$distance), 3.5)
  expect_lt(max(D$distance), 6)
  expect_is(Db, 'data.frame')
  expect_equal(dim(Db), c(26, 8))
  expect_gt(min(Db$distance), 3.5)
  expect_lt(max(Db$distance), 6)


  E <- ball(pdb = '6e7f', chain = 'A', res = 105, r = 6)
  Eb <- ball(pdb = '6e7f', chain = 'A', res = 105, r = 6, backbone = TRUE)

  expect_is(E, 'data.frame')
  expect_equal(dim(E), c(12, 8))
  expect_gt(min(E$distance), 3.3)
  expect_lt(max(E$distance), 6)
  expect_is(Eb, 'data.frame')
  expect_equal(dim(Eb), c(14, 8))
  expect_gt(min(Eb$distance), 3.3)
  expect_lt(max(Eb$distance), 6)

  H <- ball(pdb = '6e7f', chain = 'A', res = 135, r = 6)
  Hb <- ball(pdb = '6e7f', chain = 'A', res = 135, r = 6, backbone = TRUE)

  expect_is(H, 'data.frame')
  expect_equal(dim(H), c(4, 8))
  expect_gt(min(H$distance), 3.5)
  expect_lt(max(H$distance), 6)
  expect_is(Hb, 'data.frame')
  expect_equal(dim(Hb), c(11, 8))
  expect_gt(min(Hb$distance), 3.5)
  expect_lt(max(Hb$distance), 6)

  R <- ball(pdb = '6e7f', chain = 'A', res = 182, r = 6)
  Rb <- ball(pdb = '6e7f', chain = 'A', res = 182, r = 6, backbone = TRUE)

  expect_is(R, 'data.frame')
  expect_equal(dim(R), c(10, 8))
  expect_gt(min(R$distance), 3.5)
  expect_lt(max(R$distance), 6)
  expect_is(Rb, 'data.frame')
  expect_equal(dim(Rb), c(17, 8))
  expect_gt(min(Rb$distance), 3.5)
  expect_lt(max(Rb$distance), 6)

  K <- ball(pdb = '6e7f', chain = 'A', res = 204, r = 6)
  Kb <- ball(pdb = '6e7f', chain = 'A', res = 204, r = 6, backbone = TRUE)

  expect_is(K, 'data.frame')
  expect_equal(dim(K), c(6, 8))
  expect_gt(min(K$distance), 3.5)
  expect_lt(max(K$distance), 6)
  expect_is(Kb, 'data.frame')
  expect_equal(dim(Kb), c(16, 8))
  expect_gt(min(Kb$distance), 2.7)
  expect_lt(max(Kb$distance), 6)

  Phe <- ball(pdb = '6e7f', chain = 'A', res = 214, r = 6)
  Pheb <- ball(pdb = '6e7f', chain = 'A', res = 214, r = 6, backbone = TRUE)

  expect_is(Phe, 'data.frame')
  expect_equal(dim(Phe), c(23, 8))
  expect_gt(min(Phe$distance), 3.5)
  expect_lt(max(Phe$distance), 6)
  expect_is(Pheb, 'data.frame')
  expect_equal(dim(Pheb), c(39, 8))
  expect_gt(min(Pheb$distance), 3.0)
  expect_lt(max(Pheb$distance), 6)

  Y <- ball(pdb = '6e7f', chain = 'A', res = 219, r = 6)
  Yb <- ball(pdb = '6e7f', chain = 'A', res = 219, r = 6, backbone = TRUE)

  expect_is(Y, 'data.frame')
  expect_equal(dim(Y), c(5, 8))
  expect_gt(min(Y$distance), 1.5)
  expect_lt(max(Y$distance), 6)
  expect_is(Yb, 'data.frame')
  expect_equal(dim(Yb), c(19, 8))
  expect_gt(min(Yb$distance), 1.5)
  expect_lt(max(Yb$distance), 6)

  W <- ball(pdb = '6e7f', chain = 'A', res = 117, r = 6)
  Wb <- ball(pdb = '6e7f', chain = 'A', res = 117, r = 6, backbone = TRUE)

  expect_is(W, 'data.frame')
  expect_equal(dim(W), c(16, 9))
  expect_gt(min(W$distance1), 1.5)
  expect_gt(max(W$distance1), 7)
  expect_is(Wb, 'data.frame')
  expect_equal(dim(Wb), c(25, 9))
  expect_gt(min(Wb$distance2), 1.5)
  expect_gt(max(Wb$distance2), 7)
})
