library(ptm)
context("Abundance DB")

## ----------------------------------------------- ##
#     Testing the function abundance                #
## ----------------------------------------------- ##
test_that("abundance() works properly", {

  a <- abundance(id = 'A0AVT1')
  b <- abundance(id = 'A0AVT1', 'jurkat')
  c <- abundance(id = 'A0AVT1', 'hela')
  d <- abundance(id = 'P30034') # from Sus scrofa

  expect_is(a, 'numeric')
  expect_is(b, 'numeric')
  expect_is(c, 'numeric')
  expect_is(d, 'numeric')

  expect_equal(attributes(a)$units, "ppm")
  expect_equal(attributes(b)$units, "ppm")
  expect_equal(attributes(c)$units, "ppm")
  expect_equal(attributes(d)$units, "ppm")

  expect_equal(attributes(a)$species, "Homo sapiens")
  expect_equal(attributes(b)$species, "Homo sapiens")
  expect_equal(attributes(c)$species, "Homo sapiens")
  expect_equal(attributes(d)$species, "Sus scrofa")

  expect_equal(attributes(a)$string, "A0AVT1")
  expect_equal(attributes(b)$string, "A0AVT1")
  expect_equal(attributes(c)$string, "A0AVT1")
  expect_equal(attributes(d)$string, "P30034")

  expect_equivalent(a, 63.7)
  expect_equivalent(b, 27.4)
  expect_equivalent(c, 29.3)
  expect_equivalent(d, 26.9)

})


