context("misc")

test_that("bit vector", {

  #example logical index
  indx <- rep(F, 1e6)
  indx[sample.int(1e6, 1e5)] <- T
  
  #convert to bit vector
  bitVec <- toBitVec(indx)
  expect_is(bitVec, "raw")
  expect_equal(attr(bitVec, "bitlen"), 1e6)
  #convert it back
  indx1 <- toLogical(bitVec)
  expect_equal(indx, indx1)
})