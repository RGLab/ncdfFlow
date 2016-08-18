#' return nothing when pass the test
#' @param orig the original flowFrame served as the reference 
#' @param new the new flowFrame to test 
is_equal_flowFrame <- function(orig, new, exprs = TRUE, description = TRUE){
  
  if(exprs){
    expect_equal(orig@exprs, new@exprs, tol = 1e-07, check.attributes = FALSE)
  }
  
  
  expect_equal(orig@parameters, new@parameters, tol = 1e-07, check.attributes = FALSE)
  
  #keyword may have minor change
  if(description)
    expect_true(all(is.element(orig@description,new@description)))
}

#' @param fs_orig the original flowSet served as the reference 
#' @param fs_new the new flowSet to test 
is_equal_flowSet <- function(fs_orig, fs_new, ...){
  all.equal(sampleNames(fs_orig), sampleNames(fs_new))
  invisible(lapply(sampleNames(fs_orig), function(sn){
    orig <- fs_orig[[sn]]
    target <- fs_new[[sn]]
    is_equal_flowFrame(orig, target, ...)
  }))  
} 
