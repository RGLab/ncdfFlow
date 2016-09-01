fs <- read.flowSet(path = system.file("extdata","compdata","data",package="flowCore"))
comp.fs <- ncdfFlowSet(fs)
test_that("compensate", {
      
      comp.mat <- as.matrix(read.table(system.file("extdata","compdata","compmatrix",package="flowCore"),header=TRUE,skip=2,check.names=FALSE))
            
      # Compensate with single comp
      row.names(comp.mat) <- colnames(comp.mat)
      expectRes <- fsApply(compensate(fs, comp.mat), colMeans, use.exprs = TRUE)
      
      comp.fs1 <- compensate(comp.fs, comp.mat)
      expect_equal(fsApply(comp.fs1, colMeans, use.exprs = TRUE), expectRes)
      
      # list
      comp <- sapply(sampleNames(comp.fs), function(sn)comp.mat, simplify = FALSE)
      comp.fs2 <- compensate(comp.fs, comp)
      
      expect_equal(fsApply(comp.fs2, colMeans, use.exprs = TRUE), expectRes)
      expect_failure(expect_equal(comp.fs@file, comp.fs2@file))
      
      # unmatched names
      names(comp)[1] <- "dd"
      expect_error(compensate(comp.fs, comp), regexp = "must match")
      
      #unmatched length
      comp <- comp[1:3]
      expect_error(compensate(comp.fs, comp), regexp = "must match")
      
      #modify comp[5]
      comp <- sapply(sampleNames(comp.fs), function(sn)comp.mat, simplify = FALSE)
      comp[[5]][2] <- 0.001
      comp.fs4 <- compensate(comp.fs, comp)
      expect_failure(expect_equal(fsApply(comp.fs4, colMeans, use.exprs = TRUE)
              , expectRes), regexp = "8.399298e-06")
      
      #extra comp element
      comp <- sapply(sampleNames(comp.fs), function(sn)comp.mat, simplify = FALSE)
      comp[["dd"]] <- 1:10
      comp.fs5 <- compensate(comp.fs, comp)
      expect_equal(fsApply(comp.fs5, colMeans, use.exprs = TRUE), expectRes)
      
})