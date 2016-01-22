context("ncdfFlowList accessors")

nc_merge <- NULL

suppressMessages(nc1 <- ncdfFlowSet(GvHD[1:2]))
suppressMessages(nc2 <- ncdfFlowSet(GvHD[3:4]))
suppressMessages(nc3 <- ncdfFlowSet(GvHD[5:6]))

nclist <- ncdfFlowList(list(nc1,nc2,nc3))


test_that("as flowFrame", {
      #coerce(collapse) from ncdfFlowList to a single flowFrame
      collapsedData <- as(nclist, "flowFrame")
      expect_is(collapsedData, "flowFrame")
      expect_equal(nrow(collapsedData), sum(unlist(lapply(nclist, nrow))))
      
    })

test_that("rbind2", {
      
      suppressMessages(nc_merge <<- rbind2(nclist))
      is_equal_flowSet(nc_merge, GvHD[1:6])
      expect_false(getFileName(nc_merge) == getFileName(nc1))
      expect_false(getFileName(nc_merge) == getFileName(nc2))
      expect_false(getFileName(nc_merge) == getFileName(nc3))
      
      # merge by using existing cdf
      origFile <- nc_merge@file
      sub1 <- nc_merge[3]
      sub2 <- nc_merge[1]
      tmp <- ncdfFlowList(list(sub1, sub2))
      expect_error(rbind2(tmp, ncdfFile = origFile), "must be provided")
      merge1 <- rbind2(tmp, ncdfFile = origFile, samples = nc_merge@origSampleVector)
      expect_equal(merge1@file, origFile)
      expect_equal(merge1@origSampleVector, nc_merge@origSampleVector)
      is_equal_flowFrame(merge1[[1]], sub1[[1]])
      is_equal_flowFrame(merge1[[2]], sub2[[1]])
      expect_equal(sampleNames(merge1), sampleNames(tmp))
      
      
    })

test_that("lapply", {
      
      invisible(lapply(nclist, FUN= function(i)expect_is(i, "flowFrame")))
      invisible(lapply(nclist, FUN= function(i)expect_is(i, "ncdfFlowSet"), level = 1))
    })

test_that("filter", {
      
      fres1 <- filter(nclist, rectGate)
      fres2 <- filter(nc_merge, rectGate)
      expect_equal(fres1, fres2)
      
    })

test_that("length", {
      expect_equal(length(nclist), length(nc_merge))
      
    })

test_that("colnames", {
      expect_equal(colnames(nclist), colnames(nc_merge))
      
    })


test_that("sampleNames", {
      expect_equal(sampleNames(nclist), sampleNames(nc_merge))
      
    })

test_that("pData", {
      expect_equal(pData(nclist), pData(nc_merge))
      
    })

test_that("pData<-", {
      
      pd  <- pData(nclist)
      pd$g <- letters[1:6]
      pData(nclist) <- pd
      expect_equal(pData(nclist), pd)
      
    })



test_that("[[", {
      sn <- sampleNames(nc_merge)[1]
      is_equal_flowFrame(nc_merge[[1]], nclist[[1]])
      is_equal_flowFrame(nc_merge[[sn]], nclist[[sn]])
    })


test_that("[", {
      sn <- sampleNames(nc_merge)[1:2]
      
      nclist1 <- nclist[sn]
      expect_is(nclist1, "ncdfFlowList")
      
      expect_equal(sampleNames(nclist1), sn)
      
    })

test_that("subset", {
      
      nc_sub <- subset(nclist, as.integer(Visit) <= 3)
      is_equal_flowSet(nc_sub, nclist[1:3])
      
      nc_sub <- subset(nclist, as.integer(Visit) <= 3 & Days >=0)
      is_equal_flowSet(nc_sub, nclist[2:3])
      
      nc_sub <- subset(nclist, as.integer(Visit) < 3 | Days == 12)
      is_equal_flowSet(nc_sub, nclist[c(1,2,4)])
      
    })

test_that("split", {

      nclist1 <- split(nclist, rep(letters[1:2],3))
      expect_is(nclist1, "list")
      invisible(lapply(nclist1, function(i)expect_is(i, "ncdfFlowList")))
      
      is_equal_flowSet(nclist1[[1]],nclist[c(1,3,5)])
      is_equal_flowSet(nclist1[[2]],nclist[c(2,4,6)])
    })

test_that("markernames", {
  expect_true(setequal(markernames(nclist), markernames(nc_merge)))
  
  #create discrepancy within fs
  chnls <- c("FL1-H", "FL3-H")
  markers <- c("CD15", "CD14")
  names(markers) <- chnls
  markernames(nclist@data[[1]][[1]]) <- markers
  expect_warning(res <- markernames(nclist), "not unique")
  suppressWarnings(tmp <- markernames(nclist@data[[1]]))
  tmp <- lapply(tmp, sort)
  expect_equal(res, tmp)
  
  #create discrepancy across fs
  markernames(nclist@data[[1]]) <- markers
  expect_warning(res <- markernames(nclist), "not unique")
  expect_equal(res, tmp)
  
  #setter
  markernames(nclist) <- markers
  expect_equivalent(markernames(nclist)[c(2,1)], markers)
})
