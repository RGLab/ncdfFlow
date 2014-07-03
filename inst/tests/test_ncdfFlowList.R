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

test_that("split", {

      nclist1 <- split(nclist, rep(letters[1:2],3))
      expect_is(nclist1, "list")
      invisible(lapply(nclist1, function(i)expect_is(i, "ncdfFlowList")))
      
      is_equal_flowSet(nclist1[[1]],nclist[c(1,3,5)])
      is_equal_flowSet(nclist1[[2]],nclist[c(2,4,6)])
    })