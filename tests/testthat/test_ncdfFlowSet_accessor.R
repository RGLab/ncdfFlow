context("ncdfFlowSet accessors")
library(flowStats)
library(flowViz)
morphGate <- norm2Filter("FSC-H", "SSC-H", filterId = "MorphologyGate",scale = 2)
fs <- GvHD[pData(GvHD)$Patient %in% 6:7][1:4]
suppressMessages(ncfs <- ncdfFlowSet(fs))
samples <- sampleNames(ncfs)
lgcl <- logicleTransform( w = 0.5, t= 10000, m =4.5)

test_that("keyword", {
  expect_equal(keyword(ncfs, "$TOT"), keyword(fs, "$TOT"))
  keyword(fs) <- list("$FILENAME" = NA, newk = "v")
  keyword(ncfs) <- list("$FILENAME" = NA, newk = "v")
  kl <- list("$FILENAME" , "newk")
  expect_equal(keyword(ncfs, kl), keyword(fs, kl))
  })
test_that("[[", {
      
      sn <- samples[1]
      fr <- ncfs[[sn]]
      expect_is(fr, "flowFrame")
      fr1 <- fs[[sn]]
      expect_equal(fr, fr1, tol = 1e-07)
      fr <- ncfs[[1]]
      expect_equal(fr, fr1, tol = 1e-07)
      
      #without reading data
      fr <- ncfs[[sn, use.exprs = FALSE]]
      fr1@exprs <- matrix(0, nrow = 0, ncol= 0)
      expect_equal(fr, fr1)
      
      #subset by channel
      chnls <- c("FSC-H", "FL2-H")
      fr <- ncfs[[sn, chnls]]
      fr1 <- fs[[sn, chnls]]
      is_equal_flowFrame(fr1, fr)
      
      #subset by int
      chnls <- c(3,5,1)
      fr <- ncfs[[sn, chnls]]
      fr1 <- fs[[sn, chnls]]
      is_equal_flowFrame(fr1, fr)
      
      #subset by single channel 
      chnls <- c(3)
      fr <- ncfs[[sn, chnls]]
      fr1 <- fs[[sn, chnls]]
      is_equal_flowFrame(fr1, fr)
      
    })
test_that("as.flowSet", {
      #Subset by gate
      is_equal_flowSet(Subset(ncfs, rectGate), Subset(fs, rectGate))
      
    })      
test_that("as.flowSet", {
  fs1 <- as.flowSet(ncfs)
  expect_is(fs1, "flowSet")
  expect_equal(colnames(fs1), colnames(ncfs))
  expect_equal(pData(fs1), pData(ncfs))
  
  is_equal_flowSet(ncfs, fs1)
  
})


test_that("unlink", {
      suppressMessages(nc1 <- ncdfFlowSet(GvHD[1]))
      cdfFile <- getFileName(nc1)
      expect_true(file.exists(cdfFile))
      unlink(nc1)
      expect_false(file.exists(cdfFile))
      
    })


test_that("getIndices & Subset", {
      sn <- samples[1]
      
      #initial index is NA
      expect_equal(getIndices(ncfs, sn), NA) 
      
      #subset with filter
      nc1 <- Subset(ncfs, morphGate)
      ind <- getIndices(nc1, sn)
      expect_equal(sum(ind), nrow(nc1[[sn]]))
      expect_equal(length(ind), nrow(ncfs[[sn]]))
      
      #reset indices
      initIndices(nc1)
      expect_equal(getIndices(nc1, sn), NA)
      
      #update ind
      updateIndices(nc1, sn, ind)
      expect_equal(ind, getIndices(nc1, sn))
    })

test_that("[", {
      sn <- samples[1:2]
      nc1 <- ncfs[sn]
      expect_is(nc1, "ncdfFlowSet")
      expect_equal(length(nc1), 2)
      is_equal_flowSet(fs[sn], nc1)
      
      #nc1 and nc share the cdf file
      all.equal(getFileName(nc1), getFileName(ncfs))
      
    })

test_that("subset", {
      
      nc_sub <- subset(ncfs, as.integer(Visit) <= 3)
      is_equal_flowSet(nc_sub, fs[1:3])
      
      nc_sub <- subset(ncfs, as.integer(Visit) <= 3 & Days >=0)
      is_equal_flowSet(nc_sub, fs[2:3])
      
      nc_sub <- subset(ncfs, as.integer(Visit) < 3 | Days == 12)
      is_equal_flowSet(nc_sub, fs[c(1,2,4)])
      
    })



test_that("[[<-", {

      sn <- samples[1]
      suppressMessages(nc <- ncdfFlowSet(fs[sn]))
      
      #return the entire flowFrame
      fr <- nc[[sn]]
      #swap cols of fr
      origcol <- colnames(fr)
      colnames(fr)[7:8] <- origcol[8:7]
      #update nc
      expect_error(nc[[sn]] <- fr, "not consistent")
      #succeed in write after reorder cols
      nc1 <- nc[,colnames(fr)]
      nc1[[sn]] <- fr
      expect_equal(range(nc1[[sn]])[7:8], range(fr)[7:8])
      expect_equal(range(nc1[[sn]], "data")[7:8], range(fr, "data")[7:8])
      
      #transform the data
      #construct transformList first instead of 
      # trransform(fr, `FL1-H` = lgcl(`FL1-H`), `FL2-H` = lgcl(`FL2-H`))
      # because the latter only works in console mode (global envir)
      translist <- transformList(c("FL1-H", "FL2-H"), lgcl)
      
      #list of transformList
      trans.list <- sapply(sampleNames(nc), function(sn)translist)
      trans.fs1 <- transform(nc, trans.list)
      trans_range <- range(trans.fs1[[sn]], "data")
      expect_equal(trans_range[, c("FL1-H")], c(0.6312576, 4.0774226)) 
      expect_equal(trans_range[, c("FL2-H")], c(0.6312576, 3.7131872))
      
      trans.list[[1]] <- logicleTransform()
      expect_error(trans.fs1 <- transform(nc, trans.list), "a valid 'transformList'")
      
      trans.list[[1]] <- translist
      names(trans.list)[1] <- "d"
      expect_error(trans.fs1 <- transform(nc, trans.list), "consistent with flow data")
      
      fr_trans <- transform(fr, translist)
      
      #update the data
      suppressMessages(nc1[[sn]] <- fr_trans)
      trans_range <- apply(exprs(nc1[[sn]]), 2, range)
      expect_equal(trans_range[, c("FL1-H")], c(0.6312576, 4.0774226)) 
      expect_equal(trans_range[, c("FL2-H")], c(0.6312576, 3.7131872))
      
      #subset on channels
      suppressMessages(nc <- ncdfFlowSet(fs[sn]))
      expect_error(nc[[sn]] <- fr_trans[,c("FL1-H")], "colnames of the input are not consistent")
      nc1 <- nc[,c("FL1-H")]
      #only write the channels of interest (reduce disk IO)
      suppressMessages(nc1[[sn]] <- fr_trans[,c("FL1-H")])
      trans_range <- apply(exprs(nc[[sn]]), 2, range)
      #transformed channel
      expect_equal(trans_range[, c("FL1-H")], c(0.6312576, 4.0774226)) 
      #untransformed channel
      expect_equal(trans_range[, c("FL2-H")], c(1.000, 1637.104), tol = 8e-08)
      
      #update chanel colnames
      suppressMessages(nc <- ncdfFlowSet(fs[sn]))
      colnames(fr_trans)[3:4] <- c("<FL1-H>", "<FL2-H>")
       #colnames remain unchanged
      expect_equal(colnames(nc), colnames(ncfs))
      expect_error(nc[[sn]] <- fr_trans, "colnames of the input are not consistent")
      
    })

test_that("ncfsApply", {
      sn <- samples[1]
      #use ncfsApply when FUN returns a flowFrame
      translist <- transformList(c("FL1-H", "FL2-H"), lgcl)
      suppressMessages(nc1 <- ncdfFlow:::ncfsApply(ncfs, transform, translist))
      expect_is(nc1, "ncdfFlowSet")
      expect_equal(sampleNames(ncfs), sampleNames(nc1))
      expect_equal(colnames(ncfs), colnames(nc1))
      #the other channels remain the same
      is_equal_flowSet(ncfs[, -c(3:4)], nc1[, -c(3:4)], description = FALSE)
      #tow channels are tranformed
      trans_range <- apply(exprs(nc1[[sn]]), 2, range)
      expect_equal(trans_range[, c("FL1-H")], c(0.6312576, 4.0774226)) 
      expect_equal(trans_range[, c("FL2-H")], c(0.6312576, 3.7131872))
      expect_false(getFileName(nc1) == getFileName(ncfs))
      
    })    

test_that("sampleNames<-", {
      sn <- samples[1:2]
      nc <- ncfs[sn]
      newNames <- c("s1", "s2")
      sampleNames(nc) <- newNames
      expect_equal(sampleNames(nc), newNames)
      expect_equal(nc@origSampleVector, c(newNames,samples[-c(1:2)]))
      expect_equal(ls(nc@indices), newNames)
      is_equal_flowFrame(ncfs[sn][[1]], nc[[1]])
      
      newNames <- c("s01", "s2")
      sampleNames(nc) <- newNames
      expect_equal(sampleNames(nc), newNames)
      expect_equal(nc@origSampleVector, c(newNames,samples[-c(1:2)]))
      expect_equal(ls(nc@indices), newNames)
      is_equal_flowFrame(ncfs[sn][[1]], nc[[1]])
      
      newNames <- c("s2", "s2")
      expect_error(sampleNames(nc) <- newNames, "Replacement values are not unique")
      
      #replace the single subsetted fs
      nc <- nc["s2"]
      sampleNames(nc) <- "dd"
      expect_equal(sampleNames(nc), "dd")
      expect_equal(nc@origSampleVector, c("s01","dd",samples[-c(1:2)]))
      expect_equal(ls(nc@indices), "dd")
      is_equal_flowFrame(ncfs[sn][[2]], nc[[1]])
      
      #replace with the name that is conflicting with values in origSampleVector
      sampleNames(nc) <- "s01"
      expect_equal(nc@origSampleVector[-1], c("s01",samples[-c(1:2)]))
      is_equal_flowFrame(ncfs[sn][[2]], nc[[1]])
      
      })

test_that("colnames<-", {
      sn <- samples[1:2]
      coln <- colnames(ncfs)
      
      nc <- ncfs[sn, coln[1:2]]
      newColNames <- c("c1", "c2")
      colnames(nc) <- newColNames
      expect_equal(colnames(nc), newColNames)
      expect_equal(nc@origColnames, c(newColNames,coln[-c(1:2)]))
      invisible(fsApply(nc, function(fr)expect_equal(colnames(fr), newColNames)))
      is_equal_flowSet(ncfs[sn, coln[1:2]], nc)
      expect_equivalent(unlist(keyword(nc[[1]])[c("$P1N", "$P2N")]), newColNames)
      
      #change the order of colnames
      nc <- ncfs[sn, coln[2:1]]
      colnames(nc) <- newColNames
      expect_equal(nc@origColnames, c(newColNames[2:1],coln[-c(1:2)]))
      is_equal_flowSet(ncfs[sn, coln[2:1]], nc)
      expect_equivalent(unlist(keyword(nc[[1]])[c("$P1N", "$P2N")]), rev(newColNames))
    })  

test_that("xyplot ncdfFlowSet", {
      sn <- samples[1:2]
      nc <- ncfs[sn]
      ncObj <- xyplot(`SSC-H`~`FSC-H`, nc)
      
      expect_is(ncObj[["panel.args.common"]][["frames"]], "ncdfFlowSet")
      
      expect_equal(ncObj[["panel.args.common"]][["type"]], "xyplot")
      
      fsObj <- xyplot(`SSC-H`~`FSC-H`, fs[sn])
      
      expect_is(fsObj[["panel.args.common"]][["frames"]], "environment")
      
      expect_equal(sub("type", "plotType", deparse(ncObj[["call"]])), deparse(fsObj[["call"]]))
      
      ncObj[["panel.args.common"]][["frames"]] <- NULL
      fsObj[["panel.args.common"]][["frames"]] <- NULL
      ncObj[["panel.args.common"]][["type"]] <- NULL
      ncObj[["call"]] <- NULL
      fsObj[["call"]] <- NULL
      
      expect_equivalent(fsObj, ncObj)
    })

test_that("densityplot ncdfFlowSet", {
      sn <- samples[1:2]
      nc <- ncfs[sn]
      ncObj <- densityplot(~`SSC-H`, nc)
      fsObj <- densityplot(~`SSC-H`, fs[sn])
      
      expect_is(ncObj[["panel.args.common"]][["frames"]], "ncdfFlowSet")
      expect_is(fsObj[["panel.args.common"]][["frames"]], "environment")

      expect_equal(deparse(ncObj[["call"]]), paste0("flowViz:::", deparse(fsObj[["call"]])))
      
      ncObj[["panel.args.common"]][["frames"]] <- NULL
      fsObj[["panel.args.common"]][["frames"]] <- NULL
      ncObj[["call"]] <- NULL
      fsObj[["call"]] <- NULL
      
      expect_equivalent(fsObj, ncObj)
    })

test_that("split", {
      
      #split by factor
      splitBy <- factor(c("p1","p2","p1","p2"))
      
      nclist <- split(ncfs, splitBy)
      fslist <- split(fs, splitBy)
      expect_is(nclist, "list")
      
      expect_equal(names(nclist), names(fslist))
      invisible(lapply(names(nclist), function(thisPop){
                is_equal_flowSet(nclist[[thisPop]], fslist[[thisPop]])
                expect_equal(getFileName(nclist[[thisPop]]), getFileName(ncfs))
              }))
      

      #split by filter
      nclist <- split(ncfs, rectGate)
      fslist <- split(fs, rectGate)
      expect_is(nclist, "list")
      expect_equal(names(nclist), names(fslist))
      invisible(lapply(names(nclist), function(thisPop){
                is_equal_flowSet(nclist[[thisPop]], fslist[[thisPop]])
                expect_equal(getFileName(nclist[[thisPop]]), getFileName(ncfs))
              }))
      
    })


test_that("clone.ncdfFlowSet", {
  
  nc1 <- ncfs[1:2]
  ##clone the ncdfFlowSet object,by default the actual raw data is not added
  nc2 <- clone.ncdfFlowSet(nc1,"clone.nc", isEmpty = TRUE)
  expect_equal(nrow(nc2[[1]]), 0)
  expect_equal(getFileName(nc2), "clone.nc")
  
  #add the actual raw data
  suppressMessages(nc2[[1]] <- nc1[[1]])
  is_equal_flowFrame(nc1[[1]], nc2[[1]])
  
  suppressMessages(nc2 <- clone.ncdfFlowSet(nc1, "clone.nc"))
  is_equal_flowSet(nc1, nc2)
  expect_equal(getFileName(nc2), "clone.nc")
  expect_false(identical(nc2@frames, nc1@frames))
  
  unlink(nc2)
})

