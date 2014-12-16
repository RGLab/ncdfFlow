#' extract a \code{flowFrame} object from \code{ncdfFlowSet}
#' 
#' Simliar to \code{\link[=[[,flowSet-method]{[[}}, and there are cerntain ways to 
#' reduce the disk IO and optimize the speed.
#'  
#' @param x a \code{ncdfFlowSet} or \code{ncdfFlowList}
#' @param i a \code{numeric} or \code{character} used as sample index
#' @param j a \code{numeric} or \code{character} used as channel index
#' @param use.exprs a \code{logical} scalar indicating whether to read the actual data from cdf
#' @param ... other arguments. not used.
#' @export 
#' @rdname extractFlowFrame
#' @aliases [[,ncdfFlowSet,ANY-method
#' @examples 
#' data(GvHD)
#' nc <- ncdfFlowSet(GvHD[1:2])
#' samples <- sampleNames(nc)
#' sn <- samples[1]
#' #return the entire flowFrame
#' fr <- nc[[sn]]  
#' 
#' #access the flowFrame meta data without loading the raw event data from disk
#' nc[[sn, use.exprs = FALSE]]
#' 
#' #only read a subset of channels (more efficient than reading entire data set) 
#' nc[[sn, 1:2]]
#' 
setMethod("[[",
    signature=signature(x="ncdfFlowSet"),
    definition=function(x, i, j, use.exprs = TRUE, ...)
    {
      
      
      cppflag <- globalenv()[["ncdfFlow.[[.cpp"]]
      if(is.null(cppflag))
        cppflag <- TRUE  
      
      if(!cppflag){
        if(length(i) != 1)
          stop("subscript out of bounds (index must have length 1)")
        sampleName<-if(is.numeric(i)) sampleNames(x)[[i]] else i
        fr <- x@frames[[sampleName]]
        
        #get channel index 
        
        localChNames <-colnames(x)
        
        
        #subset by channel
        if(!missing(j)){
          if(is.character(j)){
            j <- match(j, localChNames)
            if(any(is.na(j)))
              stop("subscript out of bounds")
          }  
          #we don't update description slot(i.e. keywords) as flowCore does 
          fr@parameters <- fr@parameters[j, , drop = FALSE]
          localChNames <- localChNames[j]
        }
        
        if(use.exprs){
          origChNames <-x@origColnames ##    
          chIndx <- match(localChNames,origChNames)#only fetch the subset of channels
          
          Indice <- x@indices[[sampleName]]
          if(is.null(Indice))
            stop("Invalid sample name '",sampleName, "'! It is not found in 'indices' slot!")
          subByIndice <- all(!is.na(Indice))
          
          #get sample index
          samplePos <- which(x@origSampleVector==sampleName)
          if(length(samplePos) == 0)
            stop("Invalid sample name '", sampleName, "'! It is not found in 'origSampleVector' slot!")
          
          mat <- .Call(C_ncdfFlow_readSlice, x@file, as.integer(chIndx), as.integer(samplePos), localChNames)
          if(!is.matrix(mat)&&mat==FALSE) stop("error when reading cdf.")
          
          #subset data by indices if neccessary	
          if(subByIndice&&nrow(mat)>0)
            mat<-mat[getIndices(x,sampleName),,drop=FALSE]  
          
          
          
          fr@exprs <- mat
          
        }
        fr
      }else{
        if(missing(j))
          j <- NULL
#    browser()
        readFrame(x, i, j, use.exprs)
#    dd <- readFrame(x, sampleName, j, use.exprs)
#    dd
#    head(exprs(dd))
        
      }
      
    })
    
    
#' @rdname extractFlowFrame
#' @aliases [[,ncdfFlowList,numeric-method
setMethod("[[",c(x="ncdfFlowList",i="numeric"),function(x,i,j, ...){
      
      #convert non-character indices to character
      this_samples <- sampleNames(x)
      nSamples <- length(this_samples)
      if(i > nSamples){
        stop(i, " is larger than the number of samples: ", nSamples)
      }
      
      if(missing(j))#somehow missing j causes trouble for S4 dispatching
        x[[this_samples[i], ...]]
      else
        x[[this_samples[i], j, ...]]
      
    })
#' @rdname extractFlowFrame
#' @aliases [[,ncdfFlowList,logical-method
setMethod("[[",c(x="ncdfFlowList",i="logical"),function(x,i, j, ...){
      #convert non-character indices to character
      if(missing(j))
        x[[sampleNames(x)[i], ...]]
      else
        x[[sampleNames(x)[i], j, ...]]
      
    })
#' @rdname extractFlowFrame
#' @aliases [[,ncdfFlowList,character,missing-method
setMethod("[[",c(x="ncdfFlowList",i="character"),function(x,i, j, ...){
      #convert non-character indices to character
      
      res <- NULL
      object_ind <- tryCatch(
          {
            x@samples[[i]]
          }
          , error = function(cond){
            stop("'", i, "' not found in ", class(x), "!")
          }
      )
#      browser()
      
      object <- x@data[[object_ind]]
      
      if(missing(j))
        res <- object[[i, ...]]
      else
        res <- object[[i, j, ...]]
      
      
    })

    