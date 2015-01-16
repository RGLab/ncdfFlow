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
      
        if(missing(j))
          j <- NULL

        readFrame(x, i, j, use.exprs)
      
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

    