#' validity check for samples slot        
.isValidSamples<-function(samples,object){
  return (setequal(unlist(lapply(object,sampleNames, level = 1)),samples))
}

#' lapply method for ncdfFlowList
#' 
#' Depending on \code{level} parameter, loop either iterates through the list of ncdfFlowSet objects 
#' or every\code{flowFrame} objects.
#' 
#' 
#' @param X \code{ncdfFlowList} object
#' @param FUN \code{function} to apply
#' @param level \code{numeric}. It controls whether loop at `ncdfFlowSet` level or `sample` level. 
#' when level = 2 (default value),\code{FUN} is applied to each sample. When level = 1, \code{FUN} is applied to each object stored in \code{data} slot.  
#' @param ... other arguments passed to \code{FUN}
#' 
#' @rdname lapply-methods
#' @export 
#' @aliases 
#' lapply,ncdfFlowList-method
setMethod("lapply","ncdfFlowList",function(X,FUN, level = 2,...){
      if(level == 1)
        lapply(X@data,FUN,...)
      else
      {
        sapply(sampleNames(X),function(thisSample,...){
              x <- X[[thisSample]]
              FUN(x, ...)
            }, simplify = FALSE, ...)
      }
    })
#it is not exposed to user to avoid its potential huge memory usage
setMethod("fsApply",
    signature=signature(x="ncdfFlowList",
        FUN="ANY"),
    definition=function(x,FUN,..., simplify = TRUE, use.exprs=FALSE)
    {
      selectMethod("fsApply", signature = c("flowSet"))(x, FUN, ..., simplify = simplify, use.exprs = use.exprs)
      })
  
#' @rdname ncdfFlowList-class
##' @param x \code{ncdfFlowList} object
##' @param filter \code{filter} to be applied
##' @param method \code{missing} not used
##' @param sides \code{missing} not used
##' @param circular \code{missing} not used
##' @param init \code{missing} not used
#' @export 
#' @aliases 
#' filter,ncdfFlowList,filter-method
setMethod("filter",
    signature=signature(x="ncdfFlowList",
        filter="filter"),
    definition=function(x, filter, method = "missing", sides = "missing", circular = "missing", init = "missing")
    {
      selectMethod("filter", signature = c("flowSet", "filter"))(x, filter)
    })

#' @rdname ncdfFlowList-class
#' @param i \code{numeric} index
#' @param j column index
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
#' @rdname ncdfFlowList-class
#' @aliases [[,ncdfFlowList,logical-method
setMethod("[[",c(x="ncdfFlowList",i="logical"),function(x,i, j, ...){
      #convert non-character indices to character
      if(missing(j))
        x[[sampleNames(x)[i], ...]]
      else
        x[[sampleNames(x)[i], j, ...]]
      
    })
#' @rdname ncdfFlowList-class
#' @aliases [[,ncdfFlowList,character,missing-method
setMethod("[[",c(x="ncdfFlowList",i="character"),function(x,i, j, ...){
      #convert non-character indices to character
      
      fr <- NULL
      for(object in x@data){
        this_samples <- sampleNames(object)
        ind <- match(i,this_samples)
        if(!is.na(ind)){
          if(missing(j))
            fr <- object[[ind, ...]]
          else
            fr <- object[[ind, j, ...]]
        }
      }
      if(is.null(fr)){
        stop(i, " not found in ", class(x), "!")
      }else{
        return (fr)
      }
    })




#' @aliases 
#' length,ncdfFlowList-method
#' @rdname ncdfFlowList-class
setMethod("length",
    signature=signature(x="ncdfFlowList"),
    definition=function(x){
      selectMethod("length", signature = c("ncdfFlowSet"))(x) 
    })
#' @rdname ncdfFlowList-class
#' @aliases 
#' show,ncdfFlowList-method
setMethod("show",
    signature = signature(object="ncdfFlowList"),
    definition = function(object) { 
      cat("An ", class(object), " with", length(object@data), class(object@data[[1]]), "\n")
      cat("containing", length(object), " unique samples.") 
      cat("\n")
    })

#' @rdname ncdfFlowList-class
#' @aliases 
#' sampleNames,ncdfFlowList-method
setMethod("sampleNames", 
    signature = signature(object = "ncdfFlowList"),
    function(object) {
      object@samples      
    })

#' @rdname ncdfFlowList-class
setMethod("[",c(x="ncdfFlowList"),function(x,i,j,...){
      
      if(missing(i) && missing(j)) 
        return(x)
      
      samples <- sampleNames(x)
      
      if(missing(i)){
        sampInd <- NULL
        matchInd <- NULL
      }else{
        
        if(is.numeric(i) || is.logical(i)) {
          sampInd <- sampleNames(x)[i]
        }else
          sampInd <- i
        
        noFound <- is.na(sampInd)
        if(any(noFound)){
          stop("sample index out of boundary!")
        }
        matchInd <- match(sampInd,samples)
        noFound <- is.na(matchInd)
        if(length(matchInd) == 0)
          stop("no sample selected!")
        
        if(any(noFound)){
          stop(paste(i[noFound], collapse = " "), " not found in ", class(x), "!")
        }
      }
        
      if(missing(j))
        j <- NULL
      
      res <- lapply(x,function(object){
      
            this_samples <- sampleNames(object)
            if(is.null(sampInd)){
              this_subset <- this_samples
            }else{
              ind <- match(sampInd,this_samples)
              this_subset <- sampInd[!is.na(ind)]  
            }
             
            if(length(this_subset)>0){
              if(is.null(j))
                return (object[this_subset, ...])
              else
                return (object[this_subset, j, ...])
            }else{
              NULL
            }
          }, level =1)
      res <- res[!unlist(lapply(res,is.null))]
      res <- as(res, "ncdfFlowList")
      if(is.null(matchInd))
        res@samples <- samples
      else
        res@samples <- samples[matchInd]
      res
    })


#' @export 
#' @rdname ncdfFlowList-class
#' @aliases split,ncdfFlowList,factor-method
setMethod("split",signature=signature(x="ncdfFlowList",f="factor"),definition=function(x, f, ...)
    {
      
      selectMethod("split", signature = c("ncdfFlowSet", "factor"))(x, f, ...)
      
    })
#' @rdname ncdfFlowList-class
#' @aliases split,ncdfFlowList,character-method
setMethod("split", signature=signature(x="ncdfFlowList", f="character"), definition=function(x, f, ...)
    {
      selectMethod("split", signature = c("ncdfFlowSet", "character"))(x, f, ...)
    })
#' @aliases
#' phenoData,ncdfFlowList-method
#' phenoData<-,ncdfFlowList,AnnotatedDataFrame-method
#' @export 
#' @rdname ncdfFlowList-class
setMethod("phenoData","ncdfFlowList",function(object){
      res <- phenoData(object@data[[1]])
      pData(res) <- pData(object)
      res
    })
#' @exportMethod phenoData<-
setReplaceMethod("phenoData",c("ncdfFlowList","AnnotatedDataFrame"),function(object,value){
      
      if(!.isValidSamples(rownames(value),object))
        stop("The sample names in data.frame are not consistent with the ",class(x), "!")
      
      res <- lapply(object,function(fs){
            this_pd <- value[sampleNames(fs), ]
            phenoData(fs) <- this_pd
            fs
          }, level =1)
      
      res <- as(res, "ncdfFlowList")
      res
    })
#' @aliases
#' pData,ncdfFlowList-method
#' pData<-,ncdfFlowList,data.frame-method
#' @rdname ncdfFlowList-class
#' @export 
setMethod("pData","ncdfFlowList",function(object){
      
      res <- lapply(object,pData, level =1)

      res <- do.call(rbind,res)
      rownames(res) <- res[, "name"]
      res[object@samples,,drop=FALSE]
    })
#' @exportMethod pData<-
setReplaceMethod("pData",c("ncdfFlowList","data.frame"),function(object,value){
      
      if(!.isValidSamples(rownames(value),object))
        stop("The sample names in data.frame are not consistent with the ",class(x), "!")
      
      res <- lapply(object,function(gs){
            this_pd <- subset(value,name%in%sampleNames(gs))
            pData(gs) <- this_pd
            gs
          }, level =1)
      
      res <- as(res, "ncdfFlowList")
      res        
    })
#setReplaceMethod("sampleNames",
#		signature = signature(object = "ncdfFlowList"),
#		definition = function(object, value)
#		{
#			oldNames <- sampleNames(object)
#			value <- as.character(value)
#			if(length(oldNames)!=length(value) ||
#					!is.character(value))
#				stop(" replacement values must be character vector ",
#						"of length equal to number of frames in the set'",
#						call.=FALSE)
#			if(any(duplicated(value)))
#				stop("Replacement values are not unique.", call.=FALSE)
#
#            ##FIXME fix identifier info 
#            object@sampleNames <- value
#			return(object)
#		})

#' @rdname ncdfFlowList-class
#' @export 
setMethod("colnames",
        signature = signature(x = "ncdfFlowList"),
        function(x) {

            cols <- lapply(x, function(k) {
                    colnames(k)
                    }, level = 1)
            cols <- unique(cols)
            if(length(cols) > 1)
              stop("colnames not unique across ncdfFlowSets!")
            cols[[1]]

        })
