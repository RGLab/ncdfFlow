#setMethod("ncfsUnlink",
#		signature=signature(x="ncdfFlowList"),
#		definition=function(x)
#		{
#			lapply(x@datalist,function(ncfs)unlink(ncfs@file))
#			
#		}
#)
### to ncdfFlowSet

#' validity check for samples slot        
.isValidSamples<-function(samples,object){
  return (setequal(unlist(lapply(object,sampleNames, level = 1)),samples))
}


#' @param \code{ncdfFlowList} object
#' @param FUN \code{function} to apply
#' @param level \code{numeric}. When \code{X} is a \code{ncdfFlowList}, \code{level} 2 (default value)
#' \code{FUN} is applied to each element of list in \code{data} slot. When it is set to 1, \code{FUN} is applied to each list stored in \code{data} slot.  
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



setMethod("[[",c(x="ncdfFlowList",i="numeric"),function(x,i,j, ...){
      
      #convert non-character indices to character
      this_samples <- sampleNames(x)
      nSamples <- length(this_samples)
      if(i > nSamples){
        stop(i, " is larger than the number of samples: ", nSamples)
      }
      x[[this_samples[i], j, ...]]
      
    })

setMethod("[[",c(x="ncdfFlowList",i="logical"),function(x,i, j, ...){
      #convert non-character indices to character
      
      x[[sampleNames(x)[i], j, ...]]
      
    })
setMethod("[[",c(x="ncdfFlowList",i="character"),function(x,i, j, ...){
      #convert non-character indices to character
      
      fr <- NULL
      for(object in x@data){
        this_samples <- sampleNames(object)
        ind <- match(i,this_samples)
        if(!is.na(ind)){
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
#' @rdname length-methods
setMethod("length",
    signature=signature(x="ncdfFlowList"),
    definition=function(x){
      selectMethod("length", signature = c("ncdfFlowSet"))(x) 
    })

setMethod("show",
    signature = signature(object="ncdfFlowList"),
    definition = function(object) { 
      cat("An ", class(object), " with", length(object@data), class(object@data[[1]]), "\n")
      cat("containing", length(object), " unique samples.") 
      cat("\n")
    })


setMethod("sampleNames", 
    signature = signature(object = "ncdfFlowList"),
    function(object) {
      object@samples      
    })


setMethod("[",c(x="ncdfFlowList",i="missing"),function(x,i,j, ...){
      i <- sampleNames(x)
      x[i,j, ...]
    })
setMethod("[",c(x="ncdfFlowList",i="numeric"),function(x,i, j, ...){
 
        sampleInd <- sampleNames(x)[i]
        noFound <- is.na(sampleInd)
        
        if(any(noFound)){
          stop("sample ", paste(i[noFound], collapse = ""), " not found in ", class(x), "!")
        }
        x[sampleInd, j, ...]
    })

setMethod("[",c(x="ncdfFlowList",i="logical"),function(x,i, j, ...){
        x[sampleNames(x)[i], j, ...]
    })
setMethod("[",c(x="ncdfFlowList",i="character"),function(x,i,j,...){

      if(missing(j))
        j <- NULL
      
      samples <- sampleNames(x)
      matchInd <- match(i,samples)
      noFound <- is.na(matchInd)
      if(any(noFound)){
        stop(i[noFound], " not found in ", class(x), "!")
      }
      
      res <- lapply(x,function(object){

            this_samples <- sampleNames(object)
            ind <- match(i,this_samples)
            this_subset <- i[!is.na(ind)] 
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
      res@samples <- samples[matchInd]
      res
    })



setMethod("split",signature=signature(x="ncdfFlowList",f="factor"),definition=function(x, f, ...)
    {
      
      selectMethod("split", signature = c("ncdfFlowSet", "factor"))(x, f, ...)
      
    })

setMethod("split", signature=signature(x="ncdfFlowList", f="character"), definition=function(x, f, ...)
    {
      selectMethod("split", signature = c("ncdfFlowSet", "character"))(x, f, ...)
    })


#' @aliases
#' pData,ncdfFlowList-method
#' pData<-,ncdfFlowList,data.frame-method
#' @rdname pData-methods
setMethod("pData","ncdfFlowList",function(object){
      
      res <- lapply(object,pData, level =1)

      res <- do.call(rbind,res)
      rownames(res) <- res[, "name"]
      res[object@samples,,drop=FALSE]
    })

setReplaceMethod("pData",c("ncdfFlowList","data.frame"),function(object,value){
      browser()
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
