#' validity check for samples slot    
#' @param samples \code{character} vector
#' @param \code{list} of objects   
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
  
#' @rdname flowSet-accessor
#' @param filter \code{filter} to be applied
#' @param method \code{missing} not used
#' @param sides \code{missing} not used
#' @param circular \code{missing} not used
#' @param init \code{missing} not used
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




#' @aliases 
#' length,ncdfFlowList-method
#' @rdname flowSet-accessor
setMethod("length",
    signature=signature(x="ncdfFlowList"),
    definition=function(x){
      selectMethod("length", signature = c("ncdfFlowSet"))(x) 
    })
#' @rdname ncdfFlowList-class
#' @param object \code{ncdfFlowList}
#' @aliases 
#' show,ncdfFlowList-method
setMethod("show",
    signature = signature(object="ncdfFlowList"),
    definition = function(object) { 
      cat("An ", class(object), " with", length(object@data), class(object@data[[1]]), "\n")
      cat("containing", length(object), " unique samples.") 
      cat("\n")
    })

#' @rdname flowSet-accessor
#' @aliases 
#' sampleNames,ncdfFlowList-method
setMethod("sampleNames", 
    signature = signature(object = "ncdfFlowList"),
    function(object) {
      names(object@samples)      
    })



#' @export 
#' @rdname ncdfFlowSet-split
#' @aliases split,ncdfFlowList,factor-method
setMethod("split",signature=signature(x="ncdfFlowList",f="factor"),definition=function(x, f, ...)
    {
      
      selectMethod("split", signature = c("ncdfFlowSet", "factor"))(x, f, ...)
      
    })
#' @rdname ncdfFlowSet-split
#' @aliases split,ncdfFlowList,character-method
setMethod("split", signature=signature(x="ncdfFlowList", f="character"), definition=function(x, f, ...)
    {
      selectMethod("split", signature = c("ncdfFlowSet", "character"))(x, f, ...)
    })
#' @aliases
#' phenoData,ncdfFlowList-method
#' phenoData<-,ncdfFlowList,AnnotatedDataFrame-method
#' @export 
#' @rdname flowSet-accessor
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
      
      ncdfFlowList(res)
      res
    })
#' @aliases
#' pData,ncdfFlowList-method
#' pData<-,ncdfFlowList,data.frame-method
#' @param object \code{ncdfFlowList}
#' @rdname flowSet-accessor
#' @export 
setMethod("pData","ncdfFlowList",function(object){
      
      res <- lapply(object,pData, level =1)

      res <- do.call(rbind,res)
      res[sampleNames(object),,drop=FALSE]
    })

#' @rdname flowSet-accessor
#' @exportMethod pData<-
setReplaceMethod("pData",c("ncdfFlowList","data.frame"),function(object,value){
      
      if(!.isValidSamples(rownames(value),object))
        stop("The sample names in data.frame are not consistent with the ",class(x), "!")
      
      res <- lapply(object,function(fs){
            this_pd <- subset(value,rownames(value)%in%sampleNames(fs))
            pData(fs) <- this_pd
            fs
          }, level =1)
      
      ncdfFlowList(res)[sampleNames(object)]        
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

#' @rdname flowSet-accessor
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
