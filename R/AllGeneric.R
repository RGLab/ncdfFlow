#' @export 
setGeneric("updateIndices", 
		function(obj,y,z)
			standardGeneric("updateIndices"))
#' @export 
setGeneric("getIndices",function(obj,y,...){
			standardGeneric("getIndices");
		})
#' @export 
setGeneric("initIndices", 
		function(obj)
			standardGeneric("initIndices"))

#' @export 
setGeneric("ncdfFlowSet", function(x,...) standardGeneric("ncdfFlowSet"))

#' @export 
setGeneric("ncfsApply",function(x,FUN,...,use.exprs=FALSE,newNcFile=NULL)
			standardGeneric("ncfsApply"))    
#' @export 
setGeneric("unlink")