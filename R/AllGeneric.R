#' @export 
setGeneric("updateIndices", 
		function(x,y,z)
			standardGeneric("updateIndices"))
#' @export 
setGeneric("getIndices",function(obj,y,...){
			standardGeneric("getIndices");
		})
#' @export 
setGeneric("initIndices", 
		function(x)
			standardGeneric("initIndices"))

#' @export 
setGeneric("ncdfFlowSet", function(x,...) standardGeneric("ncdfFlowSet"))

#' @export 
setGeneric("ncfsApply",function(x,FUN,...,use.exprs=FALSE,newNcFile=NULL)
			standardGeneric("ncfsApply"))    
#' @export 
setGeneric("unlink")