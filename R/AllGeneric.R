setGeneric("updateIndices", 
		function(x,y,z)
			standardGeneric("updateIndices"))

setGeneric("getIndices",function(obj,y,...){
			standardGeneric("getIndices");
		})

setGeneric("initIndices", 
		function(x,y)
			standardGeneric("initIndices"))

setGeneric("ncdfFlowSet", function(x,...) standardGeneric("ncdfFlowSet"))


setGeneric("ncfsApply",function(x,FUN,...,use.exprs=FALSE,newNcFile=NULL)
			standardGeneric("ncfsApply"))    

setGeneric("ncfsUnlink",
		function(x)
			standardGeneric("ncfsUnlink"),package="ncdfFlow")