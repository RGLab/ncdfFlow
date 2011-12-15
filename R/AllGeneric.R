#setGeneric("isEmpty", 
#		function(x,y)
#			standardGeneric("isEmpty"))



setGeneric("addFrame", 
		function(ncfs,data,sampleName)
			standardGeneric("addFrame"))

setGeneric("updateIndices", 
		function(x,y,z)
			standardGeneric("updateIndices"))

#keep this function definition as same as flowWorkspace to avoid conflicts
#setGeneric("getIndices", 
#		function(x,y)
#			standardGeneric("getIndices"))


setGeneric("getIndices",function(obj,y,...){
			standardGeneric("getIndices");
		})

setGeneric("initIndices", 
		function(x,y)
			standardGeneric("initIndices"))

setGeneric("ncdfFlowSet", function(x,...) standardGeneric("ncdfFlowSet"))

setGeneric("ncdfFlowSet_open", function(x,...) standardGeneric("ncdfFlowSet_open"))

setGeneric("ncdfFlowSet_sync", function(x,...) standardGeneric("ncdfFlowSet_sync"))



setGeneric("NcdfFlowSetToFlowSet", function(x,top) standardGeneric("NcdfFlowSetToFlowSet"))


#setGeneric("addObject", function(proj, object,...) standardGeneric("addObject"))



#setGeneric("removeObject", function(proj, object,...) standardGeneric("removeObject"))

#setGeneric("ncdfxyplot", function(data, x, y, sample, ...) standardGeneric("ncdfxyplot"))

#setGeneric("ncdfdensityplot", function(data, x, sample, ...) standardGeneric("ncdfdensityplot"))

setGeneric("ncfsApply",function(x,FUN,...,use.exprs=FALSE,newNcFile=NULL)
			standardGeneric("ncfsApply"))    

setGeneric("ncfsUnlink",
		function(x)
			standardGeneric("ncfsUnlink"),package="ncdfFlow")