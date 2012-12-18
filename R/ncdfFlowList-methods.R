setMethod("ncfsUnlink",
		signature=signature(x="ncdfFlowList"),
		definition=function(x)
		{
			lapply(x@datalist,function(ncfs)unlink(ncfs@file))
			
		}
)
## to ncdfFlowSet
setMethod("[[",
		signature=signature(x="ncdfFlowList"),
		definition=function(x, i, j, ...)
		{
			if(length(i) != 1)
				stop("subscript out of bounds (index must have length 1)")
#			
			y<-x@datalist
			groupName<-if(is.numeric(i)) names(y)[[i]] else i
			ncfs <- y[[groupName]]
			return(ncfs)
		})

## to the list of ncdfFlowSet
setMethod("[",
		signature=signature(x="ncdfFlowList"),
		definition=function(x, i, j, ...)
		{
			if(length(i) != 1)
				stop("subscript out of bounds (index must have length 1)")
#			
			y<-x@datalist
			groupName<-if(is.numeric(i)) names(y)[[i]] else i
			ncfs <- y[groupName]
			return(ncfs)
		})


#ncdfFlowList  <- function(ncdfFlowList, sampleNames) {
#    nc <- new("ncdfFlowList", ncdfFlowList = ncdfFlowList, sampleNames = sampleNames)
#    nc
#}


ncdfFlowList <- function(ncdfList) {
    sampleNames <- unique(unlist(lapply(ncdfList, sampleNames)))
    nc <- new("ncdfFlowList", sampleNames = sampleNames, datalist=ncdfList)
    nc
}


setMethod("show",
        signature = signature(object="ncdfFlowList"),
        definition = function(object) { 
            cat("An ncdfFlowList with", length(object@datalist),"ncdfFlowSets\n")
            cat("containing", length(sampleNames(object)), " unique samples.") 
            cat("\n")
        })


setMethod("sampleNames", 
        signature = signature(object = "ncdfFlowList"),
        function(object) {
            object@sampleNames      
        })

setReplaceMethod("sampleNames",
		signature = signature(object = "ncdfFlowList"),
		definition = function(object, value)
		{
			oldNames <- sampleNames(object)
			value <- as.character(value)
			if(length(oldNames)!=length(value) ||
					!is.character(value))
				stop(" replacement values must be character vector ",
						"of length equal to number of frames in the set'",
						call.=FALSE)
			if(any(duplicated(value)))
				stop("Replacement values are not unique.", call.=FALSE)

            ##FIXME fix identifier info 
            object@sampleNames <- value
			return(object)
		})

#
#setMethod("stainNames",
#        signature = signature(x = "ncdfFlowList"),
#        function(x) {
#
#            sapply(x@ncdfFlowList, stainNames)
#
#        })


setMethod("colnames",
        signature = signature(x = "ncdfFlowList"),
        function(x) {

            sapply(x@datalist, function(k) {
                    k@colNames
                    } )

        })
