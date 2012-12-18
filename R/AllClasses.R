
setClass("ncdfFlowSet",                   
		representation=representation(
				file = "character",
				maxEvents="integer",
				flowSetId = "character",
				indices = "environment",
				origSampleVector="character",
				origColnames="character"),
		contains="flowSet",
		prototype=list(frames=new.env(hash=TRUE, parent=emptyenv()),
						phenoData=new("AnnotatedDataFrame",
								data=data.frame(),
								varMetadata=data.frame()),
#				colnames=character(0),
				file = character(0),
				maxEvents=integer(0),
				flowSetId = character(0),
				indices=new.env(hash=TRUE, parent=emptyenv()),
				origSampleVector=character(0),
				origColnames=character(0)
				)
#		,validity=function(object){
#		return(TRUE)
#		}
		)
		




				
setClass("ncdfFlowList",
        representation = representation(
                sampleNames = "character",
				datalist = "list"),
        prototype = list(sampleNames = character(0),
				datalist=list()),
        validity = function(object) {

        })
#
#
#setClass("projDetails",
#         representation = representation(
#                name = "character",
#                flowSet = "ncdfFlowSet",
#                processTree = "graph",
#                opsObject = "list"))
#
#setClass("flowProject", 
#          representation = representation(env = "environment"))

    
