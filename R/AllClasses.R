#' ncdfFlow: A package that provides CDF storage based flow cytometry data analysis.
#' 
#' Define important flow cytometry data classes:
#'     \code{\link[ncdfFlow:ncdfFlowSet-class]{ncdfFlowSet}}( a subclass of \code{\link[flowCore:flowSet-class]{flowSet}})
#' and \code{\link[ncdfFlow:ncdfFlowList-class]{ncdfFlowList}}(a list of ncdfFlowSet object) and their accessors.
#' 
#' Provide important compensation,transformation,filter,gating,subsetting,splitting functions for data analysis of large volumns of flow cytometry data that is too big to be held in memory.
#'
#'
#' \tabular{ll}{
#'   Package: \tab ncdfFlow\cr
#'   Version: \tab 2.9.24\cr
#'   Date:\tab 2014-04-16\cr
#'   Depends: \tab R (>= 2.8.1), flowCore\cr
#'   License: \tab Artistic-2.0\cr
#'   
#'  }
#'
#' @author
#' Mike Jiang \email{wjiang2@@fhcrc.org},
#' Greg Finak \email{gfinak@@fhcrc.org}
#'
#' Maintainer: Mike Jiang \email{wjiang2@@fhcrc.org}
#' @name ncdfFlow
#' @docType package
#' @title ncdfFlow: A package that provides CDF storage based flow cytometry data analysis.
#' @keywords package
NULL

#' 'ncdfFlowSet': a class for storing flow cytometry raw data in HDF5 format
#' 
#' This class is a subclass of
#'   \code{\link{flowSet}}. It stores the raw data in cdf file instead of memory so that the analysis tools
#'   provided by flowCore based packages can be used in the study that produces hundreds or thousands FCS files.

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
		

setClass("ncdfFlowList"
    ,representation=representation(
        data = "list"
        ,samples="character" #this determine the order of samples exposed to user
    ))

setAs(from = "list", to = "ncdfFlowList", def = function(from){
      
      new("ncdfFlowList", data = from, samples = unname(unlist(lapply(from, sampleNames))))
    })
setAs(from = "ncdfFlowList", to = "flowFrame", def = function(from){
      selectMethod("coerce", signature = c("flowSet", "flowFrame"))(from)      

    })