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

#' @import zlibbioc 
NULL

#' @useDynLib "ncdfFlow"

NULL


#' a class for storing flow cytometry raw data in HDF5 format
#' 
#' This class is a subclass of
#'   \code{\link{flowSet}}. It stores the raw data in cdf file instead of memory so that the analysis tools
#'   provided by flowCore based packages can be used in the study that produces hundreds or thousands FCS files.
#' 
#' @section Slots:
#'  \describe{
#' 
#' \item{\code{file}:}{A character containing the ncdf file name.}
#' \item{\code{maxEvents}:}{An integer containing the maximum number of events of all samples stored in this ncdfFlowSet object }
#' \item{\code{flowSetId}:}{A character for the id of ncdfFlowSet object }
#' \item{\code{indices}:}{Object of class \code{"environment"} containing events indices of each sample stored as \code{"raw"} vector. Each index value is either TURE or FALSE and the entire indices vector is used to subset the raw data.
#'                         the indices vector of each sample is NA by default when the ncdfFlowSet first created.It is assigned with actual value when ncdfFlowSet object is subsetted by \code{\link{Subset}}
#'                          or other subsetting methods.}
#' \item{\code{origSampleVector}:}{A character vector containing the sample names,
#'                                  which indicates the original order of samples physically stored in cdf format}
#' \item{\code{origColnames}:}{A character vector containing the flow channel names,
#'                              which indicates the original order of columns physically stored in cdf format}
#' 
#' \item{\code{frames}:}{Object of class \code{"environment"}, which replicates the "frame" slot in \code{\link{flowSet}},
#'                      except that  \code{\link[=exprs,flowFrame-method]{exprs}} matrix is empty and the actual data is stored in cdf file. }
#' \item{\code{phenoData}:}{see \code{\link[=phenoData,flowSet-method]{phenoData}}}
#' \item{\code{colnames}:}{see \code{\link[=colnames,flowSet-method]{colnames}}. Here it serves as the current data view which does not reflect the actual number and order of columns stored in cdf file.}
#' }
#' 
#' @section Extends: 
#' Class \code{"\linkS4class{flowSet}"}, directly.
#' 
#' @exportClass ncdfFlowSet
#' @rdname ncdfFlowSet-class
#' @importClassesFrom Biobase AnnotatedDataFrame
#' @importClassesFrom flowCore flowFrame flowSet
#' @importMethodsFrom Biobase description description<- exprs exprs<- pData pData<- phenoData phenoData<- sampleNames sampleNames<- varLabels varMetadata varMetadata<- AnnotatedDataFrame
#' @importMethodsFrom flowCore colnames colnames<- compensate filter fsApply identifier %in% ncol nrow parameters parameters<- split Subset transform
#' @import methods BiocGenerics
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
				file = character(0),
				maxEvents=integer(0),
				flowSetId = character(0),
				indices=new.env(hash=TRUE, parent=emptyenv()),
				origSampleVector=character(0),
				origColnames=character(0)
				)
		)
		
#' a class that stores multiple ncdfFlowSet objects
#' 
#' It is a list of ncdfFlowSet objects
#' 
#' @section Objects from the Class:
#' Objects can be created by coercing a list of ncdfFlowSet objects 
#' as("ncdfFlowList",nclist = .... #a list of ncdfFlowSet objects)
#' 
#' @section Slots:
#'  \describe{
#' 
#' \item{\code{data}:}{A list containing the \link{ncdfFlowSet} objects.}
#' \item{\code{samples}:}{A \code{integer} vector containing the index of the \link{ncdfFlowSet} object to which each sample belongs.
#'                        The name of the vector is the sample names that determine the order of samples exposed to the user, which can
#'                        be different from the physical storing order.}
#' }
#' @seealso \code{\link{ncdfFlowSet}}
#' @exportClass ncdfFlowList
#' @examples 
#' data(GvHD)
#' nc1 <- ncdfFlowSet(GvHD[1])
#' nc2 <- ncdfFlowSet(GvHD[2])
#' nc3 <- ncdfFlowSet(GvHD[3])
#' list1 <- list(nc1, nc2, nc3)
#' #coerce from list to ncdfFlowList
#' nclist <- ncdfFlowList(list1)
#' nclist
#' #coerce(collapse) from ncdfFlowList to a single flowFrame
#' collapsedData <- as(nclist, "flowFrame")
#' collapsedData
#' @rdname ncdfFlowList-class
setClass("ncdfFlowList"
    ,representation=representation(
        data = "list"
        ,samples="integer" # named integer this determine the order of samples exposed to user
    ))

#' create the sample index 
#' 
#' The index is a named integer vector used for fast indexing
#' 
#' @param x \code{list} of objects
.indexingSample <- function(x){
       
  unlist(lapply(seq_along(x), function(i){
            sn <- sampleNames(x[[i]])
            ind <- rep(i, length(sn))
            names(ind) <- sn
            ind
          }))
  
}
#' constuctor for \code{ncdfFlowList}
#' @param samples \code{integer} see \code{samples} slot of \code{ncdfFlowList} class.
#'                  or \code{character} that specifiy the order to samples.
#'                                    If not given then reconstruct the index.
#' @return \code{ncdfFlowList-class} 
#' @rdname ncdfFlowList-class
#' @param x \code{list} of \code{ncdfFlowSet} objects
#' @export 
ncdfFlowList <- function(x, samples = NULL){
  
      if(is.null(samples)){
        
        sampleIndex <- .indexingSample(x)
        
      }else if(is.character(samples))
      {
          sampleIndex <- .indexingSample(x)
          #reorder by the given samples vector
          thisInd <- match(samples, names(sampleIndex))
          sampleIndex <- sampleIndex[thisInd]
          
       }else if(is.integer(samples)){
         sampleIndex <- samples
      }else
        stop("invalid sampleIndex!")
      
      new("ncdfFlowList", data = x, samples = sampleIndex)
      
    }
setAs(from = "ncdfFlowList", to = "flowFrame", def = function(from){
      selectMethod("coerce", signature = c("flowSet", "flowFrame"))(from)      

    })