\name{ncdfFlowSet-class}
\docType{class}
\alias{ncdfFlowSet-class}
\alias{ncdfFlowSet}
\alias{[,ncdfFlowSet,ANY-method}
\alias{[[,ncdfFlowSet,ANY-method}
\alias{[[<-,ncdfFlowSet,ANY,ANY,flowFrame-method}
\alias{colnames<-,ncdfFlowSet-method}
\alias{compensate,ncdfFlowSet,ANY-method}
\alias{getIndices}
\alias{getIndices,ncdfFlowSet,character-method}
\alias{initIndices}
\alias{initIndices,ncdfFlowSet-method}
\alias{ncdfFlowSet_sync}
\alias{ncdfFlowSet_sync,ncdfFlowSet-method}
\alias{NcdfFlowSetToFlowSet}
\alias{NcdfFlowSetToFlowSet,ncdfFlowSet-method}
\alias{as.flowSet}
\alias{ncfsApply}
\alias{ncfsApply,ncdfFlowSet-method}
\alias{ncfsUnlink}
\alias{ncfsUnlink,ncdfFlowSet-method}
\alias{rbind2}
\alias{rbind2,ncdfFlowSet,flowFrame-method}
\alias{rbind2,ncdfFlowSet,ncdfFlowSet-method}
\alias{show,ncdfFlowSet-method}
\alias{split}
\alias{split,ncdfFlowSet,factor-method}
\alias{split,ncdfFlowSet,filter-method}
\alias{split,ncdfFlowSet,filterResultList-method}
\alias{split,ncdfFlowSet,list-method}
\alias{split,ncdfFlowSet,character-method}
\alias{Subset}
\alias{Subset,ncdfFlowSet,filter-method}
\alias{Subset,ncdfFlowSet,filterResultList-method}
\alias{Subset,ncdfFlowSet,list-method}
\alias{transform,ncdfFlowSet-method}
\alias{updateIndices}
\alias{updateIndices,ncdfFlowSet,character,logical-method}
\alias{ncdfFlowSet_open}
\alias{ncdfFlowSet_open,character-method}
\alias{ncdfFlowSet,flowFrame-method}
\alias{ncdfFlowSet,flowSet-method}
\alias{densityplot,formula,ncdfFlowSet-method}
\alias{xyplot,formula,ncdfFlowSet-method}

\title{'ncdfFlowSet': a class for storing flow cytometry raw data in netCDF format}
\description{
This class is a subclass of
  \code{\link{flowSet}}. It stores the raw data in cdf file instead of memory so that the analysis tools
  provided by flowCore based packages can be used in the study that produces hundreds or thousands FCS files.  
}
\section{Objects from the Class}{
Objects can be created by using:
\code{\link{read.ncdfFlowSet}},\code{\link{clone.ncdfFlowSet}}
 
or 
\code{ncdfFlowSet(
			x,					#x is a flowSet
			ncdfFile			#the output cdf file name
			)
 }
 
}
\section{Slots}{
  \describe{
    \item{\code{file}:}{A character containing the ncdf file name.}
    \item{\code{maxEvents}:}{An integer containing the maximum number of events of all samples stored in this ncdfFlowSet object }
    \item{\code{flowSetId}:}{A character for the id of ncdfFlowSet object }
    \item{\code{indices}:}{Object of class \code{"environment"} containing events indices of each sample stored as \code{"raw"} vector. Each index value is either TURE or FALSE and the entire indices vector is used to subset the raw data.
    						the indices vector of each sample is NA by default when the ncdfFlowSet first created.It is assigned with actual value when ncdfFlowSet object is subsetted by \code{\link{Subset}}
								 or other subsetting methods.}
    \item{\code{origSampleVector}:}{A character vector containing the sample names,
    								which indicates the original order of samples physically stored in cdf format}
	\item{\code{origColnames}:}{A character vector containing the flow channel names,
    								which indicates the original order of columns physically stored in cdf format}

    \item{\code{frames}:}{Object of class \code{"environment"} 
    	which replicates the "frame" slot in \code{\link{flowSet}},
    	except that  \code{\link[=exprs,flowFrame-method]{exprs}} matrix is empty and the actual data is stored in cdf file. }
    \item{\code{phenoData}:}{see \code{\link[=phenoData,flowSet-method]{phenoData}}}
    \item{\code{colnames}:}{see \code{\link[=colnames,flowSet-method]{colnames}}. Here it serves as the current data view which does not reflect the actual number and order of columns stored in cdf file.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{flowSet}"}, directly.
}
\section{Methods}{
Most of the other flowSet methods are not listed here,but they all work as the same due to its inheritance from flowSet.Please see for more \code{\link{flowSet}}details for these methods.
  \describe{
    \item{addFrame}{add or replace the flowFrame in ncdfFlowSet.See \code{\link{addFrame}} for more details
		
				     \emph{Usage:}
		
				      \code{   addFrame(ncfs,data,sampleName)}
      				}
      
    \item{[,[[}{Subsetting. similar to \code{\link[=[,flowSet-method]{[}},\code{\link[=[[,flowSet-method]{[[}} for flowSet.}
	\item{[[<-}{replace data with a \code{flowFrame}. 
	      
	      \emph{Usage:}
	      
	      \code{   ncfs[[thisSample]] <- value}
	
	      \code{   ncfs[[thisSample,j = c("SSC-A","FSC-A")], check.names = TRUE, only.exprs = FALSE] <- value }
	      Note that providing the channel index through \code{j} can speed up writing process since only the specified channels are written.
	      Also \code{check.names} flag indicates whether the colnames of flowFrame should be matched to ncdfFlowSet, it can be set as FALSE 
							Thus simply update the first n channels wihtin ncdfFlowSet without matching channel names
							It is useful in \code{\link{flowWorkspace:parseWorkspace}} where the flowFrame with pre-fixed colnames needs to be written to ncdfFlowSet
							 where the colnames has not yet ready to be updated in the middle of parsing
							  
	      \code{only.exprs} option allows for only updating the exprs data when it is set to TRUE.
		   otherwise, the parameters and decriptions slot within flowFrame are updated as well.
	      
	    }
    
    \item{getIndices}{extract the event indices of one or multiple samples from ncdfFlowSet,return a logical vector.
				     	
				     	\emph{Usage:}
				     	
				      	\code{   getIndices(ncfs,sampleName)}
				      }
    \item{initIndices}{initialize the event indices for the entire ncdfFlowSet with NA 
			    		
			    		\emph{Usage:}
				      	
				      	\code{   initIndices(ncfs)}
			    		}
    \item{updateIndices}{update the event indices of the target sample in ncdfFlowSet
			    		
			    		\emph{Usage:}
			    	  	
			    	  	\code{   updateIndices(ncfs,sampleName,y)}	
			    	}
    \item{ncdfFlowSet_open}{(Deprecated!)load the ncdfFlowSet object from the cdf file,return a ncdfFlowSet object
    						
					    	\emph{Usage:}
					      
					     	 \code{   ncdfFlowSet_open(filename)}
     	 
					     	Note that in order to successfully recover the 
					     	entire ncdfFlowSet object,the phenoData has to be already saved 
					     	in cdf either by explicitly calling \code{\link{ncdfFlowSet_sync}}
							or setting \code{isSaveMeta} as TRUE when ncdfFlowSet is 
							created by \code{link{read.ncdfFlowSet}} or \code{link{clone.ncdfFlowSet}}
     	 					}
    \item{ncdfFlowSet}{create ncdfFlowSet from a flowSet object 
    						
					    	\emph{Usage:}
					      
					     	 \code{   ncdfFlowSet(fs1)}
     	 
					     	Note that only flowSet can be coerced to ncdfFlowSet,attempt to apply this method to flowFrame will get an error message.
     	 					}
    \item{ncdfFlowSet_sync}{(Deprecated!)save phenoData to cdf file.
    						
    						\emph{Usage:}
      
     						 \code{   ncdfFlowSet_sync(filename)}
     	 					}
    \item{ncfsUnlink}{delete the netCDF file associated with the ncdfFlowSet object
					   	 
					   	 \emph{Usage:}
					      
					      \code{   ncfsUnlink(ncfs)}
     	 			
				     	 Note that ncdfFlowSet object is unrecoverable after cdf is deleted.
				     	 So this method is usually called when ncdfFlowSet object is no longer in need.
     	 				}

    \item{as.flowSet}{convert a subset of ncdfFlowSet to flowSet.
    
						  		  \emph{Usage:}
						      
						     	 \code{as.flowSet(ncfs,top)}
						     	 
						     	 Argument \code{top} specifies the number of samples evenly selected
						     	 from ncdfFlowSet. 
    							}
    \item{ncfsApply}{equivalent to \code{\link{fsApply}}., which could cause memory issue due to the.
    				So ncfsApply will return a ncdfFlowSet object.
    				 
			    	 \emph{Usage:}
			      
			     	 \code{ncfsApply(x="ncdfFlowSet",FUN="ANY")}
			     	 
			     	 Note that the function given by argument "FUN" 
    				should return an entire flowFrame object with the same 
    				size of the original one(such as compensate,transform...)
			     	 Otherwise,\code{\link[flowCore:fsApply]{fsApply}} should be used instead. 
					}
    \item{rbind2}{similar to \code{\link[=rbind2,flowSet,flowSet-method]{flowCore:rbind2}}. It returns a new ncdfFlowSet with a new cdf file that combines two raw datasets.
     				It is recommended to construct a ncdfFlowList and apply rbind2 to it directly when combining more than two ncdfFlowSets.
     				Because using "do.call" on a list ncdfFlowSets will create one cdf file for every pair of ncdfFlowSets,which is not efficient.      
     				}
   
    \item{split}{equivalent to \code{\link[=split,flowSet,ANY-method]{flowCore:split}}. It returns a new ncdfFlowSet object without creating new cdf raw data file but assigning logical indices to subset the original raw data.}
   
    \item{Subset}{equivalent to \code{\link[=Subset,flowSet,ANY-method]{flowCore:Subset}}. It returns a new ncdfFlowSet object without creating new cdf raw data file but assigning logical indices to subset the original raw data.}    
    
    \item{densityplot,xyplot}{equivalent to \code{\link[=densityplot,formula,flowSet-method]{flowViz:densityplot}} and
    							\code{\link[=xyplot,formula,flowSet-method]{flowViz:xyplot}}.
     						User need to be careful about applying these plot methods to ncdfFlowSet because it could be slow for large number of flow data.}
    
    
	 }
}
\author{
Mike Jiang,Greg Finak,N.Gopalakrishnan


Maintainer: Mike Jiang <wjiang2@fhcrc.org>

}


\seealso{
\code{\link{flowSet}},\code{\link{read.ncdfFlowSet}},\code{\link{ncdfFlowList}}
}
\examples{
data(GvHD)

nc1 <- ncdfFlowSet(GvHD[1:2])
nc2 <- ncdfFlowSet(GvHD[3:4])
nc3 <- ncdfFlowSet(GvHD[5:6])

##combine two ncdfFlowSets
nc4 <- rbind2(nc1,nc2)

##combine more than two ncdfFlowSets
ncfslist <- as(list(nc1,nc2,nc3),"ncdfFlowList")
nc5 <- rbind2(ncfslist)

}
\keyword{classes}
