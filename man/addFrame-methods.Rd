\name{addFrame-methods}
\docType{methods}
\alias{addFrame}
\alias{addFrame-methods}
\alias{addFrame,ncdfFlowSet,flowFrame-method}
\title{add/replace the data in ncdfFlowSet}
\description{
Add one flowFrame to the ncdfFlowSet.(it is defunct, For writing flow data to cdf, use [[<- directly)
}
\usage{
addFrame(ncfs,data,sampleName)
}
\section{Methods}{
  \describe{
  	  
	}
}
\arguments{

	\item{ncfs}{a \code{ncdfFlowSet} object}

 	\item{data}{a \code{flowFrame} to be added}


	\item{sampleName}{
					a \code{character} used as sample index of \code{ncdfFlowSet}  
			}
			

}
\details{
The dimensions of the flowFrame to be added has to match the target sample data in ncdfFlowSet.
It updates the target sample data if it already exists in ncdfFlowSet object.   
}

\author{
Mike Jiang,Greg Finak,N.Gopalakrishnan

Maintainer: Mike Jiang <wjiang2@fhcrc.org>
}

\seealso{
\code{\link{read.ncdfFlowSet}}
}
\examples{
path<-system.file("extdata","compdata","data",package="flowCore")
files<-list.files(path,full.names=TRUE)[1:3]

##create empty ncdfFlowSet from fcs 
nc1  <- read.ncdfFlowSet(files=files,ncdfFile="ncfsTest.nc",flowSetId="fs1",isWriteSlice= FALSE)
fs1<-read.flowSet(files);
#add the actual data slices afterwards
nc1[[1]] <- fs1[[1]]
nc1[[1]]##show the added flowFrame
nc1[[2]]##show empty flowFrame
}
\keyword{methods}

