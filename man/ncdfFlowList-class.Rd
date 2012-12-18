\name{ncdfFlowList-class}
\docType{class}
\alias{ncdfFlowList}
\alias{ncdfFlowList-class}
\alias{[[,ncdfFlowList,ANY-method}
\alias{[,ncdfFlowList,ANY-method}
\alias{ncfsUnlink,ncdfFlowList-method}
\alias{show,ncdfFlowList-method}
\alias{colnames,ncdfFlowList-method}

\title{'ncdfFlowList': a class that stores multiple ncdfFlowSet objects}
\description{
It is a list of ncdfFlowSet objects
}
\section{Objects from the Class}{
Objects can be created using 
new("ncdfFlowList", 
sampleNames = ...., #character of sample names
datalist = .... #a list of ncdfFlowSet objects
)

or visa the constructor,which takes a list of ncdfFlowSet objects: 
ncdfFlowList(ncdfList) 

 
}
\section{Slots}{
  \describe{
    \item{\code{sampleNames}:}{character of sample names }
    	
    \item{\code{datalist}:}{a list of ncdfFlowSet objects }
  }
}
\section{Methods}{
 \describe{
 	\item{[}{return a list of ncdfFlowSets}
 		
 	\item{[[}{return a specific ncdfFlowset object}
    
    \item{show}{display object summary.}
    	
    \item{sampleNames}{Extract and replace sample names from the phenoData object.}
    	
    \item{colnames}{
    		Extract or replace the colnames slot.}
    	
    \item{ncfsUnlink}{
    	delete the cdf files associated with ncdfFlowSet objects contained in current ncdfFlowList }
}
}

\author{
Mike Jiang,Greg Finak,N.Gopalakrishnan

Maintainer: Mike Jiang <wjiang2@fhcrc.org>

}


\seealso{
\code{\link{ncdfFlowSet}}
}
\keyword{classes}
