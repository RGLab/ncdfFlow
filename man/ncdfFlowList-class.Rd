\name{ncdfFlowList-class}
\docType{class}
\alias{ncdfFlowList}
\alias{ncdfFlowList-class}
\alias{show,ncdfFlowList-method}
\alias{rbind2,ncdfFlowList-method}
\alias{rbind2,ncdfFlowList,ANY-method}

\title{'ncdfFlowList': a class that stores multiple ncdfFlowSet objects}
\description{
It is a list of ncdfFlowSet objects
}
\section{Objects from the Class}{
Objects can be created by coercing a list of ncdfFlowSet objects 
as("ncdfFlowList",nclist = .... #a list of ncdfFlowSet objects)

}
\section{Methods}{
 \describe{
 	
    \item{show}{display object summary.}
     \item{rbind2}{combine a list of ncdfFlowSets into one.}
     }
    
}

\author{
Mike Jiang,Greg Finak,N.Gopalakrishnan

Maintainer: Mike Jiang <wjiang2@fhcrc.org>

}


\seealso{
\code{\link{ncdfFlowSet}}
}
\examples{
data(GvHD)
GvHD <- GvHD[pData(GvHD)$Patient %in% 6:7][1:4]
nc1<-ncdfFlowSet(GvHD)

##splitting by a gate
qGate <- quadGate(filterId="qg", "FSC-H"=200, "SSC-H"=400)
fr<-filter(nc1,qGate)
ncList<-split(nc1,fr)
ncList 
nc1[1]#each ncdfFlowSet in the ncdfFlowList share the same cdf file with the original nc1 
nc1[[1]]
ncList[[2]][[1]]
ncList[[1]][[1]]

##create the new cdf file in ncdfFlowList object thus seperate its data source from the original one
ncList<-split(nc1,fr,isNew=TRUE)
ncList[1]


## split the ncdfFlowSet by a factor
ncList<-split(nc1[1:3], as.factor(c(1,1,2)))
ncList
ncList[1]
ncList[[1]][[1]]
summary(ncList[[1]])

}
\keyword{classes}
