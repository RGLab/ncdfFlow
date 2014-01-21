\name{ncdfFlow-package}
\alias{ncdfFlow-package}
\alias{ncdfFlow}
\docType{package}
\title{
ncdfFlow: A package that provides netCDF storage based flow cytometry data analysis.
}
\description{
Define important flow cytometry data classes:
  \code{\link[ncdfFlow:ncdfFlowSet-class]{ncdfFlowSet}}( a subclass of \code{\link[flowCore:flowSet-class]{flowSet}})
  and \code{\link[ncdfFlow:ncdfFlowList-class]{ncdfFlowList}}(a list of ncdfFlowSet object) and their accessors.
  
  Provide important compensation,transformation,filter,gating,subsetting,splitting functions for data analysis of large volumns of flow cytometry data that is too big to be held in memory.

}
\details{
\tabular{ll}{
Package: \tab ncdfFlow\cr
Version: \tab 0.99.4\cr
Date:\tab 2011-08-02\cr
Depends: \tab R (>= 2.8.1), flowCore\cr
License: \tab Artistic-2.0\cr

}
}


\author{
Mike Jiang,Greg Finak,N.Gopalakrishnan

Maintainer: Mike Jiang <wjiang2@fhcrc.org>
}
\references{
\url{http://www.rglab.org/}
}
\keyword{ package }


