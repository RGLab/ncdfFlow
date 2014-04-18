\name{clone.ncdfFlowSet}
\alias{clone.ncdfFlowSet}
\title{
Clone a ncdfFlowSet
}
\description{
Create a new ncdfFlowSet object from an existing one
}
\usage{
clone.ncdfFlowSet(ncfs,ncdfFile=NULL,isEmpty=TRUE,isNew=TRUE)
}
\arguments{
  
  \item{ncfs}{ A \code{\link{ncdfFlowSet}}.}

  \item{isNew}{A logical scalar indicating whether the new cdf file should be created.
  					If FALSE, the original cdf file is associated with the new ncdfFlowSet object.}
  				
  \item{ncdfFile}{A character scalar giving the output file name. By
    default, It is NULL and the function will generate a random
    file name, potentially adding the \code{.cdf} suffix unless a file
    extension is already present. It is only valid when isNewNcFile=TRUE }
  
  \item{isEmpty}{A logical scalar indicating whether the raw data should also be copied.if FALSE,
  				an empty cdf file is created with the same dimensions (sample*events*channels) as the orignial one. }
  
  
}

\value{
  
  A ncdfFlowSet object
  
}

\author{
Mike Jiang,Greg Finak,N.Gopalakrishnan


Maintainer: Mike Jiang <wjiang2@fhcrc.org>
}

\seealso{\code{\link{read.ncdfFlowSet}}}

\examples{
library(ncdfFlow)

path<-system.file("extdata","compdata","data",package="flowCore")
files<-list.files(path,full.names=TRUE)[1:3]

#create ncdfFlowSet from fcs
nc1  <- read.ncdfFlowSet(files=files,ncdfFile="ncfsTest.nc",flowSetId="fs1",isWriteSlice= TRUE)


##clone the ncdfFlowSet object,by default the actual raw data is not added
nc2<-clone.ncdfFlowSet(nc1,"clone.nc")
nc2[[1]]

#add the actual raw data
fs1  <- read.flowSet(files=files)
nc2[[sampleNames(fs1)[1]]] <- fs1[[1]]
nc2[[1]]


#delete the cdf file associated with ncdfFlowSet before removing it from memory
ncfsUnlink(nc2)
rm(nc2)

ncfsUnlink(nc1)
rm(nc1)
}
\keyword{IO}
