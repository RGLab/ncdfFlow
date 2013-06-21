# TODO: Add comment
# 
# Author: mike
###############################################################################
dll<-NULL
.onLoad = function(libname, pkgname){
	dll<<- library.dynam("ncdfFlow",pkgname,lib.loc=.libPaths())
}

validFilterResultList <- function(fres, set, strict=TRUE)
{
	res <- TRUE
	flowCore:::checkClass(fres, "filterResultList")
	flowCore:::checkClass(strict, "logical", 1)
	if(!missing(set)){
		flowCore:::checkClass(set, "ncdfFlowSet")
		if(res <- !all(names(fres) == sampleNames(set)))
			warning("Sample names don't match between flowSet and ",
					"filterResultList", call.=FALSE)
	}
	if(strict){
		fTypes <- sapply(fres, function(x) class(x))
		if(length(unique(fTypes)) != 1){
			warning("Not all filterResults in the list are of equal",
					" type.", call.=FALSE)
			res <- FALSE
		}
		nrPops <- sapply(fres, function(x) length(x))
		if(length(unique(nrPops)) != 1){
			warning("Not all filterResults in the list share the",
					" same number of sub-populations.", call.=FALSE)
			res <- FALSE
		}
		return(res)
	}
}





.guid <- function(len=10){
	ltrs <- c(LETTERS,letters)
	paste(c(sample(ltrs,1),sample(c(ltrs,0:9),len-1,replace=TRUE)),collapse="")
}	


.findTimeCh <- function(x) {
	tm <- grep("^Time$", colnames(x), value = TRUE, ignore.case = TRUE)[1]
	if(!length(tm))
		stop("Unable to identify the time channel")
	tm
}


.checkVals <- function(x, name) {
	
	tmp <- sapply(x, function(k) { as.character( pData(k)[[name]]) })
	nm <- apply(tmp, 1, function(p) { all(p == p[1]) })
	if(!all(nm[!is.na(nm)])) {
		warning(paste("Parameters ", name, " do not match for ",
						paste(colnames(unique(tmp)), sep="", collapse=","), sep=""))
		return(FALSE)
	}else{
		return(TRUE)
	}
	
}

.checkPars <- function(x) {
	
	clNames <- names(pData(x[[1]]))
	pp <- sapply(clNames, function(k) { .checkVals(x, k) })
	if(!all(pp == 0)) return(FALSE)
	else return(TRUE)
	
}

.checkRanges <- function(x) {
	samp <- sampleNames(x)
	rng <- sapply(samp, function(k){
				range(x,k)
			})
	res <- apply(rng, 1, function(k){
				any(k != k[1])
			})
	if(any(res))
		warning("Parameter ranges for some files are not identical")
}

checkParameters <- function(x) {
	
	.checkPars(parameters(x))    
	
}



