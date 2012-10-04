## Split a ncdfFlowSet by a single filter, by first creating a list of
## filterResult and then working our way through that in the next
## method.
setMethod("split",
		signature=signature(x="ncdfFlowSet",
				f="filter"),
		definition=function(x, f, drop=FALSE, population=NULL,
				prefix=NULL, ...)
		{
			fres <- filter(x,f)
#			
			
			split(x, fres, population=population, prefix=prefix,...)
		})

setMethod("split",
		signature=signature(x="ncdfFlowSet",
				f="filterResultList"),
		definition=function(x, f, drop=FALSE, population=NULL,
				prefix=NULL, ...)
		{
			
			frameNames<-f@frameId
			f <- f@.Data
			names(f) <-frameNames 
			
			split(x, f, drop=drop, population=NULL, prefix=NULL, ...)
		})

setMethod("split",
		signature=signature(x="ncdfFlowSet",
				f="list"),
		definition=function(x, f,isNew=FALSE, drop=FALSE, population=NULL,
				prefix=NULL, ...)
		{
			
			
#			
				
			
			sample.name <- sampleNames(x)
			lf <- length(f)
			lx <- length(x)
			if(lf!=lx)
				stop("list of filterResults or filters must be same ",
						"length as flowSet.", call.=FALSE)
			if(!all(sapply(f, is, "filter")))
				stop("Second argument must be list of filterResults or filters,",
						call.=FALSE)
			lapply(f, flowCore:::compatibleFilters,  f[[1]])
			## split everything or just some populations
			## (if multipleFilterResult)
			if(is.null(population)){
				if(!is.null(names(f[[1]])))
					population <- names(f[[1]])
				else
					population <- c("positive", "negative")
			} else if(!all(sapply(population, is, "character")))
				stop("'population' must be a single character vector ",
						"or a list of character vectors", call.=FALSE)
			if(!is.list(population)){
				n <- population
				population <- as.list(population)
				names(population) <- n
			}
			## FIXME: Do we want to allow for different names when splitting
			## flowSets by multipleFilterResults?
			if(lf>1 && !identical(unique(as.vector(sapply(f, names))),
					names(f[[1]]))){
				for(i in 2:lf)
					names(f[[i]]) <- names(f[[1]])
				warning("Filtering operation produced non-unique population ",
						"names.\n  Using names of the first frame now.\n",
						"  Please check parameter descriptions in the ",
						"parameter slots\n  of the individual flowFrames.",
						call.=FALSE)
			}
			
			
			finalRes <- vector(mode="list", length=length(population))
			names(finalRes) <- names(population)
			
			for(p in seq_along(population)){
				tp <- population[p]
#				res <- vector(mode="list", length=lf)
				ncfs<-clone.ncdfFlowSet(x,isNew = FALSE)
				for(i in 1:lf){
					tp <- unlist(tp)
					curMultiFilterResult<-f[[i]]
					curFilterResult<-curMultiFilterResult[[tp]]
#					
					
					indice<-x[[i]]%in%curFilterResult
					curSampleName<-names(f)[i]
					##get original indice vector 
					rawIndice<-getIndices(x,curSampleName)
					##update original indice vector with the new subset indice which is shorter than original one
#					rawIndice[which(rawIndice)]<-rawIndice[which(rawIndice)]&indice
					##update original indice vector with the new subset indice which is shorter than original one
					if(all(is.na(rawIndice)))
						rawIndice<-indice
					else
						rawIndice[which(rawIndice)]<-indice
					updateIndices(ncfs,curSampleName,rawIndice)
					
				}
				np <- names(population)[p]
				if(isNew)
					ncfs<-clone.ncdfFlowSet(ncfs,isEmpty = FALSE,isNew=TRUE)
				
				finalRes[[np]] <- ncfs
				phenoData(finalRes[[np]])$population <- np
				varMetadata(finalRes[[np]])["population", "labelDescription"] <-
						"population identifier produced by splitting"
			}
#			
			res<-ncdfFlowList(finalRes)
			return(res)
		})

## Split by frames of flowSet according to a factor, character or numeric.
## Those have to be of the same length as the flowSet. We can't allow for
## drop=TRUE, because this would create invalid sets.
setMethod("split",
		signature=signature(x="ncdfFlowSet",
				f="factor"),
		definition=function(x, f,isNew=FALSE, drop=FALSE, ...)
		{
			
			
			
			if(!is.atomic(f) || length(f)!=length(x))
				stop("split factor must be same length as flowSet",
						call.=FALSE) 
			gind <- split(1:length(f), f, drop=TRUE)
			res <- vector(mode="list", length=length(gind))

			for(g in seq_along(gind)){
				ncfs<-x[sampleNames(x)[gind[[g]]]]
				if(isNew)
					ncfs<-clone.ncdfFlowSet(ncfs,isEmpty=FALSE,isNew=TRUE)
				res[[g]] <-ncfs
				
				phenoData(res[[g]])$split <- levels(f)[g]
				varMetadata(res[[g]])["split", "labelDescription"] <-
						"Split"
			}
			names(res) <- names(gind)
			res<-ncdfFlowList(res)
			return(res)
		})

setMethod("split",
		signature=signature(x="ncdfFlowSet",
				f="character"),
		definition=function(x, f, ...)
		{
			split(x,as.factor(f),...)
		})	