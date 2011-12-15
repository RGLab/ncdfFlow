## ==========================================================================
## ncdfFlowSets store the expr matrices to ncdf files and keep meta data in memory
## ==========================================================================
#addFrame method assumes the ncfs already have all the slots and meta data available
setMethod("addFrame",
		signature=c(ncfs="ncdfFlowSet",data="flowFrame"),
		definition=function(ncfs,data,sampleName){
			
#			browser()
			#write to ncdf
			#Since we don't update the indices, we have to make sure we update the correct subset
			ind<-ncdfFlow::getIndices(ncfs,sampleName)
			ddim<-dim(data)
#						browser()
			
			#sdim is either size subset or full size
			sdim<-dim(ncfs[[sampleName]])
			#if indice is defined,update the subset defined by indices
			if(!any(is.na(ind))){
				updateIndices(ncfs,sampleName,NA)
				slice<-ncfs[[sampleName]]
				tmp<-exprs(slice)
				tmp[ind,]<-exprs(data)
				exprs(slice)<-tmp
				parameters(slice)<-parameters(data)
				description(slice)<-description(data)
				data<-slice
			}
			#checks if sample dimensions matches data dimensions
			#does data size match ncdf file size?
			if(all(sdim==ddim)){
				.writeSlice(ncfs,data,sampleName)
				#No? Then check if we are doing a full replacement. Does data size match full ncdf size (length(ind))
			}else if((ddim[1]==length(ind)&ddim[2]==sdim[2])){
				#warn that we are replacing the full size data
				warning("ncdfFlowSet size ",length(ind),", view size ",dim(slice)[1]," data size ",ddim[1])
				.writeSlice(ncfs,data,sampleName)
			}else if(sdim[1]==0)##when source event is empty then add the frame
			{
				.writeSlice(ncfs,data,sampleName)
			}
			#replace the indices
			updateIndices(ncfs,sampleName,ind);

#			updateIndices(ncfs,sampleName,rep(TRUE,nrow(data)))
			
			##update all other slots to keep the whole flowFrame consistent
			eval(parse(text=paste("ncfs@frames$'",sampleName,"'@parameters<-parameters(data)",sep="")))
			eval(parse(text=paste("ncfs@frames$'",sampleName,"'@description<-description(data)",sep="")))
			
		})


setMethod("NcdfFlowSetToFlowSet",
		signature=(x="ncdfFlowSet"),
		definition=function(x,top){
			
			if(!missing(top))
			{
				indice<-round(seq(from=1,to=length(x),length.out=top))
				x<-x[indice]
			}
			frs <- structure(lapply(sampleNames(x),function(n)x[[n]])
					,names=sampleNames(x))
			fs<-as(frs,"flowSet")
			fs@phenoData<-x@phenoData
			return(fs)
		})
#create ncdfFlowSet from flowFrame
setMethod("ncdfFlowSet",
			signature=(x="flowFrame"),
			definition=function(x,fileName){
			warning("Please convert the flowFrame to flowSet before converting it to ncdfFlowSet!")

		}
)

#load ncdfFlowSet from previous saved ncdf file
setMethod("ncdfFlowSet_open",
		signature=(x="character"),
		definition=function(x){
			ret<-.Call(dll$readMeta, x)
			if(!is.raw(ret)&&!ret)stop()
			if(length(ret)==0)
			{
				warning("Loading Error:no metadata available!")
			}else
				unserialize(ret)
		})

#create ncdfFlowSet from flowSet
setMethod("ncdfFlowSet",
		signature=(x="flowSet"),
		definition=function(x,fileName){		
			if(missing(fileName))
				fileName <-tempfile(pattern = "ncfs")#ncdfFlow:::.guid() 
			flowSetId = fileName
			
			
			if (!length(grep(".", fileName, fixed = TRUE)))  
				fileName <- paste(fileName, "nc", sep = ".")
			
			e1<-new.env(hash=TRUE, parent=emptyenv())

			
			maxEvents<-0
			for(guid in sampleNames(x))
			{
				assign(guid, new("flowFrame",exprs=matrix(numeric(0),nrow=0,ncol=0),parameters(x[[guid]]),description(x[[guid]])), env=e1)
				maxEvents<-max(maxEvents,nrow(exprs(x[[guid]])))				
			}
			
			#assign the maximum number of indices to estimate the ncfs object size
			e2<-new.env(hash=TRUE, parent=emptyenv())
			for(guid in sampleNames(x))
			{
				assign(guid,rep(TRUE,maxEvents),e2)
			}
			
#			browser()
			ncfs<-new("ncdfFlowSet", file = fileName, colnames = colnames(x), 
					frames =e1 ,maxEvents=as.integer(maxEvents),flowSetId = flowSetId,
					phenoData= phenoData(x),indices=e2,origSampleVector=sampleNames(x)
					,origColnames=colnames(x))
			
			metaSize<-length(serialize(ncfs,NULL))
			#create new ncdf file			
			msgCreate <-.Call(dll$createFile, fileName, as.integer(ncfs@maxEvents), 
							as.integer(length(colnames(ncfs))), as.integer(length(ncfs)),
							as.integer(metaSize),as.logical(TRUE))
			if(!msgCreate)stop()
			
			for(guid in sampleNames(x))
			{
				.writeSlice(ncfs,x[[guid]],guid)
			}
#			browser()
			initIndices(ncfs,TRUE)
			ncdfFlowSet_sync(ncfs)
			ncfs
		})
#save ncdfFlowSet object to ncdf file
setMethod("ncdfFlowSet_sync",
		signature=(x="ncdfFlowSet"),
		definition=function(x,...){
			toWrite<-serialize(x,NULL)
#			browser()
			.Call(dll$writeMeta, x@file, toWrite , as.integer(1),as.integer(length(toWrite)))
			
			path<-list(...)$path
			if(!is.null(path))
			{
				if(!identical(dirname(x@file),path))
				{
					file.copy(x@file,path,overwrite=TRUE)
					x@file<-file.path(path,x@file)
				}
			}
			
			print(paste("ncdfFlowSet saved in", x@file))

				
		})


#setMethod("range",
#		signature=signature(x="ncdfFlowSet"),
#		definition=function(x, sample) {
##			browser()
#			if(missing(sample))
#				stop("Please specify a sample for which the range is to be returned")
#			par <- parameters(x@frames[[sample]])
#			rng <- t(par@data[c("minRange", "maxRange")])
#			rownames(rng) <- c("min","max")
#			colnames(rng) <- colnames(x)
#			rng
#			
#		})



### replace a flowFrame
#setReplaceMethod("parameters",
#		signature=signature(object="ncdfFlowSet",
#				value="AnnotatedDataFrame"),
#		definition=function(object, sample, value)
#		{
#			parameters(x,sample)
#			eval(parse(text=paste("object@frames$",sample,"@parameters<-value",sep="")))
#		}
#)



setMethod("ncfsUnlink",
		signature=signature(x="ncdfFlowSet"),
		definition=function(x)
		{
			unlink(x@file)
		}
)





				
setMethod("getIndices",
		signature=signature(obj="ncdfFlowSet",y="character"), 
		definition=function(obj,y,...)
		{
#			browser()
			ret<-get(y,obj@indices)
			if(all(!is.na(ret)))
				ret<-.getBitStatus(ret)
			ret			
		})
		
setMethod("initIndices",
		signature=signature(x="ncdfFlowSet",y="logical"), 
		definition=function(x,y)
		{
			
			for(i in sampleNames(x)){
#				browser()
				curData<-ncdfExprs(object=x,sample=i,subByIndice=FALSE)
				nEvent<-nrow(curData)
				if(is.na(y))
					updateIndices(x,i,y)
				else
					updateIndices(x,i,rep(y,nEvent))
					}
		})
		
setMethod("updateIndices",
		signature=signature(x="ncdfFlowSet",y="character",z="logical"), 
		definition=function(x,y,z)
		{
#			browser()
			if(all(!is.na(z)))
				z<-.makeBitVec(length(z),z)
			assign(y,z,x@indices)
		})
## ==========================================================================
##  we iterate over each frame and provide summaries for those.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#setMethod("summary",
#		signature=signature(object="ncdfFlowSet"), 
#		definition=function(object, ...)
#		{
#			fsApply(object, function(x) apply(exprs(x), 2, summary),
#					simplify=FALSE)
#		})


## ==========================================================================
## subsetting by sampleNames,channels(not for events) methods 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## to ncdfFlowSet
setMethod("[",
		signature=signature(x="ncdfFlowSet"),
		definition=function(x, i, j, ..., drop=FALSE)
		{
			
			
			if(missing(i) && missing(j)) 
				return(x)
			
			newChnLen<-length(if(missing(j))colnames(x) else j)
			newSampLen<-length(if(missing(i))sampleNames(x) else i)
			
			#copy ncdfFlowSet object
			ncfs<-x
			#init two environment
			ncfs@frames<-new.env(hash=TRUE, parent=emptyenv())##create new frame env
			ncfs@indices<-new.env(hash=TRUE, parent=emptyenv())##create new frame env
			
			#update frames by samples
			if(missing(i)) {
				copy<-sampleNames(x)

			} else {
				#update samples info in phenoData(can't use phenoData<- due to the validity check)
				ncfs@phenoData <- phenoData(x)[i,]
				if(is.numeric(i) || is.logical(i)) {
					copy <- sampleNames(x)[i]
				} else {
					copy <- i
					i <- match(i,sampleNames(x))
				}
			}
			
			if(any(is.na(copy)))
				stop("Subset out of bounds", call.=FALSE)
			for(nm in copy)
			{
				#copy the frames and indices for the selected samples	
				assign(nm,x@frames[[nm]],ncfs@frames)
#				assign(nm,getIndices(x,nm),ncfs@indices)

				updateIndices(x=ncfs,y=nm,z=getIndices(x,nm))
#				browser()
				#update channels info for each frame
				if(!missing(j))
				{
					##get old AnnotatedDataFrame
					pd<-parameters(ncfs@frames[[nm]])
					#update the parameter by subsetting AnnotatedDataFrame wotj parameter name
					parameters(ncfs@frames[[nm]]) <- pd[pData(pd)$name%in%j,]
				}
			}
			 
			#update channels info for ncdfFlowSet  
			if(!missing(j)){
				if(is.character(j))
					colnames(ncfs) <- colnames(x)[match(j, colnames(x))]
				else
					colnames(ncfs) <- colnames(x)[j] 
				if(any(is.na(colnames(ncfs))))
					stop("Subset out of bounds", call.=FALSE)
			}
#			browser()
#			if(isNew)
#				ncfs<-clone.ncdfFlowSet(ncfs,isEmpty=FALSE,isNew=TRUE)
			return(ncfs)
		})

setReplaceMethod("colnames",
		signature=signature(x="ncdfFlowSet",
				value="ANY"),
		definition=function(x, value)
		{
			colIndex<-which(x@origColnames%in%x@colnames)
			x@colnames <- value
			x@origColnames[colIndex]<-value

			for(i in sampleNames(x))
				x@frames[[i]]@parameters@data$name <- value
				#colnames(x@frames[[i]]) <- value
			x
		})		
		
#
#setMethod("isEmpty",
#		signature=signature(x="ncdfFlowSet"),
#		definition=function(x,y)
#		{
#			if(missing(y))
#				return(length(x@indices)==0)
#			
#			if(length(y) != 1)
#				stop("subscript out of bounds (index must have length 1)")
##			browser()
#			sampleName<-if(is.numeric(y)) sampleNames(x)[[y]] else y
#			
#			return(is.null(eval(parse(text=paste("x@indices$'",sampleName,"'",sep="")))))
#			
#		})


## ==========================================================================
## subsetting by sampleName(channels)(not for events) methods 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## to flowFrame
setMethod("[[",
		signature=signature(x="ncdfFlowSet"),
		definition=function(x, i, j, ...)
		{
			if(length(i) != 1)
				stop("subscript out of bounds (index must have length 1)")

			sampleName<-if(is.numeric(i)) sampleNames(x)[[i]] else i
			fr <- x@frames[[sampleName]]
#						browser()
			subByIndice<-all(!is.na(x@indices[[sampleName]]))
#			browser()
			
#			if(!isEmpty(x,i))
#			{
				fr@exprs <- ncdfExprs(x, sampleName,subByIndice=subByIndice)
#			}
			if(!missing(j))
				fr <- fr[,j]
			return(fr)
		})



## to flowFrame
#setMethod("$",
#          signature=signature(x="ncdfFlowSet"),
#          definition=function(x, name) 
#		  {
##			  browser()
#			  x[[name]]
#  })

## replace a flowFrame
setReplaceMethod("[[",
		signature=signature(x="ncdfFlowSet",
				value="flowFrame"),
		definition=function(x, i, j, ..., value)
		{
			if(length(i) != 1)
				stop("subscript out of bounds (index must have ",
						"length 1)")
                        cnx <- colnames(x)
                        cnv <- colnames(value)
                        if(length(cnx) != length(cnv) || !all(sort(cnv) == sort(cnx)))
                            stop("The colnames of this flowFrame don't match ",
                                 "the colnames of the flowSet.")
                        
			sel <- if(is.numeric(i)) sampleNames(x)[[i]] else i
			
#			x@frames[[sel]]@parameters <- value@parameters
#			x@frames[[sel]]@description <- value@description
#			ncdfFlow:::.writeSlice(x,value,sel)

			addFrame(x,value,sel)
			
			return(x)
		})



## ==========================================================================
## apply method for ncdfFlowSet,do not return direct results,
## instead,creating new nc file,save the results in nc file,and return new ncdfFlowSet
## use this method only when the FUN returns a flowFrame,otherwise use fsApply
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("ncfsApply",
		signature=signature(x="ncdfFlowSet",
				FUN="ANY"),
		definition=function(x,FUN,...,use.exprs=FALSE,newNcFile=NULL)
		{
			
			if(missing(FUN))
				stop("ncfsApply function missing")
			FUN <- match.fun(FUN)
			if(!is.function(FUN))
				stop("This is not a function!")
			x<-clone.ncdfFlowSet(x,newNcFile,isEmpty=FALSE)
#			browser()			
			lapply(sampleNames(x),function(n) {
								y <- as(x[[n]],"flowFrame")
								y<-FUN(if(use.exprs) exprs(y) else y,...)
								ncdfFlow:::.writeSlice(x,y,n)
							})
			x
		})


#setReplaceMethod("colnames",
#		signature=signature(x="ncdfFlowSet",
#				value="ANY"),
#		definition=function(x, value)
#		{
#			x@colnames <- value
#			for(i in sampleNames(x))
#				x@frames[[i]]@parameters@data$name <- value
#			x
#		})

## ===========================================================================
## compensate method
## ---------------------------------------------------------------------------
setMethod("compensate",
		signature=signature(x="ncdfFlowSet",
				spillover="ANY"),
		definition=function(x, spillover)
		{
#			browser()
			newNcFile<-paste(x@file,"comp",sep=".")
			ncfsApply(x, compensate, spillover,newNcFile=newNcFile)
			
		}

)

setMethod("show",
		signature=signature(object="ncdfFlowSet"),
		definition=function(object)
		{ 
			cat("An ncdfFlowSet with", length(sampleNames(object)),"samples.\n")
			cat("flowSetId :", object@flowSetId, "\n") 
			cat("NCDF file :", object@file, "\n")
#			if(any(varMetadata(object@phenoData)$labelDescription != "Name")){
				show(object@phenoData)
				cat("\n")
#			}
			cat("  column names:\n  ")
			cat(" ", paste(colnames(object), collapse = ", "))
			cat("\n")
#			cat("  maxEvents:  ")
#			cat(object@maxEvents)
#			cat("\n")
			
			
			
#			cat("  stain names:\n   ")
#			cat(" ", paste(stainNames(object), collapse = ", "))
			cat("\n")
			#checkParameters(object) 
		})

## ==========================================================================
## Transformation methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("transform",
		signature=signature(`_data`="ncdfFlowSet"),
		definition=function(`_data`,...)
		{
#			browser()
			newNcFile<-paste(`_data`@file,"trans",sep=".")
			ncfsApply(`_data`,transform,...,newNcFile=newNcFile)
		})




# TODO:
## ===========================================================================
## spillover method
## ---------------------------------------------------------------------------
setMethod("spillover",
          signature=signature(x="ncdfFlowSet"),
          definition=function(x, unstained=NULL, patt=NULL, fsc="FSC-A",
                              ssc="SSC-A", method="median", useNormFilt=FALSE)
      {
          if(is.null(unstained)) {
              stop("Sorry, we don't yet support unstained cells blended ",
                   "with stained cells", call.=FALSE)
          } else {
              ## We often only want spillover for a subset of the columns 
              allcols <- colnames(x)
              cols <- if(is.null(patt)) allcols else grep(patt, allcols,
                                                          value=TRUE)

              ## Ignore these guys if they somehow got into cols.
              ## cols <- cols[-match(c(fsc,ssc),cols)]
              cols <- cols[!(cols %in% c(fsc,ssc))]
              
              ## There has got to be a better way of doing this...
              if(!is.numeric(unstained)) {
                  unstained <- match(unstained,sampleNames(x))
                  if(is.na(unstained))
                      stop("Baseline not in this set.", call.=FALSE)
              }
              ## Check to see if the unstained sample is in the list of
              ## stains. If not, we need to add it, making it the first
              ## row and adjust the unstained index accordingly.
              ## If it is there we adjust to the appropriate index.
              
              ## pdh: you shouldn't use the nor2Filter as a default without telling people!
              if(useNormFilt){
                  if(is.numeric(fsc)) fsc <- allcols[fsc]
                  if(is.numeric(ssc)) ssc <- allcols[ssc]
                  
                  if(is.na(match(fsc,allcols)))
                      stop("Could not find forward scatter parameter. ",
                           "Please set the fsc parameter", call.=FALSE)
                  if(is.na(match(ssc,allcols)))
                      stop("Could not find side scatter parameter. ",
                           "Please set the ssc parameter", call.=FALSE)
                  n2f <- norm2Filter(fsc, ssc, scale.factor=1.5)
                  x <- Subset(x,n2f)
              }
              ## inten <- fsApply(Subset(x, n2f), each_col,method)[, cols]
              inten <- fsApply(x, each_col,method)[, cols]
              inten <- pmax(sweep(inten[-unstained,], 2,inten[unstained,]), 0)
              inten <- sweep(inten, 1,apply(inten, 1, max), "/")
              row.names(inten) <- colnames(inten)[apply(inten ,1, which.max)]
              inten[colnames(inten),]
          }
      })







