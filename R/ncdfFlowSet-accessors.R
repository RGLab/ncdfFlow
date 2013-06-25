## ==========================================================================
## ncdfFlowSets store the expr matrices to ncdf files and keep meta data in memory
## ==========================================================================
#deFunct
setMethod("addFrame",
		signature=c(ncfs="ncdfFlowSet",data="flowFrame"),
		definition=function(ncfs,data,sampleName){
			.Defunct("[[<-")
#			
			#write to ncdf
			#Since we don't update the indices, we have to make sure we update the correct subset
			ind<-ncdfFlow::getIndices(ncfs,sampleName)
			ddim<-dim(data)
			
#			
			
			#sdim is either size subset or full size
            
			sdim<-dim(origData)
			#if indice is defined,update the subset defined by indices
			if(!any(is.na(ind))){
				updateIndices(ncfs,sampleName,NA)
                slice <- ncfs[[sampleName]]	
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

	
			##update all other slots to keep the whole flowFrame consistent
			eval(parse(text=paste("ncfs@frames$'",sampleName,"'@parameters<-parameters(data)",sep="")))
			eval(parse(text=paste("ncfs@frames$'",sampleName,"'@description<-description(data)",sep="")))
			
		})

as.flowSet <- function(from,top)
    {
      if(!missing(top))
      {
        indice<-round(seq(from=1,to=length(from),length.out=top))
        from<-from[indice]
      }
      frs <- structure(lapply(sampleNames(from),function(n)from[[n]])
          ,names=sampleNames(from))
      fs<-as(frs,"flowSet")
      fs@phenoData<-from@phenoData
      return(fs)
    }
    
setMethod("NcdfFlowSetToFlowSet",
		signature=(x="ncdfFlowSet"),
		definition=function(x,top){
			.Defunct("as.flowSet")
            
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
			definition=function(x,ncdfFile){
			warning("Please convert the flowFrame to flowSet before converting it to ncdfFlowSet!")

		}
)

#load ncdfFlowSet from previous saved ncdf file
setMethod("ncdfFlowSet_open",
		signature=(x="character"),
		definition=function(x){
          .Deprecated()
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
		definition=function(x,ncdfFile){		
			if(missing(ncdfFile))
				ncdfFile <-tempfile(pattern = "ncfs")#ncdfFlow:::.guid() 
			flowSetId = ncdfFile
			
			
			if (!length(grep(".", ncdfFile, fixed = TRUE)))  
				ncdfFile <- paste(ncdfFile, "nc", sep = ".")
			
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
			
#			
			ncfs<-new("ncdfFlowSet", file = ncdfFile, colnames = colnames(x), 
					frames =e1 ,maxEvents=as.integer(maxEvents),flowSetId = flowSetId,
					phenoData= phenoData(x),indices=e2,origSampleVector=sampleNames(x)
					,origColnames=colnames(x))
                
            
            metaSize<-0
			#create new ncdf file			
			msgCreate <-.Call(dll$createFile, ncdfFile, as.integer(ncfs@maxEvents), 
							as.integer(length(colnames(ncfs))), as.integer(length(ncfs)),
							as.integer(metaSize),as.logical(FALSE))
			if(!msgCreate)stop()
                        initIndices(ncfs)			
			for(guid in sampleNames(x))
			{

				ncfs[[guid]] <- x[[guid]]
			}

			ncfs
		})
#save ncdfFlowSet object to ncdf file
setMethod("ncdfFlowSet_sync",
		signature=(x="ncdfFlowSet"),
		definition=function(x,...){
          .Deprecated()
			toWrite<-serialize(x,NULL)
#			
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
#			
			ret<-get(y,obj@indices)
			if(all(!is.na(ret)))
				ret<-.getBitStatus(ret)
			ret			
		})
		
setMethod("initIndices",
		signature=signature(x="ncdfFlowSet"), 
		definition=function(x,y)
		{
			
			for(i in sampleNames(x)){
					updateIndices(x,i,NA)
                  }
		})
		
setMethod("updateIndices",
		signature=signature(x="ncdfFlowSet",y="character",z="logical"), 
		definition=function(x,y,z)
		{
#			
			if(all(!is.na(z)))
				z<-.makeBitVec(length(z),z)
			assign(y,z,x@indices)
		})


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
#						
			
			if(any(is.na(copy)))
				stop("Subset out of bounds", call.=FALSE)
			for(nm in copy)
			{
				#copy the frames and indices for the selected samples	
				assign(nm,x@frames[[nm]],ncfs@frames)


				updateIndices(x=ncfs,y=nm,z=ncdfFlow::getIndices(x,nm))
#				
				#update channels info for each frame
				if(!missing(j))
				{
#					browser()
					##get old AnnotatedDataFrame
					pd<-parameters(ncfs@frames[[nm]])
					#update the parameter by subsetting AnnotatedDataFrame wotj parameter name
					if(is.character(j))
						parameters(ncfs@frames[[nm]]) <- pd[match(j,pData(pd)$name),]
					else
						parameters(ncfs@frames[[nm]]) <- pd[j,]
				}
			}
			 
			#update colnames slot for ncdfFlowSet  
			if(!missing(j)){
				if(is.character(j))
					ncfs@colnames <- colnames(x)[match(j, colnames(x))]

				else
					ncfs@colnames <- colnames(x)[j]
 
				if(any(is.na(colnames(ncfs))))
					stop("Subset out of bounds", call.=FALSE)
			}
#			
			return(ncfs)
		})

setReplaceMethod("colnames",
		signature=signature(x="ncdfFlowSet",
				value="ANY"),
		definition=function(x, value)
		{
#			
			if(length(value) != length(colnames(x)))
				stop("length of new colnames doesn't match with the old one",call.=FALSE)
	
			#get the index of the colnames in the original colnames vector
			colIndex<-which(x@origColnames%in%x@colnames)
			x@colnames <- value#update colnames slot
			x@origColnames[colIndex]<-value#update the original colnames baed on the inex

			##updte colnames of each flowFrames
			for(i in sampleNames(x))
				x@frames[[i]]@parameters@data$name <- value
				
			x
		})		
		

#' extract a \code{flowFrame} object from \code{ncdfFlowSet}
#' 
#' @param x a \code{ncdfFlowSet}
#' @param i a \code{numeric} or \code{character} used as sample index
#' @param j a \code{numeric} or \code{character} used as channel index
setMethod("[[",
		signature=signature(x="ncdfFlowSet"),
		definition=function(x, i, j, ...)
		{
			if(length(i) != 1)
				stop("subscript out of bounds (index must have length 1)")
			
			sampleName<-if(is.numeric(i)) sampleNames(x)[[i]] else i
			fr <- x@frames[[sampleName]]
#						
			subByIndice<-all(!is.na(x@indices[[sampleName]]))
#			browser()
			
			#get channel index 
			origChNames <-x@origColnames ##
            localChNames <-colnames(x)
            if(!missing(j))
              localChNames <- localChNames[j]
			chIndx <- match(localChNames,origChNames)
            
			#get sample index
			samplePos<-which(x@origSampleVector==sampleName)
			mat <- .Call(dll$readSlice, x@file, as.integer(chIndx), as.integer(samplePos))
			if(!is.matrix(mat)&&mat==FALSE) stop("error when reading cdf.")

			##append colnames to matrix
			colnames(mat) <- localChNames
			
			#subset data by indices if neccessary	
			if(subByIndice&&nrow(mat)>0)
				mat<-mat[getIndices(x,sampleName),,drop=FALSE]
			
			fr@exprs <- mat

			
			return(fr)
		})


#' write the flow data from a \code{flowFrame} to \code{ncdfFlowSet} 
#' @param x a \code{ncdfFlowSet}
#' @param i a \code{numeric} or \code{character} used as sample index of \code{ncdfFlowSet}
#' @param j not used
#' @param check.names a \code{logical} indicating whether the colnames of flowFrame
#' should be matched to ncdfFlowSet, it can be set as FALSE in case of parseWorkspace 
#' where the flowFrame with pre-fixed colnames needs to be written to fs
#' yet the colnames of fs is not ready to be updated in the middle of parsing    

setReplaceMethod("[[",
		signature=signature(x="ncdfFlowSet",value="flowFrame"),
		definition=function(x, i, j = "missing", check.names = TRUE, ..., value)
{
       
        #check sample index  
		if(length(i) != 1)
				stop("subscript out of bounds (index must have ",
						"length 1)")
        sampleName <- if(is.numeric(i)) sampleNames(x)[[i]] else i
       
        #validity check for channels in flowFrame
        localChNames <-colnames(x)
        frChNames <- colnames(value)
        if(check.names){
          localChIndx <- match(frChNames,localChNames)
          if(any(is.na(localChIndx)))
            stop("Colnames of the input are not consistent with ncdfFlowSet!", sampleName)  
        }else{
          localChIndx <- 1:length(frChNames)
        }
        

        
        #####################################
        #prepare the data matrix to write
        #####################################
        ncfs <- x[,localChIndx]
        #Since we don't update the indices, we have to make sure we update the correct subset
        ind<-ncdfFlow::getIndices(ncfs,sampleName)
              
        #source data to be updated
        updateIndices(ncfs,sampleName,NA)#clear indices to get the data of original size
        srcFr <- ncfs[[sampleName]]
        srcData<-exprs(srcFr)
        srcCount<-nrow(srcData)
        
        #input data
        newData <- exprs(value)
        newCount<-nrow(newData)
        
        #if indice is defined,extend newData to the original size
        if(all(!is.na(ind))){
          srcData[ind,] <- newData
          newData <- srcData 
        }
        
        if(is.na(ind)){
          origCount <- 0  
        }else{
          origCount <- length(ind) #event count in the orginal cdf
        }
        

        if(newCount == srcCount){
          #update the source with data of the same size
          message("updating ", sampleName , "...")
          
        }else if(newCount == origCount){
          #give the warning when view size doesn't match the new size
          # but matches the original cdf cell couint
          warning("ncdfFlowSet size ", length(ind)
                    , ", view size ", srcCount
                    , " data size ", newCount
                    , sampleName
                  )
        }else if(srcCount == 0)
        {
          #add the data when source event is empty
          message("write ", sampleName, "to empty cdf slot...")
        }
        
        ##################
        #write to ncdf
        ###################
        mode(newData) <- "single"
        #make sure to use origSampleVector for IO since phetaData slot may change after subsetting
        sampleInd<-which(ncfs@origSampleVector==sampleName)
        
        #get original channel ind  
        origChNames <-x@origColnames ##
        if(check.names){
          chIndx <- match(frChNames,origChNames)
          if(any(is.na(chIndx)))
          {
            stop("Colnames of the input are not consistent with ncdfFlowSet!"
                ,sampleName)    
          }
        }else{
          chIndx <- match(colnames(ncfs),origChNames)
        }
        
        
        #write to disk
        msgWrite <- .Call(dll$writeSlice, ncfs@file, newData, as.integer(chIndx), as.integer(sampleInd))
        
        if(!msgWrite)
        {
          stop("Writing to CDF file failed!",sampleName)
        }
        #restore the indices
        updateIndices(ncfs,sampleName,ind);
        
        ##update all other slots to keep the whole flowFrame consistent
        eval(parse(text=paste("x@frames$'",sampleName,"'@parameters<-parameters(value)",sep="")))
        eval(parse(text=paste("x@frames$'",sampleName,"'@description<-description(value)",sep="")))
        
		
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
#						
			lapply(sampleNames(x),function(n) {
								y <- as(x[[n]],"flowFrame")
								y<-FUN(if(use.exprs) exprs(y) else y,...)
    							x[[n]]<-y
							})
			x
		})



## ===========================================================================
## compensate method
## ---------------------------------------------------------------------------
setMethod("compensate",
		signature=signature(x="ncdfFlowSet",
				spillover="ANY"),
		definition=function(x, spillover)
		{
#			
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

				show(object@phenoData)
				cat("\n")
#			}
			cat("  column names:\n  ")
			cat(" ", paste(colnames(object), collapse = ", "))
			cat("\n")
			cat("\n")
 
		})

## ==========================================================================
## Transformation methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("transform",
		signature=signature(`_data`="ncdfFlowSet"),
		definition=function(`_data`,...)
		{
#			
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







