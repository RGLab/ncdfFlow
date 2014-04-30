#' convert from a \code{ncdfFlowSet} to a \code{flowSet}
#' 
#' The main purpose of this API is to convert the archived data (stored in \code{ncdfFlowSet}) to \code{flowSet}
#' when the speed is more concerned than memory effieciency. 
#' Although \code{ncdfFlowSet} is designed to minimize the disk-IO cost, so usually it is not necessary to do such coersion.  
#'  
#' @param from a \code{ncdfFlowSet}
#' @param top \code{integer} specifies a certain number of samples are evenly selected for the coersion.
#'                            If this argument is missing, then coerce all the samples within the \code{ncdfFlowSet}.
#'                            It is to be used with caution because it can incur the huge memory consumption given  the \code{flowSet} is all-in-memory data structure.    
#' @export 
#' @examples 
#' data(GvHD)
#' nc1 <- ncdfFlowSet(GvHD[1:4])
#' fs <- as.flowSet(nc1)
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
    
#' create ncdfFlowSet from flowFrame (not supported)
#' @rdname read.ncdfFlowSet 
setMethod("ncdfFlowSet",
			signature=(x="flowFrame"),
			definition=function(x,ncdfFile){
			warning("Please convert the flowFrame to flowSet before converting it to ncdfFlowSet!")

		}
)

#' create ncdfFlowSet from flowSet
#' 
#' Normally the \code{ncdfFlowSet} is constructed by loading raw FCS files using \code{read.ncdfFlowSet}.
#' In case there is a legacy \code{flowSet} object, we can convert it to \code{ncdfFlowSet} with this constructor.
#'
#' @param x \code{flowSet}
#' @param ncdfFile \code{character} specifies the file name of cdf file
#' @param dim \code{integer} see details in \link{read.ncdfFlowset}.
#' @param compress \code{integer} see details in \link{read.ncdfFlowset}.
#' @export 
#' @examples 
#' data(GvHD)
#' fs <- GvHD[1:2]
#' ncfs <- ncdfFlowSet(fs)
setMethod("ncdfFlowSet",
		signature=(x="flowSet"),
		definition=function(x,ncdfFile, dim = 2, compress = 0){
          
            dim <- as.integer(match.arg(as.character(dim), c("2","3")))
          
			if(missing(ncdfFile))
				ncdfFile <-tempfile(pattern = "ncfs") 
			flowSetId = ncdfFile
			
			
			if (!length(grep(".", ncdfFile, fixed = TRUE)))  
				ncdfFile <- paste(ncdfFile, "nc", sep = ".")
			
			e1<-new.env(hash=TRUE, parent=emptyenv())

			
			maxEvents <- 0L
            
            
            for(guid in sampleNames(x))
            {
              assign(guid, new("flowFrame",exprs=matrix(numeric(0),nrow=0,ncol=0),parameters(x[[guid]]),description(x[[guid]])), env=e1)
              if(dim == 3)
              {    maxEvents<-max(maxEvents,nrow(exprs(x[[guid]])))				
              }  
            }
			
			
			#assign the maximum number of indices to estimate the ncfs object size
			e2<-new.env(hash=TRUE, parent=emptyenv())
			for(guid in sampleNames(x))
			{
                assign(guid, NA, e2)
#				assign(guid,rep(TRUE,maxEvents),e2)
			}
			
#			
			ncfs<-new("ncdfFlowSet", file = ncdfFile, colnames = colnames(x), 
					frames =e1 ,maxEvents=as.integer(maxEvents),flowSetId = flowSetId,
					phenoData= phenoData(x),indices=e2,origSampleVector=sampleNames(x)
					,origColnames=colnames(x))
                
            
            
			#create new ncdf file			
			msgCreate <-.Call(C_ncdfFlow_createFile, ncdfFile, as.integer(ncfs@maxEvents), 
							as.integer(length(colnames(ncfs))), as.integer(length(ncfs)), dim, as.integer(compress))
			if(!msgCreate)stop()
                        initIndices(ncfs)			
			for(guid in sampleNames(x))
			{

				ncfs[[guid, compress = compress]] <- x[[guid]]
			}

			ncfs
		})

#' delete the cdf file associated with the ncdfFlowSet object
#'         
#' ncdfFlowSet object is unrecoverable after cdf is deleted.
#' So this method is usually called when ncdfFlowSet object is no longer in need.
#' @param x \code{ncdfFlowSet}
#' @param recursive see \link[base:unlink]{unlink}
#' @param force see \link[base:unlink]{unlink}
#' @export 
#' @examples
#' data(GvHD)
#' nc <- ncdfFlowSet(GvHD[1:2])
#' nc[[1]] # data is loaded from cdf file
#' unlink(nc)
#' nc[[1]] # now events since the underlining cdf file is gone        
setMethod("unlink",
		signature=signature(x="ncdfFlowSet"),
		definition=function(x, recursive = FALSE, force = FALSE)
		{
			unlink(x@file, recursive = recursive, force = force)
		}
)

#' extract the event indices of one or multiple samples from ncdfFlowSet
#' 
#' For internal use.
#' 
#' @return a logical vector.
#' @export 
#' @examples 
#' data(GvHD)
#' nc <- ncdfFlowSet(GvHD[1:2])
#' sn <- sampleNames(nc)[1]
#' nrow(nc[[sn]])
#' getIndices(nc, sn) #initial index is NA 
#' #subset with filter
#' morphGate <- norm2Filter("FSC-H", "SSC-H", filterId = "MorphologyGate",scale = 2)
#' nc1 <- Subset(nc, morphGate)
#' ind <- getIndices(nc1, sn)
#' all.equal(sum(ind), nrow(nc1[[sn]]))
#' initIndices(nc1)
#' getIndices(nc1, sn) #reset indices
setMethod("getIndices",
		signature=signature(obj="ncdfFlowSet",y="character"), 
		definition=function(obj,y,...)
		{
		
			ret<-get(y,obj@indices)
			if(all(!is.na(ret)))
				ret<-.getBitStatus(ret)
			ret			
		})
#' initialize the event indices for the entire ncdfFlowSet with NA
#' @export 
setMethod("initIndices",
		signature=signature(x="ncdfFlowSet"), 
		definition=function(x,y)
		{
			
			for(i in sampleNames(x)){
					updateIndices(x,i,NA)
                  }
		})
#' update the event indices of the target sample in ncdfFlowSet
#' @export
#' @param z a \code{logical} vector		
setMethod("updateIndices",
		signature=signature(x="ncdfFlowSet",y="character",z="logical"), 
		definition=function(x,y,z)
		{
		
			if(all(!is.na(z)))
				z<-.makeBitVec(length(z),z)
			assign(y,z,x@indices)
		})

#' get the cdf file name associated with ncdfFlowSet object
#' 
#' @param ncfs \code{ncdfFlowSet}
#' @return \code{character} 
#' @export 
getFileName <- function(ncfs){
  ncfs@file
}

#' subsetting by sampleNames,channels(not for events) methods 
#' 
#' similar to \code{\link[=[,flowSet-method]{[}}.
#'  
#' @export
#' @examples 
#' data(GvHD)
#' nc <- ncdfFlowSet(GvHD[1:2])
#' samples <- sampleNames(nc)
#' nc[1] 
#' nc1 <- nc[samples[1]]
#' #nc1 and nc share the cdf file
#' all.equal(getFileName(nc1), getFileName(nc))
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


				updateIndices(x=ncfs,y=nm,z=getIndices(x,nm))
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


#' extract a \code{flowFrame} object from \code{ncdfFlowSet}
#' 
#' Simliar to \code{\link[=[[,flowSet-method]{[[}}, and there are cerntain ways to 
#' reduce the disk IO and optimize the speed.
#'  
#' @param x a \code{ncdfFlowSet}
#' @param i a \code{numeric} or \code{character} used as sample index
#' @param j a \code{numeric} or \code{character} used as channel index
#' @param use.exprs a \code{logical} scalar indicating whether to read the actual data from cdf
#' @export 
#' @aliases [[,ncdfFlowSet,ANY-method
#' @examples 
#' data(GvHD)
#' nc <- ncdfFlowSet(GvHD[1:2])
#' samples <- sampleNames(nc)
#' sn <- samples[1]
#' #return the entire flowFrame
#' fr <- nc[[sn]]  
#' 
#' #access the flowFrame meta data without loading the raw event data from disk
#' nc[[sn, use.exprs = FALSE]]
#' 
#' #only read a subset of channels (more efficient than reading entire data set) 
#' nc[[sn, 1:2]]
#' 
setMethod("[[",
		signature=signature(x="ncdfFlowSet"),
		definition=function(x, i, j, use.exprs = TRUE, ...)
		{
			if(length(i) != 1)
				stop("subscript out of bounds (index must have length 1)")
			
			sampleName<-if(is.numeric(i)) sampleNames(x)[[i]] else i
			fr <- x@frames[[sampleName]]

            #get channel index 
            origChNames <-x@origColnames ##
            localChNames <-colnames(x)
            
            #subset by channel
            if(!missing(j)){
              if(is.character(j)){
                j <- match(j, localChNames)
                if(any(is.na(j)))
                  stop("subscript out of bounds")
              }
             
              #we don't update description slot(i.e. keywords) as flowCore does 
              fr@parameters <- fr@parameters[j, , drop = FALSE]
              localChNames <- localChNames[j]
            }
             
            

            
            if(use.exprs){
                
                chIndx <- match(localChNames,origChNames)#only fetch the subset of channels
                
                Indice <- x@indices[[sampleName]]
                if(is.null(Indice))
                  stop("Invalid sample name '",sampleName, "'! It is not found in 'indices' slot!")
    			subByIndice <- all(!is.na(Indice))
                
    			#get sample index
    			samplePos <- which(x@origSampleVector==sampleName)
                if(length(samplePos) == 0)
                  stop("Invalid sample name '", sampleName, "'! It is not found in 'origSampleVector' slot!")
                
    			mat <- .Call(C_ncdfFlow_readSlice, x@file, as.integer(chIndx), as.integer(samplePos), localChNames)
    			if(!is.matrix(mat)&&mat==FALSE) stop("error when reading cdf.")
    			
                #subset data by indices if neccessary	
                if(subByIndice&&nrow(mat)>0)
                  mat<-mat[getIndices(x,sampleName),,drop=FALSE]  
                
    			
    			
    			fr@exprs <- mat
          }			
		return(fr)
	})

#' write the flow data from a \code{flowFrame} to \code{ncdfFlowSet}
#'  
#' flowFrame can have less channels than ncdfFlowSet,which is used for partial updating(useful for \code{normalization}) 
#'
#' @name replacement method for ncdfFlowSet
#'  
#' @param x a \code{ncdfFlowSet}
#' @param i a \code{numeric} or \code{character} used as sample index of \code{ncdfFlowSet}
#' @param j not used
#' @param only.exprs a \code{logical} Default is FALSE. which will update the parameters and decriptions slot as well as the raw data.
#'                                  Sometime it is more efficient ti set it to TRUE skip the overhead of colnames matching and updating
#'                                  when user is only concerned about raw data instead of the entire flowFrame.   
#' @param compress \code{integer} It is only relevant to writing slice to '2d' format because the compression is set during the creation of hdf5 file for '3d' format. see details in \link{read.ncdfFlowset}.
#' 
#' @exportMethod [[<-
#' @aliases [[<-,ncdfFlowSet,flowFrame-method 
#' @examples 
#' data(GvHD)
#' nc <- ncdfFlowSet(GvHD[1:2])
#' samples <- sampleNames(nc)
#' sn <- samples[1]
#' #return the entire flowFrame
#' fr <- nc[[sn]]  
#'  
#' apply(exprs(nc[[sn]]), 2, range)
#' 
#' #transform the data
#' lgcl <- logicleTransform( w = 0.5, t= 10000, m =4.5)
#' fr_trans <- transform(fr, `FL1-H` = lgcl(`FL1-H`), `FL2-H` = lgcl(`FL2-H`))
#' 
#' #update the data
#' nc[[sn]] <- fr_trans
#' apply(exprs(nc[[sn]]), 2, range)
#' 
#' #subset on channels
#' nc1 <- nc[,2:3]
#' #only write the channels of interest (reduce disk IO)
#' nc1[[sn]] <- fr_trans[,2:3]
#' 
#' #chanel colnames
#' colnames(fr_trans)[3:4] <- c("<FL1-H>", "<FL2-H>")
#' 
#' #write data without matching up the colnames  
#' nc[[sn, only.exprs = TRUE]] <- fr_trans
setReplaceMethod("[[",
		signature=signature(x="ncdfFlowSet",value="flowFrame"),
		definition=function(x, i, j = "missing", only.exprs = FALSE, compress = 0, ..., value)
{
       
        #check sample index  
		if(length(i) != 1)
				stop("subscript out of bounds (index must have ",
						"length 1)")
        sampleName <- if(is.numeric(i)) sampleNames(x)[[i]] else i
       
        #validity check for channels in flowFrame
        localChNames <-colnames(x)
        frChNames <- colnames(value)
        if(only.exprs){
          localChIndx <- 1:length(frChNames)  
        }else{
          #when need to update other slots in flowFrame
          #make sure the channel names are the same as the ones in ncfs
          if(!setequal(frChNames, localChNames))
            stop("Can't update the entire flowFrame because colnames of the input are not consistent with ncdfFlowSet!"
                , "\n To only update raw data,set only.exprs = TRUE"
            )
          localChIndx <- match(frChNames,localChNames)
#          if(any(is.na(localChIndx)))
#            stop("Not all colnames of the input are present in ncdfFlowSet!", sampleName)          
        }
        
        #####################################
        #prepare the data matrix to write
        #####################################
        ncfs <- x[,localChIndx]
        #Since we don't update the indices, we have to make sure to update the correct subset
        ind <- getIndices(ncfs,sampleName)
              
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
          message("write ", sampleName, " to empty cdf slot...")
        }
        
        ##################
        #write to ncdf
        ###################
#        mode(newData) <- "single"
        #make sure to use origSampleVector for IO since phetaData slot may change after subsetting
        sampleInd <- which(ncfs@origSampleVector==sampleName)
        
        #get original channel ind  
        origChNames <-x@origColnames ##
        if(only.exprs){
          chIndx <- match(colnames(ncfs),origChNames)
        }else{
          chIndx <- match(frChNames,origChNames)
          if(any(is.na(chIndx)))
          {
            stop("Colnames of the input are not consistent with ncdfFlowSet! "
                ,sampleName)    
          }
          
        }
        #write to disk
        msgWrite <- .Call(C_ncdfFlow_writeSlice, ncfs@file, newData, as.integer(chIndx), as.integer(sampleInd), as.integer(compress))
        
        if(!msgWrite)
        {
          stop("Writing to CDF file failed! ",sampleName)
        }
        #restore the indices
        updateIndices(ncfs,sampleName,ind);
        
        ##update all other slots of flowFrame
        ##This is valid only when value has the same colnames as x
        if(!only.exprs){
          x@frames[[sampleName]]@description<-description(value)
          x@frames[[sampleName]]@parameters<-parameters(value)
        }
		
		return(x)
})



#' apply method for ncdfFlowSet (for internal use)
#' 
#' It is equivalent to \code{\link{fsApply}}. But the latter could cause memory issue 
#' when \code{FUN} returns a \code{flowFrame}. \code{ncdfApply} writes to a new cdf file instead of memory. 
#' Thus it will return a ncdfFlowSet object.
#' 
#' When the function given by argument "FUN" does not return the entire flowFrame object with the same 
#' size of the original one (such as compensate,transform...), \code{\link[flowCore:fsApply]{fsApply}} should be used instead.
#' @export 
#' @examples 
#' data(GvHD)
#' nc <- ncdfFlowSet(GvHD[1:2])
#' 
#' #use fsApply when FUN does not return a flowFrame 
#' fsApply(nc, nrow)
#' fsApply(nc, range)
#' 
#' #use ncfsApply when FUN returns a flowFrame
#' lgcl <- logicleTransform( w = 0.5, t= 10000, m =4.5)
#' nc1 <- ncfsApply(nc, transform, `FL1-H` = lgcl(`FL1-H`), `FL2-H` = lgcl(`FL2-H`))
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



#' @rdname ncdfFlowSet-class
#' @export
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

#' @rdname ncdfFlowSet-class
#' @export
setMethod("transform",
    signature=signature(`_data`="ncdfFlowSet"),
    definition=function(`_data`,...)
    {
#			
      newNcFile<-paste(`_data`@file,"trans",sep=".")
      ncfsApply(`_data`,transform,...,newNcFile=newNcFile)
    })




#' @aliases 
#' show,ncdfFlowSet-method
#' @rdname ncdfFlowSet-class
#' @importMethodsFrom methods show
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


# .Note that the replacement method also replaces the GUID for each flowFrame)
# Besides what \code{\link[flowCore:sampleNames<-]{sampleNames<-}} does, it also
# needs to take care of the \code{origSampleVector} and \code{indices} slot.
#
#' @rdname ncdfFlowSet-class
#' @exportMethod sampleNames<-
#' @name sampleNames<-
setReplaceMethod("sampleNames",
    signature=signature(object="ncdfFlowSet"),
    definition=function(object, value)
    {

      oldSampleNames <- sampleNames(object)
      
      #update pData and flowFrame
      object <- callNextMethod()
      
      #update origSampleVector slot
      origSampleVector <- object@origSampleVector
      origSampleVector[match(oldSampleNames, origSampleVector)] <- value
      object@origSampleVector <- origSampleVector
      
   
      #update indices slot
      indEnv <- object@indices
      mapply(oldSampleNames, value, FUN = function(old, new){
            assign(new, indEnv[[old]], indEnv) # copy from old to enw  
            eval(substitute(rm(v, envir = indEnv), list(v = old))) # del the old
            
          })
      
            
      object
    })

# channel names replacement method
# 
# Besides what \code{\link[flowCore:colnames<-]{colnames<-}} does, it also
# needs to update the \code{origColnames} slot.
#' @rdname ncdfFlowSet-class
#' @exportMethod colnames<-
#' @name colnames<-
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
