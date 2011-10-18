# TODO: Add comment
# 
# Author: mike
###############################################################################
############################################################################################
# three low-level indice IO routines deal with bitvector,no intended to be exposed to user
#ncIndice is a list that contains the filename and node list
#TODO:high-level wrapper is to be written to call these routines to deal with logical vector
##############################################################################################

#len is the eventCount,nlist is the list of node/gate
createIndiceFile<-function(ncFile,eventCount,nlist)
{
	
	msgCreate <- .Call("_createIndiceFile", ncFile,as.integer(eventCount),length(nlist))
			
	if(!msgCreate)stop()
}


##indiceMat is a eventCount x N matrix which represents N contiguous indices vectors
#Node is a gate/node character
writeIndice<-function(ncIndice,indiceMat,startNode) 
{
	if(!is.matrix(indiceMat)||!all(c("bitlen","nbitset")%in%names(attributes(indiceMat))))
	{
		warning("the input needs to have attributes of bitLen and nbitset!")
		return(NULL)
	}
	nodeStartInd<-which(ncIndice$nlist==startNode)
	
	msgWrite <- .Call("_writeIndice",ncIndice$ncFile,indiceMat,as.integer(nodeStartInd))
	if(!msgWrite)stop()
	
	msgWrite
}

#Node is a gate/node character,nodeCount indicates how many contiguous indices vectors is to be read
#return value is a matrix,when single column is retrieved,it is coeced to vector
readIndice<-function(ncIndice,Node,nodeCount=1) 
{
	
	nodeStartInd<-which(ncIndice$nlist==Node)
	
	#	browser()
	indiceMat <- .Call("_readIndice",ncIndice$ncFile,nodeStartInd,as.integer(nodeCount))
	if(!is.matrix(indiceMat)&&!indiceMat)stop()
	
	indx <- nodeStartInd:(nodeStartInd+nodeCount-1)
	clNames <-ncIndice$nlist[indx]
#	colnames(indiceMat) <- clNames
	
	##wrap it to a list of bitvector
	indices<-list()
	for(i in 1:ncol(indiceMat))
	{
#		browser()
		curCol<-indiceMat[,i]
		eval(parse(text=paste("indices$`",clNames[i],"`<-structure(curCol,bitlen=attr(indiceMat,\"bitlen\"),nbitset=attr(indiceMat,\"nbitset\")[i])",sep="")))
		
	}
	
	if(length(indices)==1)
		indices<-indices[[1]]
	
	indices
}

############################################################################################
# the following routines are IO for ncdfFlowSet object to cdf file
##############################################################################################
read.ncdfFlowSet <- function(files = NULL,ncdfFile,flowSetId="",isWriteSlice= TRUE,isSaveMeta=FALSE,phenoData) 
{
	compress<-TRUE##no need for this argument anymore 
	
	if(missing(ncdfFile)) 
		ncdfFile<-tempfile(pattern = "ncfs")
#				 
	
	if (!length(grep(".", ncdfFile, fixed = TRUE)))  
		ncdfFile <- paste(ncdfFile, "nc", sep = ".")
	file.names<-basename(files)
	## obtain event counts and  number of parameters
	bigFile <- files[which.max(file.info(files)[,"size"])]
	tmp <- read.FCSheader(bigFile)[[1]]
	maxEvents <- as.integer(tmp[["$TOT"]])
	maxChannels <- as.integer(tmp[["$PAR"]])
	maxSamples <- length(files)
	### read the big fcs file to get the channel names 
	tmp  <- read.FCS(bigFile)
	channelNames <- as.character(as.character(parameters(tmp)@data[,"name"])) 
	
	#get metaData from each fcs 
	pars <- lapply(seq_len(maxSamples), function(i, verbose){
				tmp <- read.FCS(files[i])
				
				list("description" = description(tmp), "parameters" = parameters(tmp), 
						"guids" = identifier(tmp))
				
			}, verbose = TRUE)
	
	guids <- list()
	#assign metaData to two environment slots
	e1<-new.env(hash=TRUE, parent=emptyenv())
	e2<-new.env(hash=TRUE, parent=emptyenv())
	for(i in seq_len(maxSamples)) {
		
		assign(x=file.names[i],value=new("flowFrame", exprs=matrix(numeric(0),nrow=0,ncol=0), description=pars[[i]][["description"]],
						parameters=pars[[i]][["parameters"]]),envir=e1) 
		assign(file.names[i],rep(TRUE,maxEvents),e2)
		guids[[i]] <- pars[[i]][["guids"]]
		
	}
	if(!missing(phenoData)){
		guids <- sampleNames(phenoData)
	}else{
		guids <- basename(files)
#		browser()
		phenoData = new("AnnotatedDataFrame", data = data.frame(name=guids,row.names=guids), 
				varMetadata = data.frame(labelDescription="Name",row.names="name"))
	}
	
	if(any(duplicated(guids)))
		guids <- make.unique(guids)
	
	#create ncdf ncdf object 
	ncfs<-new("ncdfFlowSet", file = ncdfFile, colnames = channelNames, 
			frames = e1,maxEvents=maxEvents,flowSetId = flowSetId,phenoData=phenoData
			,indices=e2,origSampleVector=guids,origColnames=channelNames)
#		browser()
	
	#get the meta size
	metaSize<-ifelse(isSaveMeta,length(serialize(ncfs,NULL)),0)
	
	

	#create empty cdf file
#	browser()
	
	msgCreate <- .Call("createFile", ncdfFile, as.integer(maxEvents), 
					as.integer(maxChannels), as.integer(maxSamples),
					as.integer(metaSize),as.logical(compress))
	if(!msgCreate)stop()
#	browser()
	##remove indicies to keep the slot as empty by default for memory and speed issue
	initIndices(ncfs,NA)
	#############################################################
	#when isWriteSlice is False,keep the ncdf matrix empty(Mike)
	#############################################################
	if(isWriteSlice)
	{
		lapply(seq_len(maxSamples), function(i, verbose)
				{
					addFrame(ncfs,read.FCS(files[i]),guids[i])
				}, verbose = TRUE)
	}else
	{
		
#		initIndices(ncfs,TRUE)##when isWriteSlice is FALSE,indices have to been explictly initialized
		
	}
#	browser()
	
	
	#write metaData to cdf
	if(isSaveMeta)
		ncdfFlowSet_sync(ncfs)
	
	return(ncfs)
}

##################################################################
#if isEmpty is set as FAlSE,then simply copy the orginal cdf file including the data
#by default,create empty cdf file and then add the data later on 
##################################################################
clone.ncdfFlowSet<-function(ncfs,newNcFile=NULL,isEmpty=TRUE,isNewNcFile=TRUE,isSaveMeta=FALSE)
{
#		browser()
	
	##if flag isNewNcFile=TRUE then creat new ncdf file
	##otherwise keep using the original ncdf file
	if(isNewNcFile)
	{
		if(is.null(newNcFile))
			newNcFile<-tempfile(pattern = "ncfs")

		if (!length(grep(".", newNcFile, fixed = TRUE)))  
			newNcFile <- paste(newNcFile, "nc", sep = ".")
		if(isEmpty)
		{
			metaSize<-ifelse(isSaveMeta,length(serialize(ncfs,NULL)),0)
			msgCreate <- .Call("createFile", newNcFile, as.integer(ncfs@maxEvents), 
							as.integer(length(colnames(ncfs))), as.integer(length(ncfs)),
							as.integer(metaSize),as.logical(TRUE))
			if(!msgCreate)stop()
		}
		
		else
			file.copy(ncfs@file,newNcFile,overwrite=TRUE)
		#update file info
		ncfs@file<-newNcFile	
	}
	
	#update frames info
	orig<-ncfs@frames
	ncfs@frames<-new.env(hash=TRUE, parent=emptyenv())
	for(i in ls(orig))
	{
		assign(i,orig[[i]],ncfs@frames)
	}
	
	if(isEmpty)#when empty init the indices with 
	{
		orig<-ncfs@indices
		ncfs@indices<-new.env(hash=TRUE, parent=emptyenv())
		#copy indices info
		for(i in sampleNames(ncfs))
		{
			
#			updateIndices(ncfs,i,rep(TRUE,length(eval(parse(text=paste("orig$'",i,"'",sep=""))))))
			updateIndices(ncfs,i,NA)
		}
	}else
	{
		orig<-ncfs@indices
		ncfs@indices<-new.env(hash=TRUE, parent=emptyenv())
		#copy indices info
		for(i in ls(orig))
		{
			if(class(orig[[i]])=="raw"){
				updateIndices(x=ncfs,y=i,z=.getBitStatus(orig[[i]]))
			}else{
				updateIndices(x=ncfs,y=i,z=orig[[i]])
			}
#			assign(i,orig[[i]],ncfs@indices)
		}	
	}
	ncfs
}

.writeSlice <- function(ncfs,data,sampleName)
{ 
	
	mat <- exprs(data)
	mode(mat) <- "single"
	#make sure to use origSampleVector for IO since phetaData slot may change after subsetting
	i<-which(ncfs@origSampleVector==sampleName)
	############################################################################################
	#make sure the data matrix contains the same channels as original one before write it back
	#since the C function only write the entire ogirinal slice intead of a subset of it
	#and the current iput mat could be a subset through channels or events
	############################################################################################
	
#	browser()
	#if writing the data slice with the exact size and colnames of original one
	#then simply write the input matrix
	#get original slice
	origMat <- .Call("readSlice", ncfs@file, as.integer(c(1,length(ncfs@origColnames))),as.integer(i))
	if(nrow(origMat)>0)#if ncfs not empty
	{
		#get index of current channels in orignal slice
		colIndex<-which(ncfs@origColnames%in%colnames(ncfs))
		#update original slice
		origMat[,colIndex]<-mat	
	}else
	{
		if(ncol(mat)==length(ncfs@origColnames))
		{
			origMat<-mat
		}else
		{
			stop("You are trying to write a subset of columns to an empty ncdfFlowSet!")	
		}
		
	}
	#write it back to disk
	msgWrite <- .Call("writeSlice", ncfs@file, origMat , as.integer(i))
	
	if(!msgWrite)
	{
#		browser()
		stop("Writing to CDF file failed!")
	}
	msgWrite

}

.retNcdfMat <- function(object, chIndx, sampleName,subByIndice){        
	## chIndx is start and end positions and "readSlice" always get a consecutive chunk 
	#to optimize the reading process
	samplePos<-which(object@origSampleVector==sampleName)
#		browser()	
	mat <- .Call("readSlice", object@file, as.integer(chIndx),as.integer(samplePos))

	
	if(!is.matrix(mat)&&mat==FALSE) stop()
	
	indx <- seq.int(chIndx[1], chIndx[2], 1)
	##get colnames from frame slot instead of flowSet colnames slot since it may not be synchronized
#	clNames <- colnames(eval(parse(text=paste("object@frames$'",sampleName,"'",sep=""))))[indx]  
	clNames <- object@origColnames[indx]
	colnames(mat) <- clNames
#	browser()
	if(subByIndice&&nrow(mat)>0)
		
		mat[getIndices(object,sampleName),,drop=FALSE]
	else
		mat
	
}

