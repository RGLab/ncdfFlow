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
	
	msgCreate <- .Call(dll$"_createIndiceFile", ncFile,as.integer(eventCount),length(nlist))
			
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
	
	msgWrite <- .Call(dll$"_writeIndice",ncIndice$ncFile,indiceMat,as.integer(nodeStartInd))
	if(!msgWrite)stop()
	
	msgWrite
}

#Node is a gate/node character,nodeCount indicates how many contiguous indices vectors is to be read
#return value is a matrix,when single column is retrieved,it is coeced to vector
readIndice<-function(ncIndice,Node,nodeCount=1) 
{
	
	nodeStartInd<-which(ncIndice$nlist==Node)
	
	#	
	indiceMat <- .Call(dll$"_readIndice",ncIndice$ncFile,nodeStartInd,as.integer(nodeCount))
	if(!is.matrix(indiceMat)&&!indiceMat)stop()
	
	indx <- nodeStartInd:(nodeStartInd+nodeCount-1)
	clNames <-ncIndice$nlist[indx]
#	colnames(indiceMat) <- clNames
	
	##wrap it to a list of bitvector
	indices<-list()
	for(i in 1:ncol(indiceMat))
	{
#		
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
read.ncdfFlowSet <- function(files = NULL
								,ncdfFile,flowSetId=""
								,isWriteSlice= TRUE
								,isSaveMeta=FALSE
								,phenoData
								,channels=NULL
								,...) #dots to be passed to read.FCS
{
#	
	#remove nonexisting files
	fileInd<-file.exists(files)
	missingFiles<-files[!fileInd]
	if(length(missingFiles)>0)
		message(paste(missingFiles,"is missing\n"))
	files<-files[fileInd]
	if(length(files)==0)
	{
		print("No file is found!")
		return(NULL)
	}
		
	compress<-FALSE##no need for this argument anymore 
	
	if(missing(ncdfFile)) 
		ncdfFile<-tempfile(pattern = "ncfs")
#				 
	
	if (!length(grep(".", ncdfFile, fixed = TRUE)))  
		ncdfFile <- paste(ncdfFile, "nc", sep = ".")
	
	ncdfFile<-path.expand(ncdfFile)
	
	file.names<-basename(files)
#	browser()
	nFile<-length(files)
	events<-vector("integer",nFile)
	channels.all<-vector("character",nFile)
	message("Determine the max total events by reading FCS headers ...")
	for(i in 1:nFile){
		curFile<-files[i]
		txt<-flowCore:::read.FCSheader(curFile)[[1]]
		nChannels <- as.integer(txt[["$PAR"]])
		channelNames <- unlist(lapply(1:nChannels,function(i)flowCore:::readFCSgetPar(txt,paste("$P",i,"N",sep="")))) 
		channelNames<- unname(channelNames)
		if(!is.null(channels))#check if channel names contains the specified one 
		{
			channel.notFound<-!is.element(channels,channelNames)
			if(any(channel.notFound))
				stop(channels[channel.notFound],"not Found in ",basename(curFile)," !")
		}else
			channels.all[i]<-paste(channelNames,collapse="|") #save the channel names to list to find common ones later
		
		events[i]<-as.integer(txt[["$TOT"]])
	}
	#try to find common channels among fcs files
	if(is.null(channels))
	{
		
		chnls_unique<-names(table(channels.all))
		chnls_list<-lapply(chnls_unique,function(x){
					strsplit(x,split="\\|")[[1]]
				})
		chnls_common<-chnls_list[[1]]
		
		if(length(chnls_list)==1)
			message("All FCS files have the same following channels:\n"
					,paste(chnls_common,collapse="\n")
					)
		else
		{
			message("Not all FCS files have the same channels\n")
					
			for(i in 2:length(chnls_list))
			{
				chnls_common<-intersect(chnls_common,chnls_list[[i]])
			}	
			message("Only load the following common channels:\n"
					,paste(chnls_common,collapse="\n")
			)
		}
	}
#	bigFile <- files[which.max(file.info(files)[,"size"])]
	tmp <- read.FCSheader(files[1])[[1]]
#	maxEvents <- as.integer(tmp[["$TOT"]])
#	maxChannels <- as.integer(tmp[["$PAR"]])
	maxEvents<-max(events)
	message("Maximum total events: ",maxEvents)
#	channelNames <- unlist(lapply(1:maxChannels,function(i)flowCore:::readFCSgetPar(tmp,paste("$P",i,"N",sep="")))) 
#	channelNames <- unname(channelNames)
	#make a dummy parameters slot for every frames to pass the validity check of flowSet class
	params <- flowCore:::makeFCSparameters(chnls_common,tmp, transformation=F, scale=F,decades=0, realMin=-111)		
	#assign metaData to two environment slots
	e1<-new.env(hash=TRUE, parent=emptyenv())
	e2<-new.env(hash=TRUE, parent=emptyenv())
	for(i in seq_len(nFile)) {
		assign(x=file.names[i]
				,value=new("flowFrame",parameters=params)
				,envir=e1
				) 
		assign(file.names[i],rep(TRUE,maxEvents),e2)
	}
	
	if(!missing(phenoData)){
		guids <- sampleNames(pd)
		pd<-phenoData
	}else{
		guids <- basename(files)
#		
		pd <- AnnotatedDataFrame(data = data.frame(name=guids
											,row.names=guids
											,stringsAsFactors=FALSE
											)
						,varMetadata = data.frame(labelDescription="Name"
													,row.names="name")
						)
						
	}
	
	if(any(duplicated(guids)))
		guids <- make.unique(guids)
	
	#create ncdf ncdf object 
	ncfs<-new("ncdfFlowSet"
				,frames = e1
				, colnames = chnls_common 
				,flowSetId = flowSetId
				, file = ncdfFile
				,maxEvents=maxEvents
				,phenoData=pd
				,indices=e2
				,origSampleVector=guids
				,origColnames=chnls_common
			)

	
	#get the meta size
	metaSize<-ifelse(isSaveMeta,length(serialize(ncfs,NULL)),0)

	#create empty cdf file
#	
	msgCreate <- .Call(dll$createFile, ncdfFile, as.integer(maxEvents), 
					length(chnls_common), as.integer(nFile),
					as.integer(metaSize),as.logical(compress))
	if(!msgCreate)stop()
#	
	##remove indicies to keep the slot as empty by default for memory and speed issue
	initIndices(ncfs,NA)
	#############################################################
	#when isWriteSlice is False,keep the ncdf matrix empty(Mike)
	#############################################################
	if(isWriteSlice)
	{
		lapply(seq_len(nFile), function(i, verbose)
				{
					curFile<-files[i]
					message("writing data:",basename(curFile))
					addFrame(ncfs
							,read.FCS(curFile
										,column.pattern=paste(chnls_common,collapse="|")
										,...)
							,guids[i]
							)
				}, verbose = TRUE)
	}else
	{
		
#		initIndices(ncfs,TRUE)##when isWriteSlice is FALSE,indices have to been explictly initialized
		
	}
#	
	
	
	#write metaData to cdf
	if(isSaveMeta)
		ncdfFlowSet_sync(ncfs)
	
	message("done!")
	return(ncfs)
}

##################################################################
#this function is to be deprecated due to its copy of the entire cdf repository
#if isEmpty is set as FAlSE,then simply copy the orginal cdf file including the data
#by default,create empty cdf file and then add the data later on 
#has the bug of  dimensions (sample*colnames) are not consistent with original cdf
##################################################################
clone.ncdfFlowSet.old<-function(ncfs,newNcFile=NULL,isEmpty=TRUE,isNewNcFile=TRUE,isSaveMeta=FALSE)
{
#		
	
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
			msgCreate <- .Call(dll$createFile, newNcFile, as.integer(ncfs@maxEvents), 
							as.integer(length(colnames(ncfs))), as.integer(length(ncfs)),
							as.integer(metaSize),as.logical(FALSE))
			if(!msgCreate)stop("make sure the file does not exist already or your have write permission to the folder!")
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

clone.ncdfFlowSet<-function(ncfs,ncdfFile=NULL,isEmpty=TRUE,isNew=TRUE,isSaveMeta=FALSE)
{
	
#	
	
	##when isNew==TRUE, the actual data reflected by the current view of ncdfFlowSet is created 
	##and the new ncfs is no longer associated to the orginal one. which mean it is no longer a view of subset
	## of the original cdf repository
	
	
	if(isNew)
	{
		orig<-ncfs#

		#update frames info
		ncfs@frames<-new.env(hash=TRUE, parent=emptyenv())
		for(i in sampleNames(orig))
		{
			assign(i,orig[[i]],ncfs@frames)
		}
		
		if(is.null(ncdfFile))
			ncdfFile<-tempfile(pattern = "ncfs")
		
		if (!length(grep(".", ncdfFile, fixed = TRUE)))  
			ncdfFile <- paste(ncdfFile, "nc", sep = ".")
		#update file info
		ncfs@file<-ncdfFile	
		
		#sync the view info of samplenames and colnames
		ncfs@origSampleVector<-sampleNames(orig)
		ncfs@origColnames<-colnames(orig)
		
		#init view info of indices
		ncfs@indices<-new.env(hash=TRUE, parent=emptyenv())
		for(i in sampleNames(orig))
		{
			updateIndices(ncfs,i,NA)
		}

		
		
		
		metaSize<-ifelse(isSaveMeta,length(serialize(ncfs,NULL)),0)
		msgCreate <- .Call(dll$createFile, ncdfFile, as.integer(ncfs@maxEvents), 
				as.integer(length(colnames(ncfs))), as.integer(length(ncfs)),
				as.integer(metaSize),as.logical(FALSE))
		if(!msgCreate)stop("make sure the file does not exist already or your have write permission to the folder!")
#		
		if(!isEmpty)##write the actual data 
		{
			for(i in sampleNames(orig))
			{
				print(paste("copying data slice:",i))
				addFrame(ncfs,orig[[i]],i)
				
			}
			
		}
	}else
	{
		#when isNew==FALSE, then keep the exact 3D view info of the original ncfs
	
		#copy frames info
		orig<-ncfs@frames
		ncfs@frames<-new.env(hash=TRUE, parent=emptyenv())
		for(i in ls(orig))
		{
			assign(i,orig[[i]],ncfs@frames)
		}
		
		#copy indices info if isEmpty==FALSE
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

			}	
		}
	}
		
	
	ncfs
}

##NOTE:It is for internal usage,do not use this API directly,try to use addFrame always!
.writeSlice <- function(ncfs,data,sampleName)
{ 
	if(class(data)=="flowFrame")
	{
		mat <- exprs(data)
	}else
	{
		mat<-data
	}
	mode(mat) <- "single"
	#make sure to use origSampleVector for IO since phetaData slot may change after subsetting
	i<-which(ncfs@origSampleVector==sampleName)
	############################################################################################
	#make sure the data matrix contains the same channels as original one before write it back
	#since the C function only write the entire ogirinal slice intead of a subset of it
	#and the current iput mat could be a subset through channels or events
	############################################################################################
	
#	
	#if writing the data slice with the exact size and colnames of original one
	#then simply write the input matrix
	#get original slice
	origMat <- .Call(dll$readSlice, ncfs@file, as.integer(c(1,length(ncfs@origColnames))),as.integer(i))
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
			stop("Colnames of the input must be consistent with ncdfFlowSet when adding data to the empty slice of ncdfFlowSet "
					,sampleName)	
		}
		
	}
	#write it back to disk
	msgWrite <- .Call(dll$writeSlice, ncfs@file, origMat , as.integer(i))
	
	if(!msgWrite)
	{
#		
		stop("Writing to CDF file failed!",sampleName)
	}
	msgWrite

}

#to deprecated due to the merging of this routine to [[ for speed
.retNcdfMat <- function(object, chIndx, sampleName,subByIndice){        
	## chIndx is start and end positions and "readSlice" always get a consecutive chunk 
	#to optimize the reading process
	samplePos<-which(object@origSampleVector==sampleName)
#			
	mat <- .Call(dll$readSlice, object@file, as.integer(chIndx),as.integer(samplePos))

	
	if(!is.matrix(mat)&&mat==FALSE) stop()
	
	indx <- seq.int(chIndx[1], chIndx[2], 1)
	##get colnames from frame slot instead of flowSet colnames slot since it may not be synchronized
#	clNames <- colnames(eval(parse(text=paste("object@frames$'",sampleName,"'",sep=""))))[indx]  
	clNames <- object@origColnames[indx]
	colnames(mat) <- clNames
#	
	if(subByIndice&&nrow(mat)>0)
		
		mat[getIndices(object,sampleName),,drop=FALSE]
	else
		mat
	
}

