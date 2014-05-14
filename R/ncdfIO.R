#' create ncdfFlowSet from FCS files
#' 
#' read FCS files from the disk and load them into a ncdfFlowSet object
#' 
#' @param files A character vector giving the source FCS raw file paths.
#' @param ncdfFile A character scalar giving the output file name. Default is NULL and the function will generate a random
#'                  file in the temporary folder, potentially adding the \code{.cdf} suffix unless a file
#'                  extension is already present. It is sometimes useful to specify this file path to avoid the failure of writing large flow data set to cdf file 
#'                  due to the the shortage of disk space in system temporary folder. 
#'                  It is only valid when isNewNcFile=TRUE 
#' @param flowSetId  A character scalar giving the unique ncdfFlowSet ID.
#' @param isWriteSlice A logical scalar indicating whether the raw data should also be copied.if FALSE,
#'                      an empty cdf file is created with the dimensions (sample*events*channels) supplied by raw FCS files. 
#' @param phenoData An object of \code{AnnotatedDataFrame} providing a way to manually set the phenotyoic data for the whole data set in ncdfFlowSet.
#' @param channels A character vector specifying which channels to extract from FCS files.
#'                  It can be useful when FCS files do not share exactly the same channel names.
#'                  Thus this argument is used to select those common channels that are of interests.
#'                  Default value is NULL and the function will try to scan the FCS headers of all files
#'                  and determine the common channels.
#' @param dim \code{integer} the number of dimensions that specifies the physical storage format of hdf5 dataset.
#'                            Default is 2, which stores each FCS data as a seperate 2d dataset. 
#'                            Normally, user shouldn't need to change this but dim can also be set to 3, which stores all FCS data as one single 3d dataset. 
#' @param compress \code{integer} the HDF5 compression ratio (from 0 to 9). Default is 0, which does not compress the data and is recommeneded (especially for 2d format) because the speed loss usually outweights the disk saving.                              
#' @param ... extra arguments to be passed to \code{\link{read.FCS}}.
#' @return   A ncdfFlowSet object
#' @seealso \code{\link{clone.ncdfFlowSet}}
#' @aliases read.ncdfFlowset
#' @examples 
#' library(ncdfFlow)
#' 
#' path<-system.file("extdata","compdata","data",package="flowCore")
#' files<-list.files(path,full.names=TRUE)[1:3]
#' 
#' #create ncdfFlowSet from fcs with the actual raw data written in cdf
#' nc1  <- read.ncdfFlowSet(files=files,ncdfFile="ncfsTest.nc",flowSetId="fs1",isWriteSlice= TRUE)
#' nc1
#' nc1[[1]]
#' unlink(nc1)
#' rm(nc1)
#' 
#' #create empty ncdfFlowSet from fcs and add data slices afterwards
#' nc1  <- read.ncdfFlowSet(files=files,ncdfFile="ncfsTest.nc",flowSetId="fs1",isWriteSlice= FALSE)
#' fs1<-read.flowSet(files)
#' nc1[[1]] <- fs1[[1]]
#' nc1[[1]]
#' nc1[[2]]
#' @export 
#' @importFrom flowCore read.FCS read.FCSheader
read.ncdfFlowSet <- function(files = NULL
								,ncdfFile,flowSetId=""
								,isWriteSlice= TRUE
								,phenoData
								,channels=NULL
                                ,dim = 2
                                ,compress = 0
								,...) #dots to be passed to read.FCS
{
    dim <- as.integer(match.arg(as.character(dim), c("2","3")))
	#remove nonexisting files
	fileInd<-file.exists(files)
	missingFiles<-files[!fileInd]
	if(length(missingFiles)>0)
		message(paste(missingFiles,"is missing\n"))
	files<-files[fileInd]
	if(length(files)==0)
	{
		message("No file is found!")
		return(NULL)
	}
		
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
    
    maxEvents <- 0L
#	message("Determine the max total events by reading FCS headers ...")
	for(i in 1:nFile){
		curFile<-files[i]
		txt <- read.FCSheader(curFile)[[1]]
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
    if(dim == 3)
      maxEvents <- max(events)
#    message("Maximum total events: ",maxEvents)
    
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

	tmp <- read.FCSheader(files[1])[[1]]
    
	
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
        assign(file.names[i], NA, e2)
#		assign(file.names[i],rep(TRUE,maxEvents),e2)
	}
	
	if(!missing(phenoData)){
		pd<-phenoData
		guids <- sampleNames(pd)
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


	#create empty cdf file
	msgCreate <- .Call(C_ncdfFlow_createFile, ncdfFile, as.integer(maxEvents), 
					length(chnls_common), as.integer(nFile), dim, as.integer(compress))
	if(!msgCreate)stop()
#	
	##remove indicies to keep the slot as empty by default for memory and speed issue
	initIndices(ncfs)
	#############################################################
	#when isWriteSlice is False,keep the ncdf matrix empty(Mike)
	#############################################################
	if(isWriteSlice)
	{
        #escape all the meta characters within channal names
        colPattern <- .escapeMeta(chnls_common)
        colPattern <- paste(colPattern,collapse="|")
		lapply(seq_len(nFile), function(i, verbose)
				{
					curFile<-files[i]
                    this_fr <- read.FCS(curFile
                        ,column.pattern = colPattern
                        ,...)
                    #we need to reorder columns in order to make them identical across samples
                    this_fr <- this_fr[,chnls_common]
					ncfs[[guids[i], compress = compress]] <- this_fr
				}, verbose = TRUE)
	}
	
	
	message("done!")
	return(ncfs)
}
#' escape some special character for a given character vector
#' 
#' @param metaCharacters \code{character} vector specifying the special characters to escape, default is all the mete characters defined in regexprs
#' @param x \code{character} vector the original character vector that contains the special characters
#' @return the modified character vector by padding "\\" before all the special characters 
.escapeMeta <- function(x, metaCharacters = c("\\", "|", "(", ")", "[", "{", "^", "$", "*", "+", "?", ".")){
  
  #pad the \\ before each meta character
  metaCharacters <- paste(sapply(metaCharacters, function(thisMeta)paste0("\\", thisMeta)), collapse = "")
  #wrap into parenthesized pattern 
  metaCharacters <- paste0("([",metaCharacters, "])")
  #escape all the meta characters
  gsub(metaCharacters, "\\\\\\1", x)
}

#' Clone a ncdfFlowSet
#' 
#' Create a new ncdfFlowSet object from an existing one
#' 
#' @param ncfs A \code{\link{ncdfFlowSet}}.
#' @param isNew  A logical scalar indicating whether the new cdf file should be created.
#'   					If FALSE, the original cdf file is associated with the new ncdfFlowSet object.
#' @param ncdfFile A character scalar giving the output file name. By
#'                  default, It is NULL and the function will generate a random
#'                  file name, potentially adding the \code{.cdf} suffix unless a file
#'                  extension is already present. It is only valid when isNewNcFile=TRUE
#' @param isEmpty A logical scalar indicating whether the raw data should also be copied.if FALSE,
#'   				an empty cdf file is created with the same dimensions (sample*events*channels) as the orignial one.
#' @param dim \code{integer} see details in \link{read.ncdfFlowset}.
#' @param compress \code{integer} see details in \link{read.ncdfFlowset}.
#' @return A ncdfFlowSet object
#' @seealso \code{\link{read.ncdfFlowSet}}
#' @examples 
#' 
#' path<-system.file("extdata","compdata","data",package="flowCore")
#' files<-list.files(path,full.names=TRUE)[1:3]
#' 
#' #create ncdfFlowSet from fcs
#' nc1  <- read.ncdfFlowSet(files=files,ncdfFile="ncfsTest.nc",flowSetId="fs1",isWriteSlice= TRUE)
#' 
#' ##clone the ncdfFlowSet object,by default the actual raw data is not added
#' nc2<-clone.ncdfFlowSet(nc1,"clone.nc")
#' nc2[[1]]
#' 
#' #add the actual raw data
#' fs1  <- read.flowSet(files=files)
#' nc2[[sampleNames(fs1)[1]]] <- fs1[[1]]
#' nc2[[1]]
#' 
#' #delete the cdf file associated with ncdfFlowSet before removing it from memory
#' unlink(nc2)
#' rm(nc2)
#' 
#' unlink(nc1)
#' rm(nc1)
#' @export 
clone.ncdfFlowSet<-function(ncfs,ncdfFile=NULL,isEmpty=TRUE,isNew=TRUE, dim = 2, compress = 0)
{
    dim <- as.integer(match.arg(as.character(dim), c("2","3")))
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
			assign(i,orig@frames[[i]],ncfs@frames)
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

		
		msgCreate <- .Call(C_ncdfFlow_createFile, ncdfFile, as.integer(ncfs@maxEvents), 
				as.integer(length(colnames(ncfs))), as.integer(length(ncfs)), dim, as.integer(compress))
		if(!msgCreate)stop("make sure the file does not exist already or your have write permission to the folder!")
#		
		if(!isEmpty)##write the actual data 
		{
			for(i in sampleNames(orig))
			{
				message("copying data slice:",i)
				ncfs[[i, compress = compress]] <- orig[[i]]
				
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



