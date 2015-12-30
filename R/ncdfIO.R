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
#' @param mc.cores \code{numeric} passed to parallel::mclapply. Default is NULL, which read FCS files in serial mode.
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
								,ncdfFile
                                ,flowSetId = flowCore:::guid()
								,isWriteSlice= TRUE
								,phenoData
								,channels=NULL
                                ,dim = 2
                                ,compress = 0
                                , mc.cores = NULL
								,...) #dots to be passed to read.FCS
{
    dots <- list(...)
    emptyValue <- ifelse("emptyValue" %in% names(dots), dots[["emptyValue"]], TRUE)
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
	
	
#	browser()
	nFile<-length(files)
    maxEvents <- 0L
    fileIds <- seq_len(nFile)
    
    getChnlEvt <- function(curFile){
      txt <- read.FCSheader(curFile, emptyValue = emptyValue)[[1]]
      nChannels <- as.integer(txt[["$PAR"]])
      channelNames <- unlist(lapply(1:nChannels,function(i)flowCore:::readFCSgetPar(txt,paste("$P",i,"N",sep="")))) 
      channelNames<- unname(channelNames)
      if(!is.null(channels))#check if channel names contains the specified one 
      {
        channel.notFound<-!is.element(channels,channelNames)
        if(any(channel.notFound))
          stop(channels[channel.notFound],"not Found in ",basename(curFile)," !")
      }
      
      c(chnl = paste(channelNames,collapse="|"), evt = txt[["$TOT"]])  
    }
    
    #channel names check
    if(is.null(mc.cores)){
      chnlEvt <- lapply(files, getChnlEvt)  
      
    }else{
      #loadedNamespace may show the package that is loaded but not attached
      if(!any(grepl("parallel", search())))
        require("parallel")
      #parallel at higher granule level to reduce the overhead of dispatching 
      #since each taskgetChnlEvt) is relative small 
      groups <- split(files,ceiling(fileIds/(nFile/mc.cores)))
      
      chnlEvt <- mclapply(groups, function(group){
                                lapply(group, getChnlEvt)
                        }, mc.cores = mc.cores)
        
      chnlEvt <- unlist(chnlEvt, recursive = FALSE)
    }
    
    channels.all <- sapply(chnlEvt, "[[", "chnl")
    events <- as.integer(sapply(chnlEvt, "[[", "evt"))
    
    if(dim == 3)
      maxEvents <- max(events)

    
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

	tmp <- read.FCSheader(files[1], emptyValue = emptyValue)[[1]]
    
    #make a dummy parameters slot for every frames to pass the validity check of flowSet class
	params <- flowCore:::makeFCSparameters(chnls_common,tmp, transformation=F, scale=F,decades=0, realMin=-111)
    
    #determine guids
    if(!missing(phenoData)){
      pd<-phenoData
      guids <- sampleNames(pd)
    }else{
      guids <- basename(files)
      
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
    
    
	#assign metaData to two environment slots
	e1<-new.env(hash=TRUE, parent=emptyenv())
	e2<-new.env(hash=TRUE, parent=emptyenv())
	for(i in seq_len(nFile)) {
		assign(x=guids[i]
				,value=new("flowFrame",parameters=params)
				,envir=e1
				) 
        assign(guids[i], NA, e2)
	}
	
	
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
        
        my.read.FCS <- function(i)
        {
          curFile<-files[i]
          this_fr <- read.FCS(curFile
              ,column.pattern = colPattern
              ,...)
          #we need to reorder columns in order to make them identical across samples
          this_fr[,chnls_common]
        }
        
        if(is.null(mc.cores)){
          #serial read
          for(i in fileIds)
            ncfs[[guids[i], compress = compress]] <- my.read.FCS(i)
            
        }else{
          #parallel read
          #split into groups of files with group size = mc.ores
          groups <- split(fileIds,ceiling(fileIds/mc.cores))
          #read each group with multicore
          for(group in groups){
            
              frList <-  mclapply(group, my.read.FCS, mc.cores = mc.cores)

              #then serial write (because phdf5 is not implemented yet)
              for(j in seq_along(group)){
                
                i <- group[j]
                ncfs[[guids[i], compress = compress]] <- frList[[j]]
              }
                 
              
              
          }
           
          
        }
        
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
                  assign(i, orig[[i]], ncfs@indices)
					
				}else{
					updateIndices(ncfs,y=i,z=orig[[i]])
				}

			}	
		}
	}
		
	
	ncfs
}

#' @title save/load a ncdfFlowSet object to/from disk.
#'
#' @description
#' The \code{ncdfFlowSet} object contains two parts: R object and cdf file.
#' Save/load a ncdfFlowSet mainly involves the R part using saveRDS/readRDS. 
#'
#' @param ncfs A \code{ncdfFlowSet}
#' @param path A character scalar giving the path to save/load the ncdfFlowSet to/from.
#' @param overwrite A logical scalar specifying whether to overwrite the existing folder.
#' @param cdf a character scalar. The valid options are :"copy","move","skip","symlink","link" specifying what to do with the cdf data file.
#'              Sometime it is more efficient to move or create a link of the existing cdf file to the archived folder.
#'
#'
#' @return
#' \code{load_ncfs} returns a ncdfFlowSet object
#'
#' @seealso \code{\link{ncdfFlowSet-class}}
#'
#' @examples
#' \dontrun{
#' 	#ncfs is a ncdfFlowSet
#' 	save_ncfs(fs, path = "tempFolder")
#' 	fs1 <- load_ncfs(path = "tempFolder")
#'
#' }
#' @rdname save_ncfs
#' @export
#' @aliases save_ncfs load_ncfs
save_ncfs <- function(ncfs, path, overwrite = FALSE
    , cdf = c("copy", "move","link", "skip","symlink")){
  
  cdf <- match.arg(cdf)
  
  id <- ncfs@flowSetId
  rds_toSave <- paste(id,"rds",sep=".")
#  browser()
  if(file.exists(path)){
    path <- normalizePath(path,mustWork = TRUE)
    if(overwrite){
      this_files <- list.files(path)
      #validity check for non-empty folder
      if(length(this_files)!=0)
      {
        rds_ind <- grep("\\.rds$",this_files)
        
        if(length(rds_ind)!=1){
          stop("Not a valid ncdfFlowSet archiving folder!")
        }else{
          this_rds <- this_files[rds_ind]

          if(this_rds!=rds_toSave){
            stop("The ncdfFlowSet doesn't match the archived files in: ", path)
          }
        }
      }
      
      #validity check for cdf
      
      if(length(this_files)!=0){
        cdf_ind <- grep("\\.nc$",this_files)
        if(length(cdf_ind) != 1){
          stop("Not a valid ncdfFlowSet archiving folder!")
        }
      }
      
      
      if(length(this_files)!=0)
      {
        #start to delete the old files in path
        file.remove(file.path(path,rds_toSave))
        

        #check if the target path is the same as current cdf path
        this_cdf <- file.path(path,this_files[cdf_ind])
        if(normalizePath(ncfs@file) == this_cdf){
          cdf <- "skip"
        }
        if(cdf != "skip"){
          file.remove(this_cdf)
        }

      }
      
    }else{
      stop(path,"' already exists!")
    }
    
  }else{
    dir.create(path = path)
    #do the dir normalization again after it is created
    path <- normalizePath(path,mustWork = TRUE)
    
  }
  
  
  rds.file <- file.path(path,paste(id,"rds",sep="."))
  
  #save ncdf file
  if(cdf != "skip")
  {
    from <- ncfs@file
  #      browser()
    if(cdf == "move"){
      message("moving ncdf...")
      ncFile <- file.path(path,basename(from))
      res <- file.rename(from,ncFile)
    }else{
      
      ncFile <- tempfile(tmpdir = path, fileext = ".nc")
      
      if(cdf == "copy"){
        message("saving ncdf...")
        res <- file.copy(from=from,to=ncFile)
      }
      else if(cdf == "symlink"){
        message("creating the symbolic link to ncdf...")
        res <- file.symlink(from=from,to=ncFile)
      }else if(cdf == "link"){
        message("creating the hard link to ncdf...")
        res <- file.link(from=from,to=ncFile)
      }
    }
    if(!res){
      stop("failed to ",cdf," ",from,"!")
    }
  }
  
  message("saving R object...")
  saveRDS(ncfs, rds.file)
  
  message("Done\nTo reload it, use 'load_ncfs' function\n")
  
}

#' @rdname save_ncfs
#' @export
#' @aliases load_ncfs
load_ncfs<-function(path){
#  browser()
  path <- normalizePath(path,mustWork = TRUE)
  if(!file.exists(path))
    stop(path,"' not found!")
  files <- list.files(path)
#   browser()
  
  rds.file <- file.path(path,files[grep(".rds$",files)])
  nc.file<-file.path(path,files[grep(".nc$|.nc.trans$",files)])
  if(length(rds.file)==0)
    stop(".rds file missing in ",path)
  if(length(rds.file)>1)
    stop("multiple .rds files found in ",path)
  
  message("loading R object...")
  ncfs <- readRDS(rds.file)
  
  if(length(nc.file)==0)
    stop(".nc file missing in ",path)
  ncfs@file <- nc.file
  
  message("Done")
  ncfs
  
}
