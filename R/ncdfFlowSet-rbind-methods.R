# TODO: Add comment
# 
# Author: mike
###############################################################################





setMethod("rbind2",
		signature=signature(x="ncdfFlowSet",
				y="ncdfFlowSet"),
		definition=function(x, y)
		{
			env <- new.env(hash=TRUE, parent=emptyenv())
			env2<-new.env(hash=TRUE, parent=emptyenv())
			lx <- sampleNames(x)
			ly <- sampleNames(y)
			if(any(lx %in% ly))
				stop("These ncdfFlowSets contain overlapping samples.")
			pd1 <- phenoData(x)
			pd2 <- phenoData(y)
			if(!all(varLabels(pd1) == varLabels(pd2)))
				stop("The phenoData of the two frames doesn't match.",
						call.=FALSE)
			
			#new ncdf object 
			#make sure we put the new nc file in the same path as the old ncfile
			newNcFile = ncdfFlow:::.guid()
			newNcFile<-paste(dirname(x@file),newNcFile,sep="/")
			if (!length(grep(".", newNcFile, fixed = TRUE)))  
				newNcFile <- paste(newNcFile, "nc", sep = ".")
			ncfs<-x
			ncfs@file<-newNcFile## assign the new fileid
			ncfs@frames<-env
			ncfs@indices<-env2
			ncfs@maxEvents<- max(x@maxEvents,y@maxEvents)
#			browser()
			
			
			pData(pd1) <- rbind(pData(pd1), pData(pd2))
			ncfs@phenoData <- pd1
			#create new ncdf file		
			#NOTE: rbind2 will not save the metadata in the new file..	
			msgCreate <- try(.Call(dll$createFile, newNcFile, as.integer(ncfs@maxEvents), 
							as.integer(length(colnames(ncfs))), as.integer(length(ncfs)),
							as.integer(0),as.logical(FALSE)),silent = TRUE)
			
			##need to assign the sample vector before add the actual frame 
			ncfs@origSampleVector=c(lx,ly)
#			browser()
			#add frames to env and ncdf file
			for(i in lx)
			{
				assign(i, x@frames[[i]], env=env)
				assign(i, NA, env=env2)
				addFrame(ncfs,x[[i]],i)
			}
			for(i in ly)
			{
							
				
				assign(i, y@frames[[i]], env=env)
				assign(i, NA, env=env2)
				addFrame(ncfs,y[[i]],i)
			}				
			#Comment out sync for metadata
			#ncdfFlowSet_sync(ncfs)
			
			return(ncfs)
		})
##this method is removed because it is pointless to add a flowFrame 
#without specifying the phenoData with the same structure as the existing ncdfFlowSet 
setMethod("rbind2",
		signature=signature(x="ncdfFlowSet",
				y="flowFrame"),
		definition=function(x,y)
		{
			## create dummy phenoData
#			pd <- phenoData(x)[1,]
#			sampleNames(pd)
#			pData(pd)[1,] <- NA
#			browser()
#			ncfs <- as(y,"flowSet")
#			sampleNames(pd) <- sampleNames(tmp) <- "anonymous frame"
#			phenoData(tmp) <- pd
#			rbind2(x, ncfs)
			warning("Please convert the flowFrame to ncdfFlowSet with the appropriate phenoData and then use rbind2 to combine the two ncdfFlowSets!")
			return(NULL)
		})


