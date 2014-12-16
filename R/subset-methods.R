#' subsetting by sampleNames,channels(not for events) methods 
#' 
#' similar to \code{\link[=[,flowSet-method]{[}}.
#'  
#' @param x \code{ncdfFlowSet}
#' @param i sample index(or name)
#' @param j column(or channel) index (or name)
#' @param ... other arguments not used
#' @param drop \code{logical} not used.
#' @rdname subset-methods
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
        #it is too expenstive to call [ on annotationDataFrame 
        #so we simply subset data slot(a data.frame) directly
        # which is sufficient in this context (i.e. no column subsetting involved)
        ncfs@phenoData@data <- ncfs@phenoData@data[i, , drop = FALSE]
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
        
        
        updateIndices(ncfs,y=nm,z=getIndices(x,nm))
#				
        #update channels info for each frame
        if(!missing(j))
        {
#					browser()
          ##get old AnnotatedDataFrame
          pd<-parameters(ncfs@frames[[nm]])
          #update the parameter by subsetting AnnotatedDataFrame wotj parameter name
          if(is.character(j)){
            matchInd <- match(j,pData(pd)$name)
            misMatch <- is.na(matchInd)
            if(any(misMatch)){
              stop("'", paste(j[misMatch], collapse = "','"), "' not found in flow data!")
            }
            parameters(ncfs@frames[[nm]]) <- pd[matchInd,]
          }else{
            parameters(ncfs@frames[[nm]]) <- pd[j,]
          }
          
        }
      }
      
      #update colnames slot for ncdfFlowSet  
      if(!missing(j)){
        if(is.character(j))
          ncfs@colnames <- colnames(x)[match(j, colnames(x))]
        
        else
          ncfs@colnames <- colnames(x)[j]
        
        if(any(is.na(colnames(ncfs))))
          stop("Subset out of bounds")
      }
#			
      return(ncfs)
    })


    
#' @rdname subset-methods
setMethod("[",c(x="ncdfFlowList"),function(x,i,j,...){
      
      if(missing(i) && missing(j)) 
        return(x)
      
      samples <- sampleNames(x)
      
      if(missing(i)){
        sampInd <- NULL
        matchInd <- NULL
      }else{
        
        if(is.numeric(i) || is.logical(i)) {
          sampInd <- sampleNames(x)[i]
        }else
          sampInd <- i
        
        noFound <- is.na(sampInd)
        if(any(noFound)){
          stop("sample index out of boundary!")
        }
        matchInd <- match(sampInd,samples)
        noFound <- is.na(matchInd)
        if(length(matchInd) == 0)
          stop("no sample selected!")
        
        if(any(noFound)){
          stop(paste(i[noFound], collapse = " "), " not found in ", class(x), "!")
        }
      }
      
      if(missing(j))
        j <- NULL
      
      res <- lapply(x,function(object){
            
            this_samples <- sampleNames(object)
            if(is.null(sampInd)){
              this_subset <- this_samples
            }else{
              ind <- match(sampInd,this_samples)
              this_subset <- sampInd[!is.na(ind)]  
            }
            
            if(length(this_subset)>0){
              if(is.null(j))
                return (object[this_subset, ...])
              else
                return (object[this_subset, j, ...])
            }else{
              NULL
            }
          }, level =1)
      res <- res[!unlist(lapply(res,is.null))]
      
      
      if(is.null(matchInd))
        ncdfFlowList(res, x@samples)
      else
        ncdfFlowList(res, names(x@samples[matchInd]))
    })
