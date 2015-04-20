#override flowSet-version methods to pass data instead of data@frames
#' @aliases densityplot,formula,ncdfFlowSet-method
#' @rdname plot
#' @export 
#' @importMethodsFrom flowViz densityplot xyplot
#' @importFrom flowViz flowViz.par.get panel.densityplot.flowset
setMethod("densityplot",
		signature(x = "formula", data = "ncdfFlowSet"),
		function(x, data, ...)
		{
          
          #construct lattice object
          thisTrellisObj <- flowViz:::.densityplot.adapor(x, data, ...)
          #subset data on channel
          chnl <- thisTrellisObj[["panel.args.common"]][["channel.name"]]
          thisData <- thisTrellisObj[["panel.args.common"]][["frames"]]
          thisData <- thisData[,chnl]
          #update frames
          thisTrellisObj[["panel.args.common"]][["frames"]] <- thisData
          thisTrellisObj
		})

#' @aliases densityplot,formula,ncdfFlowList-method
#' @rdname plot
setMethod("densityplot",
    signature(x = "formula", data = "ncdfFlowList"),
    function(x, data, ...)
    {
      
      selectMethod("densityplot", signature = c("formula", "ncdfFlowSet"))(x, data, ...)
    })

#override flowSet-version methods to pass data instead of data@frames
#' @aliases histogram,formula,ncdfFlowSet-method
#' @rdname plot
#' @export 
#' @importMethodsFrom flowViz histogram xyplot
setMethod("histogram",
    signature(x = "formula", data = "ncdfFlowSet"),
    function(x, data, ...)
    {
      
      #construct lattice object
      thisTrellisObj <- flowViz:::.histogram.adapor(x, data, ...) 
      
      #subset data on channel
      chnl <- thisTrellisObj[["panel.args.common"]][["channel.name"]]
      thisData <- thisTrellisObj[["panel.args.common"]][["frames"]]
      thisData <- thisData[,chnl]
      #update frames
      thisTrellisObj[["panel.args.common"]][["frames"]] <- thisData
      thisTrellisObj
      
    })

#' @aliases histogram,formula,ncdfFlowList-method
#' @rdname plot
setMethod("histogram",
    signature(x = "formula", data = "ncdfFlowList"),
    function(x, data, ...)
    {
      
      selectMethod("histogram", signature = c("formula", "ncdfFlowSet"))(x, data, ...)
    })


#' flowViz plot methods.
#' 
#' @aliases xyplot,formula,ncdfFlowSet-method
#' @rdname plot
#' @param x \code{formula}
#' @param data \code{ncdfFlowSet} or \code{ncdfFlowList}
#' @param ... other arguments passed to \code{flowViz}
#' @export 
setMethod("xyplot",
		signature=signature(x="formula",
				data="ncdfFlowSet"),
		definition=function(x, data, ...)
		{
          #construct lattice object
          thisTrellisObj <- flowViz:::.xyplot.flowSet(x, data, type = "xyplot", ...)
          
          #subset data on channels
          channel.x.name <- thisTrellisObj[["panel.args.common"]][["channel.x.name"]]
          channel.y.name <- thisTrellisObj[["panel.args.common"]][["channel.y.name"]]
          thisData <- thisTrellisObj[["panel.args.common"]][["frames"]]
          thisData <- thisData[,c(channel.x.name,channel.y.name)]
          
          #update frames
          thisTrellisObj[["panel.args.common"]][["frames"]] <- thisData
          
          
          thisTrellisObj
		})

#' @aliases xyplot,formula,ncdfFlowList-method
#' @rdname plot
setMethod("xyplot",
    signature=signature(x="formula",
        data="ncdfFlowList"),
    definition=function(x, data, ...)
    {
      selectMethod("xyplot", signature = c("formula", "ncdfFlowSet"))(x, data, ...)
    })

