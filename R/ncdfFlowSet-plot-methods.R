#override flowSet-version methods to pass data instead of data@frames
#' @aliases densityplot,formula,ncdfFlowSet-method
#' @rdname ncdfFlowSet-class
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


#' @aliases xyplot,formula,ncdfFlowSet-method
#' @rdname ncdfFlowSet-class
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
#' @rdname ncdfFlowList-class
setMethod("xyplot",
    signature=signature(x="formula",
        data="ncdfFlowList"),
    definition=function(x, data, ...)
    {
      selectMethod("xyplot", signature = c("formula", "ncdfFlowSet"))(x, data, ...)
    })

#' @aliases densityplot,formula,ncdfFlowList-method
#' @rdname ncdfFlowList-class
setMethod("densityplot",
    signature(x = "formula", data = "ncdfFlowList"),
    function(x, data, ...)
    {
      
      selectMethod("densityplot", signature = c("formula", "ncdfFlowSet"))(x, data, ...)
    })