# TODO: Add comment
# 
# Author: mike
###############################################################################
## Dedicated prepanel function to set up dimensions
prepanel.densityplot.ncdfFlowset <- 
		function(x, y, darg=list(n=50), frames, 
				overlap=0.3, subscripts, ...,
				which.channel)
{
	channel.name <- unique(which.channel[subscripts])
	stopifnot(length(channel.name) == 1)
	xl <- range(eapply(frames@frames, range, channel.name), finite=TRUE)
	list(xlim=xl + c(-1,1)*0.07*diff(xl))   
}
setMethod("densityplot",
		signature(x = "formula", data = "ncdfFlowSet"),
		function(x, data, channels, xlab,
				as.table = TRUE, overlap = 0.3,
				prepanel = prepanel.densityplot.ncdfFlowset,
				panel = panel.densityplot.flowset,
				filter=NULL, scales=list(y=list(draw=F)), ...)
		{
#		  browser()
			ocall <- sys.call(sys.parent())
			ccall <- match.call(expand.dots = FALSE)
			ccall <- flowViz:::manipulate.call(ocall, ccall)
			pd <- pData(phenoData(data))
			uniq.name <- flowViz:::createUniqueColumnName(pd)
			## ugly hack to suppress warnings about coercion introducing
			## NAs (needs to be `undone' inside prepanel and panel
			## functions):
			pd[[uniq.name]] <-
					factor(sampleNames(data),
							levels=unique(sampleNames(data))) 
			if (missing(channels))
				channels <-
						setdiff(colnames(data),
								flowCore:::findTimeChannel(data))
			formula.struct <- flowViz:::analyzeDensityFormula(x, dot.names = channels)
			## we want to add a column to pd for each channel, repeating
			## pd as necessary.  We might want to skip this if there is
			## only one channel, but for now we'll use it for
			## conditioning even then.
			
			channel.name <-
					sapply(formula.struct$right.comps, flowViz:::expr2char)
			pd <- rep(list(pd), length(channel.name))
			names(pd) <- channel.name
			pd <- do.call(lattice::make.groups, pd)
			## FIXME: this won't work if pd already has a column named
			## 'which'.  Should deal with that case somehow.
			
			## Next task is to manipulate the formula.  The details of
			## the transformation depends on whether there is a
			## conditioning variable alread.
			## y ~ channel ==> y ~ sample | which
			## y ~ channel | var ==> y ~ sample | which + var
			
			new.x <- d1 ~ d2 | d3
			new.x[[2]] <- ## d1
					if (formula.struct$left) formula.struct$left.symbol
					else as.name("name")
			new.x[[3]][[2]] <- ## d2
					as.name(uniq.name)
			new.x[[3]][[3]] <- ## d3
					if (formula.struct$conditioned) {
						ans <- (~.+.)[[2]]
						ans[[3]] <- as.name("which")
						ans[[2]] <- formula.struct$cond.symbol
						## probably not the ideal order, but I don't see how
						## to easily ake 'which' the first conditioning
						## variable (in case there is more than one
						## conditioning variable to begin with)
						ans
					}
					else as.name("which")
			if (missing(xlab))
				xlab <- ""
			
			gp <- list(...)[["par.settings"]]
			gpar <- flowViz.par.get()
			if(!is.null(gp))
				gpar <- lattice:::updateList(gpar, gp)
			ccall$gpar <- gpar$gate.density
			ccall$x <- new.x
			ccall$data <- pd
			ccall$prepanel <- prepanel
			ccall$panel <- panel
			ccall$frames <- data
			ccall$channel <- formula.struct$right.comps ## channel
			## That is super ugly!!! How do we get to the channel name
			## from the formula???
			ccall$channel.name <- gsub("^.*\\(`|`\\).*$", "", channel.name)
			ccall$as.table <- as.table
			ccall$overlap <- overlap
			ccall$xlab <- xlab
			ccall$horizontal <- TRUE
			ccall$subscripts <- TRUE
			ccall$default.scales <- list(x = list(relation = "free"))
			ccall$which.channel <-
					gsub("^.*\\(`|`\\).*$", "", as.character(pd$which))
			ccall$lattice.options <-
					list(axis.padding =
									list(factor = c(0.6, 1 + 2 * overlap)))
			ccall[[1]] <- quote(lattice::bwplot)
			ans <- eval.parent(ccall)
			ans$call <- ocall
#			browser()
			ans
		})


## xyplot method for flowSets with formula. 
setMethod("xyplot",
		signature=signature(x="formula",
				data="ncdfFlowSet"),
		definition=function(x,
				data,
				smooth=TRUE,
				filter=NULL,
				as.table=TRUE,
				prepanel=flowViz:::prepanel.xyplot.flowset,
				panel=flowViz:::panel.xyplot.flowset,
				xlab=channel.x.name,
				ylab=channel.y.name,
				par.settings=NULL,
				...)
		{
			## no conditioning variable, we chose 'name' as default
			if (length(x[[3]]) == 1){
				tmp <- x[[3]]
				x[[3]] <- (~dummy | name)[[2]]
				x[[3]][[2]] <- tmp
			}
			## par.settings will not be passed on to the panel functions, so
			## we have to fetch it from ... and stick the gate relevant stuff
			## back it in there manually
			gp <- par.settings
			## ugly hack to suppress warnings about coercion introducing
			## NAs (needs to be `undone' inside prepanel and panel
			## functions):
			pd <- pData(data)
			uniq.name <- flowViz:::createUniqueColumnName(pd)
			pd[[uniq.name]] <- factor(sampleNames(data))
			## deparse the formula structure
			channel.y <- x[[2]]
			channel.x <- x[[3]]
			if (length(channel.x) == 3)
			{
				channel.x <- channel.x[[2]]
				x[[3]][[2]] <- as.name(uniq.name)
				x[[2]] <- NULL
			}
			else
			{
				x[[3]] <- as.name(uniq.name)
				x[[2]] <- NULL
			}
			channel.x.name <- flowViz:::expr2char(channel.x)
			channel.y.name <- flowViz:::expr2char(channel.y)
			channel.x <- as.expression(channel.x)
			channel.y <- as.expression(channel.y)
			## use densityplot method with dedicated panel and prepanel
			## functions to do the actual plotting
#		  browser()
			densityplot(x, data=pd, prepanel=prepanel, panel=panel,
					frames=data, channel.x=channel.x,
					channel.y=channel.y, channel.x.name=channel.x.name,
					channel.y.name=channel.y.name, xlab=xlab, ylab=ylab,
					smooth=smooth, gp=gp, as.table=as.table, filter=filter,
					par.settings=par.settings, ...)
		})
