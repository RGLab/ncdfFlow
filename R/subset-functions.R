#' subset the ncdfFlowSet/ncdfFlowList based on 'pData'
#' 
#' @param x \code{ncdfFlowSet} or \code{ncdfFlowList}
#' @param subset logical expression(within the context of pData) indicating samples to keep. see \code{\link[base:subset]{subset}}
#' @param ... other arguments. (not used)
#' @return a subset of code{ncdfFlowSet} or \code{ncdfFlowList} object
#' @rdname subset-functions
#' @export 
subset.ncdfFlowSet <- function (x, subset, ...) 
{
  
  pd <- pData(x)
  r <- if (missing(subset)) 
        rep_len(TRUE, nrow(x))
      else {
        e <- substitute(subset)
        r <- eval(e, pd, parent.frame())
        if (!is.logical(r)) 
          stop("'subset' must be logical")
        r & !is.na(r)
      }
  
  x[as.character(rownames(pd[r, ]))]
}


#' @rdname subset-functions
#' @export 
subset.ncdfFlowList <- function (x, subset, ...) 
{
  pd <- pData(x)
  r <- if (missing(subset)) 
        rep_len(TRUE, nrow(x))
      else {
        e <- substitute(subset)
        r <- eval(e, pd, parent.frame())
        if (!is.logical(r)) 
          stop("'subset' must be logical")
        r & !is.na(r)
      }
  
  x[as.character(rownames(pd[r, ]))]
}
