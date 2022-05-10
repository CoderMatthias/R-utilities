#' DESeq plotPCA for all PCs
#'
#' This function expands on the plotPCA function of DESeq2 that only returns the data for the first two principal components. This function will all PCs
#' @param object a DESeqTransform object, with data in assay(x), produced for example by either rlog or varianceStabilizingTransformation.
#' @param intgroup interesting groups: a character vector of names in colData(x) to use for grouping
#' @param ntop number of top genes to use for principal components, selected by highest row variance
#' @keywords DESeq2
#' @export
#' @examples
#' plotPCAmulti()

plotPCAmulti = function(object, intgroup="condition", ntop=500)
{
  # calculate the variance for each gene
  rv <- matrixStats::rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(pca$x, group=group, intgroup.df, name=colnames(object))
  
  # add percent of variance explained by each PC as attribute
  attr(d, "percentVar") <- percentVar
  attr(d, "plotPCAmulti") <- TRUE
  
  return(d)
}


#' Summarize all PCs from plotPCAmulti
#'
#' Function takes the plotPCAmulti dataframe as input and will show all PCs together.
#' @param data Output dataframe from plotPCAmulti
#' @param PCs Number of principal components to include (default: All)
#' @keywords DESeq2
#' @export
#' @examples
#' summarize_plotPCAmulti()

summarize_plotPCAmulti <- function(data, PCs = Inf) {
  
  # Check to see if plotPCAmulti dataframe
  if (is.null(attr(data, "plotPCAmulti"))) {
    stop("This is not an object from the plotPCAmulti function.")
  }
  
  percentVar <- round(100 * attr(data, "percentVar"))
  
  # If user input to number of PCs, shorten percentVar vector
  percentVar <- percentVar[1:min(length(percentVar), PCs)]
  
  GGally::ggpairs(data,
                  upper = list(continuous = "blank"),
                  columns = 1:length(percentVar),
                  columnLabels = paste0(colnames(data[1:length(percentVar)]), " (", percentVar, "%)"),
                  aes(alpha = .5,
                      color = group))
}