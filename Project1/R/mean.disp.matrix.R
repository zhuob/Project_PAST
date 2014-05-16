##' @title Creating a mean and dispersion matrix for a dataset
##' 
##' @param data
##' @param group   a vector of group ids 
##' 
##' @return matrix containing rel.mean freqencies and dispersion
##' 
##' @export
##' 
##' @author 
##' 

mean.disp.matrix <- function(data, group)
{
  # get dispersion using edgeR cox-reid profile likelihood, accounting for design structure
  disp <- dispersion.edgeR(data, group)
  # create a column that contains gene info to facilitate merging by gene
  disp1 <- create.geneColumn(disp)
  ## keep colunms with positive disp
  disp2 <- disp1[which(disp1[, 2]>0),]
  
  # get rel.mean frequencies from NBPSeq
  nf <- estimate.norm.factors(data)	
  y <- prepare.nb.data(data, norm.factors= nf)
  # average of rel. mean across samples
  rel.freq <- data.frame(rel.freq=apply(y$rel.frequencies, 1, mean))
  rel.freq1 <- create.geneColumn(rel.freq)
  ## keep colunms with positive mean rel.freq
  rel.freq2 <- rel.freq1[which(rel.freq1[,2]>0), ]
  
  # combining the two datasets
  target <- merge(disp2, rel.freq2, by = colnames(disp2)[1])
  return(target)
}

