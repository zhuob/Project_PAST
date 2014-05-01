
##' @title Estimate dispersion by Cox-Reid Profile Likelihood (edgeR)
##' 
##' @param counts  read counts matrix
##' @param g  number of groups 
##' @param r  number of replica per group
##' 
##' @return a vector of dispersion for each gene
##' 
##' @export
##' 
##' @author 
##' 


dispersion.edgeR <- function(counts, group.id)
{
  group <- factor(group.id)
  design <- model.matrix(~group)
  y <- DGEList(counts=counts,group=group)
  y <- calcNormFactors(y)
  y <- estimateGLMCommonDisp(y ,design)  # necessary for estimating tagwise dispersion
  y <- estimateGLMTagwiseDisp(y, design)
  disp <- y$tagwise.dispersion
  disp1<- data.frame(disp, row.names= row.names(counts))  
  return(disp1)
}

