

##' @title Estimate dispersion by Method of Moments (Assuming NB distribution)
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

dispersion.mom <- function(counts, g, r)  

{
  disper <- matrix(0, ncol= g, nrow= dim(counts)[1])
  for ( i in 1:g)
  {
    mean1 = rowMeans(counts[, (1+(i-1)*r):(i*r)])
    var1 <- rowSums((counts[, (1+(i-1)*r):(i*r)]-mean1)^2)/(r-1) 
    disper[, i] <- (var1-mean1)/mean1^2
  }
  disp <- data.frame(disper, row.names= row.names(counts))
  ## return rows with non-negative estimates
  d <- disp[apply(disp,1,function(disp)all(disp>0)),] 
  ## remove NA rows
  d1 <- d[complete.cases(d), ]  
  d2 <- apply(d1, 1, mean)
  d2 <- data.frame(d2)
  return(d2)
}
