
###  estimating sigma^2  for this particular regression


sigma.est <- function(obj, method ="gamma.reg")
{
  
  count <- obj$counts
  group <- obj$group
  x <- model.matrix(~as.factor(group)) # design matrix
  
  normFactor <- estimate.norm.factors(as.matrix(count), method="AH2010")
  nb.data <- prepare.nb.data(as.matrix(count), norm.factor=normFactor)
  
  if (method =="gamma.reg")
  {
    dispersion <- obj$phi
    sigma <- estimate.dispersion.var.edgeR(nb.data, dispersion, x = x)
  }
  
  else if (method== "NBQ")
  {
    
   dispersion <- estimate.dispersion(nb.data = nb.data, x = x, model = "NBQ", method = "MAPL")
   sigma <- estimate.dispersion.var(nb.data, dispersion, x= x) 
  }
  
  else if(method == "edgeR")
    {
    dispersion <- dispersion.edgeR(nb.data$counts, group)
    dispersion1 <- expandAsMatrix(as.vector(dispersion$disp), 
                                  c(dim(count)[1], dim(count)[2]))
    
    sigma <- estimate.dispersion.var.edgeR(nb.data, dispersion1, x = x)
    }
    
    return (sigma)
}