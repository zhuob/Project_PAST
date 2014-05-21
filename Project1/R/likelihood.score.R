
##' @title Calculating Likelihood Score and estimating the sigma^2
##' 
##' @param data  read counts matrix
##' @param phi  fitted disperison from regression
##' @param dat1 the combined dataset, used to merge with data
##' @param group a vector, indicating the group ids in \code{data}
##' 
##' @return a likelihood score
##' 
##' @export
##' 
##' @author 
##' 



likelihood.score <- function(data, phi, dat1, group)
{  
  phi.hat <- phi
  arab <- data
  arab <- create.geneColumn(arab)
  arab <- merge(arab, dat1, by= "Gene", all.y=T)[, 1:dim(arab)[2]]
  arab <- data.frame(arab[, -1], row.names= arab$Gene )
  
  
  x <- model.matrix(~as.factor(group)) # design matrix
  phi <- expandAsMatrix(phi.hat, c(dim(arab)[1], dim(arab)[2]))
  
  
  y <- prepare.nb.data(arab)
  y$counts <- as.matrix(y$counts) 
  # otherwise Error: INTEGER() can only be applied to a 'integer', not a 'NULL'
  
  getmu <- irls.nb(y$counts, y$eff.lib.sizes, x,  phi= phi,
                   beta0 = NA, mustart = NULL, maxit = 50,
                   tol.mu = 0.001/length(y), print.level = 0)
  
  likelihood <- ll.nb(kappa = 1/phi, mu=getmu$mu,  y=y$counts)
  
  
  
  obj <- list()
  obj$likelihood <- likelihood
  obj$counts <- arab
  obj$phi <- phi
  obj$group <- group
  return (obj) 
  
}  
