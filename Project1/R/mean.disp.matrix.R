##' @title Creating a mean and dispersion matrix for a dataset
##' 
##' @param data
##' @param group   a vector of group ids 
##' 
##' @return matrix containing rel.mean freqencies, dispersion and corresponding gene names 
##' 
##' @export
##' 
##' @author 
##' 

pi.disp.matrix <- function(data, group, tol.mu = 1, tol.phi=1e-7)
{
	nf = estimate.norm.factors(data) # norm factors
	nb.dat <- prepare.nb.data(data, norm.factors=nf)

	Gene <- row.names(nb.dat$counts)  # keep the gene names, for convenient merging
	# if counts are not in matrix, error would pop up.
	nb.dat$counts <-  as.matrix(nb.dat$counts)

	x = model.matrix(~ as.factor(group))
 	
	# for each gene, estimate dispersion by MLE  
	d <- fit.nb.glm(nb.dat, x) # NBPSeq 0.3.5
	# mean rel.freq for each gene
	rel.freq <- apply(nb.dat$rel.frequencies, 1, mean)
		
	id <- apply(d$mu>tol.mu, 1, all) & d$phi > tol.phi; 
	# keep genes with dipsersion at some level
	
	disp <- d$phi[id]
	rel.freq <- rel.freq[id]
	Gene <- Gene[id]
	obj <- data.frame(Gene=Gene, dispersion = disp, pi = rel.freq, row.names=NULL)
	return(obj);
}

mean.disp.matrix <- function(data, group) # this one is not used at all.
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

