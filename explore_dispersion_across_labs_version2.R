
setwd("/home/zhuob/Project2014/arabidopsis_data/")
library(edgeR)
library(NBPSeq)
#  read a bulk of data simultaneously 
# series <- c(1, 2, 6, 7, 8, 9, 10 ,12)
# datalist= list()
# 
# for ( i in 1: length(series))
# {
#   datalist[[i]] <- read.table(paste("arab", series[i],".txt", sep=""), header=T)
#   
# }


arab1 <- read.table("arab1.txt", header = T) 
## make sure at least one count on average for each row
row.names(arab1) <- arab1[,1]
arab1 <- arab1[, -1]
arab1 <- arab1[ which(rowSums(arab1)> dim(arab1)[2]),] 
arab1 <- arab1[3:(dim(arab1)[1]-3),] 
# delete "unbiguous" and "alignment_not_unique" "no_feature", "no_aligned", "too_low_aQual"

arab2 <- read.table("arab2.txt", header = T)
row.names(arab2) <- arab2[,1]
arab2 <- arab2[, -1]
arab2 <- arab2[ which(rowSums(arab2)> dim(arab2)[2]),] 
arab2 <- arab2[c(-1, -dim(arab2)[1]), ]

## for arab6, if I check the website, then sample 5- 18 are used below
## Col-0 and Col-0 A are compared, and so C24 to C24 A
arab6 <- read.table("arab6.txt", header = T)
row.names(arab6) <- arab6[,1]
arab6 <- arab6[, -1]
arab6 <- arab6[ which(rowSums(arab6)> dim(arab6)[2]),]

arab7 <- read.table("arab7.txt", header = T)
row.names(arab7) <- arab7[,1]
arab7 <- arab7[, -1]
arab7 <- arab7[ which(rowSums(arab7)> dim(arab7)[2]),] 


arab8 <- read.table("arab8.txt", header = T)
row.names(arab8) <- arab8[,1]
arab8 <- arab8[, -1]
arab8 <- arab8[ which(rowSums(arab8)> dim(arab8)[2]),] 


arab9 <- read.table("arab9.txt", header = T)
row.names(arab9) <- arab9[,1]
arab9 <- arab9[, -1]
arab9 <- arab9[ which(rowSums(arab9)> dim(arab9)[2]),] 

arab10 <- read.table("arab10.txt", header = T)
arab10 <- arab10[ which(rowSums(arab10)> dim(arab10)[2]),] 

## multiple columns for a particular gene, eg AT1G01080
arab12 <- read.table("arab12.txt", header = T)
# row.names(arab12) <- arab12[,1]
# arab12 <- arab12[, -1]
arab12 <- arab12[ which(rowSums(arab12[, -1])> dim(arab12)[2]),] 


## Method of moments:  var / mean^2
arab111 <- arab1
dim(arab1)
row.names(arab111) <- arab1[,1]
arab111 <- arab111[, -1]

arab222<- arab2
row.names(arab222) <- arab2[,1]
arab222 <- arab222[,-1]
head(arab222)
head(arab111)


dispersion.mom <- function(data, g, r)  
  # Input: need data, number of group (g), and number of replica per group(r)  
{
  disper <- matrix(0, ncol= g, nrow= dim(data)[1])
  for ( i in 1:g)
  {
    mean1 = rowMeans(data[, (1+(i-1)*r):(i*r)])
    var1 <- rowSums((data[, (1+(i-1)*r):(i*r)]-mean1)^2)/(r-1) 
    disper[, i] <- (var1-mean1)/mean1^2
  }
    disp <- data.frame(disper, row.names= row.names(data))
  ## return rows with non-negative estimates
    d <- disp[apply(disp,1,function(disp)all(disp>=0)),] 
    return(d)
}

disp.lab1 <- dispersion.mom(data= arab1, g= 2, r =3 )
disp.lab2 <- dispersion.mom(data=arab2, g=4, r=2)
disp.lab6 <- dispersion.mom(data=arab6, g=6, r=2)
disp.lab7 <- dispersion.mom(data=arab7, g=4, r=3)
disp.lab8 <- dispersion.mom(data=arab8, g=4, r=3)
disp.lab9 <- dispersion.mom(data=arab9, g=2, r=3)


arab12$mean = apply(as.matrix(arab12[,2:dim(arab12)[2]]), 1,  mean)
arab12$var = apply(as.matrix(arab12[,2:dim(arab12)[2]]), 1,  var)
arab12$var.devideby.squaredmean.lab12 = arab12$var/arab12$mean^2
arab12.disp <- arab12[, c(1, dim(arab12)[2])]

merge1 <- merge(arab1.disp, arab2.disp, by = "Gene")
merge2 <- merge(merge1, arab6.disp, by = "Gene")
merge3 <- merge(merge2, arab7.disp, by = "Gene")
merge4 <- merge(merge3, arab8.disp, by = "Gene")
merge5 <- merge(merge4, arab9.disp, by = "Gene")
disp.7labs <- merge(merge5, arab12.disp, by = "Gene")


# -- -------------------------  Cox-Reid Profile Likelihood ------------------------

dispersion.edgeR <- function(data, group.id)
{
  group <- factor(group.id)
  design <- model.matrix(~group)
  y <- DGEList(counts=data,group=group)
  y <- calcNormFactors(y)
  y <- estimateGLMCommonDisp(y ,design)  # necessary for estimating tagwise dispersion
  y <- estimateGLMTagwiseDisp(y, design)
  disp <- y$tagwise.dispersion
  disp1<- data.frame(disp, row.names= row.names(data))  
  return(disp1)
}

group.id <- c(1, 1, 1, 2, 2, 2)
disp.lab1.cox_reid <- dispersion.edgeR(arab1, group.id)

group <- c(1,1,2,2,3,3, 4, 4)
disp.lab2.cox_reid <- dispersion.edgeR(arab2, group)

group = c(1,1, 2,2,3, 3, 4, 4, 5, 5, 6, 6, 7, 7)
disp.lab6.cox_reid <- dispersion.edgeR(arab6[, 2:15], group)

group = c(1,1, 1,2,2, 2, 3, 3, 3, 4, 4, 4)
disp.lab7.cox_reid <- dispersion.edgeR(arab7, group)

group = c(1,1, 1,2,2, 2, 3, 3, 3, 4, 4, 4)
disp.lab8.cox_reid <- dispersion.edgeR(arab8, group)

group = c(1,1, 1,2,2, 2)
disp.lab9.cox_reid <- dispersion.edgeR(arab9, group)


## group factors see http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42957

arab12.edgeR <- arab12[, 2:7]
group = factor (c(1,1, 2, 2, 3,3))  # ti should actually be 1 rep per group. !!
design <- model.matrix(~group)
y <- DGEList(counts=arab12.edgeR,group=group)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y ,design)  # necessary for estimating tagwise dispersion
y <- estimateGLMTagwiseDisp(y, design)
disp.lab12.cr <- y$tagwise.dispersion
disp.lab12.cox_reid <- data.frame(arab12$Gene, disp.lab12.cr)
colnames(disp.lab12.cox_reid)[1]= "Gene"


merge1 <- merge(disp.lab1.cox_reid, disp.lab2.cox_reid, by = "row.names")
## notice after Merge, the row.names become a new column
colnames(merge1)=c("Gene","disp.lab1", "disp.lab2")
disp.lab6.cox_reid[, 2]= row.names(disp.lab6.cox_reid) 
colnames(disp.lab6.cox_reid)=c( "disp.lab6","Gene")
merge2 <- merge(merge1, disp.lab6.cox_reid, by = "Gene")
head(merge2)

disp.lab7.cox_reid[, 2]= row.names(disp.lab7.cox_reid) 
colnames(disp.lab7.cox_reid)=c( "disp.lab7","Gene")
merge3 <- merge(merge2, disp.lab7.cox_reid, by = "Gene")

head(merge3)
disp.lab8.cox_reid[, 2]= row.names(disp.lab8.cox_reid) 
colnames(disp.lab8.cox_reid)=c( "disp.lab8","Gene")
merge4 <- merge(merge3, disp.lab8.cox_reid, by = "Gene")

disp.lab9.cox_reid[, 2]= row.names(disp.lab9.cox_reid) 
colnames(disp.lab9.cox_reid)=c( "disp.lab9","Gene")
disp.6labs.cox_reid <- merge(merge4, disp.lab9.cox_reid, by = "Gene")

dim(disp.6labs.cox_reid)


#  dispersion matrix estimated by qCML in edgeR

dim(disp.7labs)
dim(disp.7labs.cox_reid)



##################
##  plot the results

lab <- c(1, 2, 6, 7, 8, 9)

n <- 100  ## number of genes to be displayed

set.seed(120)
index <- sample(1:dim(disp.6labs.cox_reid)[1],n)

## create any label I want in x-axis
xval <- seq(1:length(lab))
#-------------------Cox-Reid  profile likelihood ----------
dim(disp.6labs.cox_reid)
plot(xval, log( disp.6labs.cox_reid[1,-1]), xaxt="n", xlab= "lab",
     ylab= "dispersion in log scale", pch = 20, ylim= c(-8,1),
     main=paste("# of genes displayed:", n, "(cox-reid (edgeR))"))
# change x-axis with specified lab labels
axis(1, at=1:length(lab), labels=lab)

lines(xval, log(disp.6labs.cox_reid[1,-1]), col= 6)
for ( i in 1: n)
{
  points(xval, log(disp.6labs.cox_reid[index[i], -1]), pch =20)
  lines(xval, log(disp.6labs.cox_reid[index[i], -1]), col= i) 
}

###
colMeans(disp.6labs.cox_reid[, 2:dim(disp.7labs.cox_reid)[2]])

colSums(arab1[, 2:7])
colSums(arab2[, 2:9])
colSums(arab6[, 2:13])
colSums(arab7[, 2:13])
colSums(arab8[, 2:13])
colSums(arab9[, 2:7])
colSums(arab12[, 2:7])

rank(c(
     mean(colSums(arab1[, 2:7])),
     mean(colSums(arab2[, 2:9])),
     mean(colSums(arab6[, 2:15])),
     mean(colSums(arab7[, 2:13])),
     mean( colSums(arab8[, 2:13])),
     mean(colSums(arab9[, 2:7])),
     mean( colSums(arab12[, 2:7]))  
)
  )

round(cor(log(disp.6labs.cox_reid[, -1])),2)
## pairs(log(disp.6labs.cox_reid[, -1]))
head(disp.6labs.cox_reid)
 
NBPSeq:::smart.plot(log(disp.6labs.cox_reid[,c(4,2)]), clip = 16, pch = 20
                    , xlab="log-dispersion (lab1)", ylab="log-dispersion (lab 6)",
        main= paste( "cor =", round(cor(log(disp.6labs.cox_reid[, c(2, 4)]))[2, 1], 2)))

head(arab6)
arab6$rowsum.lab6 = rowSums(arab6[,1:14])
arab6$Gene <- row.names(arab6)
arab6.rowsum.disp = merge(disp.lab6.cox_reid, arab6[, c(dim(arab6)[2]-1, dim(arab6)[2])], by= "Gene")
head(arab6[, c(1, dim(arab6)[2]-1)])
head(arab6.rowsum.disp)


head(disp.6labs.cox_reid)

round(cor(log(disp.6labs.cox_reid[, c(2, 4)]))[2, 1], 2)
round(cor(log(arab6.rowsum.disp[, 2:3]))[2, 1], 2)
NBPSeq:::smart.plot(log(arab6.rowsum.disp[, 2:3]), clip= 16, pch = 20, 
  main= paste( "both in log scale, cor= ",
               round(cor(log(arab6.rowsum.disp[, 2:3]))[2, 1], 2)))
head(arab6.rowsum.disp)
head(arab6)
sum(arab6[1, 1:14])

