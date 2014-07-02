setwd("/Users/Bin/Disk1/Stat/Research/Project2014/arabidopsis_data/GSE")

remov <-c("no_feature", "not_aligned", "too_low_aQual","alignment_not_unique",
          "ambiguous", "no_featur")

GSE24198 <- read.table("GSE24198.txt", header = T) # This is way different!
colSums(GSE24198[, -1]) # too low

GSE35288 <- read.table("GSE35288.txt", header = T) # included (arab1)  
tr <- which(GSE35288$Gene %in% remov)
GSE35288 <- GSE35288[-tr, ]
# colSums(GSE35288[,-1])  # too much
# tail(GSE35288)

GSE37159 <- read.table("GSE37159.txt", header = T) # included (arab2)
tr <- which(GSE37159$Gene %in% remov)
GSE37159 <- GSE37159[-tr, ]
# tail(GSE37159)
# colSums(GSE37159[,-1])  # too much


GSE38879 <- read.table("GSE38879.txt", header = T) # included (arab7)
colnames(GSE38879)[1] <- "Gene"
#tr <- which(GSE38879$Gene %in% remov)
#GSE38879 <- GSE38879[-tr, ]
GSE38879$Gene <- substr(GSE38879$Gene, 1, 9)
GSE38879 <- GSE38879[complete.cases(GSE38879),] ## remove NAs
#tail(GSE38879)
# colSums(GSE38879[, -1]) # reasonable

GSE39463 <- read.table("GSE39463.txt", header = T) # to be considered
# colSums(GSE39463[, -1])
head(GSE39463)
dim(GSE39463)

# removed the last two columns, because colSums are too small
GSE42957 <- read.table("GSE42957.txt", header = T) # Included (arab6)
GSE42957 <- GSE42957[, 1:(dim(GSE42957)[2]-2)]
# tr <- which(GSE42957$Gene %in% remov)
# GSE42957 <- GSE42957[-tr, ]
# tail(GSE42957)
# colSums(GSE42957[, -1]) # reasonable

GSE42968 <- read.table("GSE42968.txt", header = T) # included (arab9)
# tr <- which(GSE42968$Gene %in% remov)
# GSE42968 <- GSE42968[-tr, ]
# tail(GSE42968)
colSums(GSE42968[,-1]) # reasonable

# I proceeded this data myself
GSE48767 <- read.table("GSE48767.txt", header = T) # included
GSE48767 <- GSE48767[1:(dim(GSE48767)[1]-5),]
colSums(GSE48767[,-1])

GSE51772 <- read.table("GSE51772.txt", header = T) # included
# head(GSE51772)
# colSums(GSE51772[,-1]) # reasonable


GSE53078 <- read.table("GSE53078.txt", header = T) # included
GSE53078$Gene <- substr(GSE53078$Gene, 1, 9)
GSE53078 <- GSE53078[1:(dim(GSE53078)[1]-5),]
# colSums(GSE53078[, -1])

GSE53952 <- read.table("GSE53952.txt", header = T) # included
GSE53952 <- GSE53952[1:(dim(GSE53952)[1]-5),]
# seems OK, but don't need all columns
## I choose the first stage (col 1, 2-4, 11-13, 20-22 )
col.num <- c(1:4, 11:13, 20:22)
GSE53952 <- GSE53952[, col.num]
colSums(GSE53952[, -1])

GSE56922_10 <- read.table("GSE56922_10.txt", header=T) # contains exons 
tail(GSE56922_10)
dim(GSE56922_10)



library(plyr)
ls <- list(data.frame(GSE35288),data.frame(GSE37159),
           data.frame(GSE38879),data.frame(GSE42957),
           data.frame(GSE42968),data.frame(GSE48767),
           data.frame(GSE51772),data.frame(GSE53078),
           data.frame(GSE53952)
)

stable <- join_all(ls, "Gene")
stable <- stable[complete.cases(stable),]
dim(stable)
as.numeric(colSums(stable[, -1]))



library(NBPSeq)
stable2 <- as.matrix(stable[rowSums(stable[,-1])>= dim(stable)[2], -1])
norm.factors <- estimate.norm.factors(stable2)
names(norm.factors)
nb.data <- prepare.nb.data(stable2, norm.factors=norm.factors)
off.set <- as.numeric(nb.data$eff.lib.sizes)

group <- c(rep(1, 6), rep(2, 8), rep(3, 12), rep(4, 12), rep(5, 6),
           rep(6, 6), rep(7, 8), rep(8, 4), rep(9, 9))
group <- as.factor(group)

trt.1 <- c(rep(1, 3), rep(2, 3))
trt.2 <- c(3, 3, 4, 4, 5, 5, 6, 6)
trt.3 <- c(rep(7:10, each =3))
trt.4 <- c(rep(11:16, each=2))
trt.5 <- c(rep(17:18, each=3))
trt.6 <- c(rep(19:20, each=3))
trt.7 <- rep(21:24, each=2)
trt.8 <- rep(25:26, each=2)
trt.9 <- rep(27:29, each=3)

trt <- c(trt.1, trt.2, trt.3, trt.4, trt.5, trt.6,
         trt.7, trt.8, trt.9)
trt <- as.factor(trt)
