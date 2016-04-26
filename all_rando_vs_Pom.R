########## RANDO ANALYSIS INCLUDING P38 INHIBITOR #########
library(RSkittleBrewer)
library(Biobase)
library(gplots)
library(devtools)
library(dplyr)
library(dendextend)
library(gplots)
library(preprocessCore)
library(qvalue)
library(edgeR)
library(ggplot2)

# Setup the palette
trop = RSkittleBrewer("tropical")
palette(trop)
par(pch=19)

# Loading data into expression set
getwd()
setwd("~/Code/Rando_vs_JD1291/")
list.files()

samples <- c("A1","A2","Q1","Q2","rectus","vastus", "p38_1","p38_2")

# a function to import each .txt file with the "htseq_counts" suffix into R:
read.sample <- function(sample.name) {
  file.name <- paste("./htseq-counts/", sample.name, "_htseq_counts.txt", sep="")
  result <- read.delim(file.name, col.names=c("gene", "count"), sep="\t", colClasses=c("character", "numeric"), header = TRUE)
}

sample.1 <- read.sample(samples[1]) #import the first sample file
head(sample.1)
sample.2 <- read.sample(samples[2])

nrow(sample.1) == nrow(sample.2) #check row number is identical
all(sample.1$gene == sample.2$gene) # are all the values there?

#Now combine all the samples into one table
all.data <- sample.1
all.data <- cbind(sample.1, sample.2$count) #add the second row of sample.2 to the table

# now append the other two sample counts to the table
for (c in 3:length(samples)) {
  temp.data <- read.sample(samples[c])
  all.data <- cbind(all.data, temp.data$count)
}

#We now have a data frame with all the data in it:
colnames(all.data)[2:ncol(all.data)] <- samples #give the columns better names
head(all.data)
dim(all.data)

# Use the first column as row names instead
rownames(all.data) <- all.data$gene
all.data <- all.data[,-1]
head(all.data)

all.matrix <- data.matrix(all.data, rownames.force =TRUE )

# get the pheno data
pData <- read.table("pData_with_p38.txt", row.names=1, header=TRUE, sep="\t")
dim(pData)
pData
names(pData)
sapply(pData, class)

all(rownames(pData)==colnames(all.matrix))

phenoData <- new("AnnotatedDataFrame", data=pData)
annotation <- "hg38"

experimentData <- new("MIAME",name="Joanna Dreux",lab="Pomerantz lab",contact="joanna.dreux@ucsf.edu",
                      title="Muscle satellite cell RNA-Seq", abstract="Comparison of two RNA-Seq experiments",
                      url="https://github.com/Joannacodes/RRrrr-RNA", other=list(
                        notes="Created from text files"))

# let's assemble the ExpressionSet
expSet <- ExpressionSet(assayData=all.matrix,phenoData=phenoData,experimentData=experimentData,
                        annotation=annotation)
head(expSet)

pdata=pData(expSet)
edata=as.data.frame(exprs(expSet))
fdata = fData(expSet)
ls()

##A.--------------EXPLORATORY ANALYSIS--------------------
# Remove rows that are mostly zeroes
edata = as.data.frame(edata)
filt_edata <- edata[rowMeans(edata)>2,]
head(filt_edata)

##B.------------DATA TRANSFORMS + CLUSTERING------------
#filter and log transform
edata_hi = filt_edata[rowMeans(filt_edata) > 500,]
head(edata_hi)
log_edata = log2(edata_hi + 1)
head(log_edata)

# Calculate euclidian distances and cluster the data
dist1 = dist(t(log_edata))
hclust1 = hclust(dist1)
par(mfrow=c(1,1))
plot(hclust1, xlab=NULL)

#use more lenients filtering for normlization
edata_hi <- filt_edata[rowMeans(filt_edata)>100,]
log_edata=log2(edata_hi+1)

##C.---------------QUANTILE NORMALIZATION--------------
#Quantile normalization
norm_edata <- normalize.quantiles(as.matrix(log_edata))
rownames(norm_edata) <- rownames(log_edata)
colnames(norm_edata) <- colnames(log_edata)
head(norm_edata)

#check out the clustering again
dist2 = dist(t(norm_edata))
hclust2 = hclust(dist2)
par(mfrow=c(1,1))
plot(hclust2, xlab=NULL)

##D.-----------DIMENSION REDUCTION FOR GENOMICS----------------
#now PCA analysis to check whether we have any batch effects remaining
svd1 = svd(norm_edata - rowMeans(norm_edata)) #remove the mean variation
plot(svd1$v[,1],svd1$v[,2],xlab="PC1",ylab="PC2",col=as.numeric(pdata$Type),
     main="PC1 vs PC2, all samples", ylim=c(-0.6,1), xlim=c(-0.6,0.8))
text(svd1$v[, 1], svd1$v[, 2], labels = colnames(edata_hi), pos=1)

#Boxplot comparing the PC for different levels of known covariates
par(mfrow=c(1,1))
boxplot(svd1$v[,1] ~ pdata$Type,border=c(1,2), main="PC1 as a function of group(Q/A)",
        ylim=c(-0.6,0.8)) 
points(svd1$v[,1] ~ jitter(as.numeric(pdata$Type)),col=as.numeric(pdata$Type))
text(svd1$v[,1]~ jitter(as.numeric(pdata$Type)), labels = colnames(norm_edata), pos=1)

#Now look at PC 2 - where is it from?
# as a function of lab
par(mfrow=c(1,1))
boxplot(svd1$v[,2] ~ pdata$Lab, border=c(1,2,3), main="PC2 as a function of Lab of origin",
        ylim=c(-0.5,0.9))
points(svd1$v[,2] ~ jitter(as.numeric(pdata$Lab)),col=as.numeric(pdata$Lab))
text(svd1$v[,2] ~ jitter(as.numeric(pdata$Lab)), labels = colnames(norm_edata), pos=1)

# as a function of type Q/A
boxplot(svd1$v[,2] ~ pdata$Type, border=c(1,2,3), main="PC2 as a function of Type Q/A",
        ylim=c(-0.5,0.9))
points(svd1$v[,2] ~ jitter(as.numeric(pdata$Type)),col=as.numeric(pdata$Type))
text(svd1$v[,2] ~ jitter(as.numeric(pdata$Type)), labels = colnames(norm_edata), pos=1)

#as a function of treatment (None/p38i)
boxplot(svd1$v[,2] ~ pdata$Tx, border=c(1,2,3), main="PC2 as a function of Tx",
        ylim=c(-0.5,0.9))
points(svd1$v[,2] ~ jitter(as.numeric(pdata$Tx)),col=as.numeric(pdata$Tx))
text(svd1$v[,2] ~ as.numeric(pdata$Tx), labels = colnames(norm_edata), pos=1)

#both labs and type
boxplot(svd1$v[,2] ~ (pdata$Lab +pdata$Type), border=c(1,2,3), 
        main="PC2 as a function of Lab AND Type", ylim=c(-0.5,0.9))

#what about individual differences?
pdata$indiv <- c(1,2,1,2,3,4,1,2)

boxplot(svd1$v[,2] ~ (pdata$indiv), border=c(1,2,3,4), main="PC2 as a function of individual",
        ylim=c(-0.5,0.9))
points(svd1$v[,2] ~ jitter(as.numeric(pdata$indiv)),
       col=as.numeric(pdata$indiv))
text(svd1$v[,2] ~ as.numeric(pdata$indiv), labels = colnames(norm_edata), pos=1)


##G.--------------------GENE CLUSTERING--------------------------------------------------
#select the interesting genes
my_genes <- c("DLK1","FOSB", "PAX7","MYOD1","RYR3","CHODL","NOTCH1","NOTCH2","LRIG1",
              "CALCR","CHRDL2","TRDN","PTCHD1","NR4A3","RHOJ","CXCR4","EPHB6","NCK2","ABLIM1",
              "EFNA1","EFNB1","EFNB2","NTN4","NFAT5","NFATC1","PLXNA2","SEMA3B","SEMA4A",
              "MYOG","MEF2C","MYH7","MYH3","CD24","SPARCL1","PDK4","VCAN","BUB1","MYBPH",
              "ALPL","ACTC1","F13A1","UNC45B","MYLPF","TNNC2","LUM","CKM","CDKN2A")

plot.me <- norm_edata[my_genes,]
length(my_genes)
dim(plot.me) # yay!

heatmap(plot.me)

#now make a nice heatmap
#breaks for the core of the distribution
breaks=seq(-4, 4, by=0.1)
mycol <- colorpanel(n=length(breaks)-1,low="yellow",mid="orange", high="red")

heatmap.2(as.matrix(plot.me), trace="none", scale="none", col=mycol, density.info = "none",
          symm=F, dendrogram="row", key=T, symkey=F, symbreaks =F)
