## Charville (Stem Cell Reports, 2015) ##

library(Biobase)
library(RSkittleBrewer)
library(dendextend)
library(gplots)
library(devtools)
library(dplyr)
library(gplots)
library(preprocessCore)
library(qvalue)
library(edgeR)
library(ggplot2)

getwd()
setwd("~/Code/Rando_ex_vivo/")
list.files()

# Setup the palette
trop = RSkittleBrewer("tropical")
palette(trop)
par(pch=19)

# Setup expression set
samples <- c("A1","A2","Q1","Q2","p38_1","p38_2")

# a function to import each .txt file with the "htseq_counts" suffix into R:
read.sample <- function(sample.name) {
  file.name <- paste("./htseq_counts/", sample.name, "_htseq_counts.txt", sep="")
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
pData <- read.table("pData_p38.txt", row.names=1, header=TRUE, sep="\t")
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

##A.---------------------- EXPLORATORY ANALYSIS --------------------------------

# Check genomic data for missing data
sum(is.na(edata))

#Make sure dimensions match up
dim(fdata)
dim(pdata)
dim(edata)

# Calculate euclidian distances and cluster the data
edata = as.data.frame(edata)
log_edata <- log2(edata[rowMeans(edata)>500,]+1)
dim(log_edata) # down to 11k
dist1 = dist(t(log_edata))
hclust1 = hclust(dist1) # dendogram plot
par(mfrow=c(1,1))
plot(hclust1, xlab=NULL)

#Look at overall distributions - see a lot of outliers
par(mfrow=c(1,1))
boxplot(log2(edata+1),col="darkorange", main="raw distribution")

#We can also look at this sample by sample with histograms
par(mfrow=c(2,3))
hist(log2(edata[,1]+1),col="darkorange", main="Activated 2", ylim=c(0, 30000),
     xlab="Gene exp distribution (log2)")
hist(log2(edata[,2]+1),col="darkorange", main="Activated 2", ylim=c(0, 30000),
     xlab="Gene exp distribution (log2)")
hist(log2(edata[,3]+1),col="darkorange", main="Quiescent 1", ylim=c(0, 30000),
     xlab="Gene exp distribution (log2)")
hist(log2(edata[,4]+1),col="darkorange", main="Quiescent 2", ylim=c(0, 30000),
     xlab="Gene exp distribution (log2)")
hist(log2(edata[,5]+1),col="darkorange", main="p38 inhibitor 1", ylim=c(0, 30000),
     xlab="Gene exp distribution (log2)")
hist(log2(edata[,5]+1),col="darkorange", main="p38 inhibitor 2", ylim=c(0, 30000),
     xlab="Gene exp distribution (log2)")

##B.------------DATA TRANSFORMS + CLUSTERING------------

## FILTERING 
#look at all genes with a CPM greater than 2 and lower than 5000
filt_edata <- edata
cpm_edata <- cpm(filt_edata)
summary(cpm_edata)
dim(cpm_edata) #38K

# keep features > 1CPM in at least 2 samples
countCheck <- cpm_edata > 2
head(countCheck)
keep <- which(rowSums(countCheck) >= 2) 
filt_edata <- filt_edata[keep,]
dim(filt_edata) # down to 13K

# keep features  < 5000 CPM in at least 2 samples
cpm_edata <- cpm(filt_edata) # recompute cpm 
countCheck <- cpm_edata < 5000
keep <- which(rowSums(countCheck) >= 2)
filt_edata <- filt_edata[keep,]
summary(filt_edata)
dim(filt_edata)

# boxplot following filtering
par(mfrow=c(1,1))
boxplot(as.matrix(log2(filt_edata+1)),col="dodgerblue", main="After removing low expression genes")

# with histograms again
par(mfrow=c(2,3))
hist(log2(filt_edata[,1]+1),col="dodgerblue", main="Activated 2", ylim=c(0, 4000),
     xlab="Gene exp distribution (log2)")
hist(log2(filt_edata[,2]+1),col="dodgerblue", main="Activated 2", ylim=c(0, 4000),
     xlab="Gene exp distribution (log2)")
hist(log2(filt_edata[,3]+1),col="dodgerblue", main="Quiescent 1", ylim=c(0, 4000),
     xlab="Gene exp distribution (log2)")
hist(log2(filt_edata[,4]+1),col="dodgerblue", main="Quiescent 2", ylim=c(0, 4000),
     xlab="Gene exp distribution (log2)")
hist(log2(filt_edata[,5]+1),col="dodgerblue", main="p38 inhibitor 1", ylim=c(0, 4000),
     xlab="Gene exp distribution (log2)")
hist(log2(filt_edata[,6]+1),col="dodgerblue", main="p38 inhibitor 2", ylim=c(0, 4000),
     xlab="Gene exp distribution (log2)")


##C.---------------QUANTILE NORMALIZATION--------------

norm_edata <- normalize.quantiles(as.matrix(log2(filt_edata+1)))
rownames(norm_edata) <- rownames(filt_edata)
colnames(norm_edata) <- colnames(filt_edata)
head(norm_edata)

par(mfrow=c(2,3))
hist(norm_edata[,1],col="hotpink", main="A1", ylim=c(0, 2000),xlab=NULL, xlim=c(0,25))
hist(norm_edata[,2],col="hotpink", main="A2", ylim=c(0, 2000),xlab=NULL,xlim=c(0,25))
hist(norm_edata[,3],col="hotpink", main="Q1", ylim=c(0, 2000),xlab=NULL,xlim=c(0,25))
hist(norm_edata[,4],col="hotpink", main="Q2", ylim=c(0, 2000),xlab=NULL,xlim=c(0,25))
hist(norm_edata[,5],col="hotpink", main="p38i_1", ylim=c(0, 2000),xlab=NULL,xlim=c(0,25))
hist(norm_edata[,6],col="hotpink", main="p38i_2", ylim=c(0, 2000),xlab=NULL,xlim=c(0,25))
#yay!


##D.-----------DIMENSION REDUCTION ----------------

svd1 = svd(norm_edata - rowMeans(norm_edata)) #remove the mean variation
par(mfrow=c(1,1))
plot(svd1$v[,1],svd1$v[,2],xlab="PC1",ylab="PC2",col=c(1,1,2,2,3,3),
     main="PC1 vs PC2", ylim=c(-0.7,0.6))
legend("top", legend=c("Activated", "Quiescent", "Treated"), col=c(1,2,3), pch=19)
text(svd1$v[, 1], svd1$v[, 2], labels = colnames(norm_edata), pos=1)

#pc1 ~ study
boxplot(svd1$v[,1] ~ pdata$Type,border=c(1,2), main="PC1 as a function of group") 
points(svd1$v[,1] ~ jitter(as.numeric(pdata$Type)),col=as.numeric(pdata$Type))

#pc2 ~ donor
boxplot(svd1$v[,2] ~ pdata$Donor,border=c(1,2), main="PC2 as a function of Donor")
points(svd1$v[,2] ~ jitter(as.numeric(pdata$Donor)),col=as.numeric(pdata$Donor))


##E.--------------------------------------EDGER-------------------------------------------
pdata
head(filt_edata)

#define 3 groups
group <- factor(paste(pdata$Type,pdata$Tx, sep="."))
pdata$group <- group
pdata$group <- relevel(pdata$group, ref="Act.None")

#design matrix
design <- model.matrix(~Donor+group, data=pdata)
design

#create dge object
cds1 <- DGEList(filt_edata)
dim(cds1)
names(cds1)
head(cds1$counts)
cds1$samples # a summary of the samples

#estimate disp
cds1<- estimateDisp(cds1, design=design)
cds1 <- calcNormFactors(cds1)

#fit a glm
fit <- glmFit(cds1,design)
colnames(fit)
##[1] "(Intercept)"     "Donor"    "groupAct.p38i"   "groupQuies.None"

# lrt <- glmLRT(fit) Quies - Act
# lrt <- glmLRT(fit, coef=3) p38i - Act

##--------------- Q vs A --------------------------------------
lrt <- glmLRT(fit)
summary(de <- decideTestsDGE(lrt))
detags <- rownames(cds1)[as.logical(de)]
plotSmear(lrt, de.tags=detags, main="DE genes between Q and A cells")
abline(h=c(-1, 1), col="blue")

##look at CPMs for top DE genes
o <- order(lrt$table$PValue)
top_de_cpm <- cpm(cds1)[o[1:10],1:4]
write.csv(top_de_cpm,file="top_de_cpm.csv", eol="\r",row.names = TRUE)

##Volcano plot
lrt.Tbl <- lrt$table 
de.genes <- lrt$table[as.logical(de),]
sum(abs(de.genes$logFC)>2) #the number of genes above this cutoff

##Save results in a table
write.csv(de.genes, file="QvsA_GLM.csv", eol="\r", row.names = TRUE)

# BH correction
lrt.Tbl$PValue_adj = p.adjust(lrt.Tbl$PValue, method='BH')
#set threshold for significance
lrt.Tbl$threshold = as.factor(abs(lrt.Tbl$logFC)>2 & lrt.Tbl$PValue_adj<0.05)
table(lrt.Tbl$threshold) # consistent with DE genes!
#Construct the plot
ggplot(data=lrt.Tbl, aes(x=logFC, y=-log10(PValue), colour=threshold)) + 
  geom_point(alpha=0.4, size=1.75) + xlab("log2 fold change") + ylab("-log10 p-value")

## -----------------------Tx vs no Tx in Act----------------------------
colnames(fit)
lrt_Tx <- glmLRT(fit, coef=3)

## get results
summary(de_Tx <- decideTestsDGE(lrt_Tx))
detags_Tx <- rownames(cds1)[as.logical(de_Tx)]
plotSmear(lrt_Tx, de.tags=detags_Tx, main="DE genes between A and A+p38i cells")
abline(h=c(-1, 1), col="blue")

lrt_Tx.Tbl <- lrt_Tx$table 
de.genes_Tx <- lrt_Tx$table[as.logical(de_Tx),]
sum(abs(de.genes_Tx$logFC)>2) #the number of genes above this cutoff

##Save results in a table
write.csv(de.genes_Tx, file="p38vsA_GLM.csv", eol="\r", row.names = TRUE)

##look at CPMs for top DE genes
o_Tx <- order(lrt_Tx$table$PValue)
top_de_cpm_Tx <- cpm(cds1)[o[1:10],]
write.csv(top_de_cpm_Tx,file="top_de_cpm_Tx.csv", eol="\r",row.names = TRUE)