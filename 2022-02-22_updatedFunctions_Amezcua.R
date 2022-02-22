# Project 1 template
df <- read.table(file, header = FALSE, sep = "\t", quote = "", check.names=FALSE)
# important -- this makes sure our runs are consistent and reproducible
set.seed(0)

file <- "TCGA_breast_cancer_ERpositive_vs_ERnegative_PAM50.tsv"
first10 <- c('NAT1','BIRC5','BAG1','BCL2','BLVRA','CCNB1','CCNE1','CDC6','CDC20','CDH3')
nfold=5

header <- scan(file, nlines = 1, sep="\t", what = character())
data <- read.table(file, skip = 2, header = FALSE, sep = "\t", quote = "", check.names=FALSE)
names(data) <- header

header2 <- scan(file, skip = 1, nlines = 1, sep="\t", what = character())

pos <- data[data$id %in% first10,header2=="Positive"]
neg <- data[data$id %in% first10,header2=="Negative"]


# define function cross_valid so we can rerun the cross validataion with various parameters
#cross_validation <- function (nfold, alg="centroid") {

# split each cancer type samples into nfold groups
pos_groups <- split(sample(colnames(pos)), 1+(seq_along(colnames(pos)) %% nfold))
neg_groups <- split(sample(colnames(neg)), 1+(seq_along(colnames(neg)) %% nfold))

result <- array()

# iterate from 1 to nfold groups -- to choose test group
for (test_group in 1:nfold) {
  
  # return all samples in the chosen test group
  testPos <- pos[,colnames(pos) %in% unlist(pos_groups[test_group])]
  testNeg <- neg[,colnames(neg) %in% unlist(neg_groups[test_group])]
  
  # return all samples *not* in the chosen test group 
  trainingPos <- pos[,!(colnames(pos) %in% unlist(pos_groups[test_group]))]
  trainingNeg <- neg[,!(colnames(neg) %in% unlist(neg_groups[test_group]))]
  
  # compute centroid for each cancer type -- mean for each gene based on all samples
  # note -- rows are gene
  centroidPos <- rowMeans(trainingPos)
  centroidNeg <- rowMeans(trainingNeg)
  
  # For each sample in the test set decide whether it will be classified
  # distance from centroid Lum A: sum(abs(x-centroidLumA))
  # distance from centroid Basal: sum(abs(x-centroidBasal))
  # distance is a sum of distances over all genes 
  # misclassification if when the distance is greater from centroid associated with known result
  misclassifiedPos <- sum(sapply(testPos, function(x) { sum(abs(x-centroidPos))>sum(abs(x-centroidNeg)) }))
  misclassifiedNeg <- sum(sapply(testNeg, function(x) { sum(abs(x-centroidPos))<sum(abs(x-centroidNeg)) }))
  
  result[test_group] <- (misclassifiedPos+misclassifiedNeg)/(ncol(testPos)+ncol(testNeg))
}

c(mean(result), sd(result))
#}

#x<-data.frame(three=cross_validation(3), five=cross_validation(5), ten=cross_validation(10))
#rownames(x) <- c('mean','sd')
#x

