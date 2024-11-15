##### Variance Analysis of Gene Expression #####

library(DESeq2)
library(ICC)# to use icc function ICCbare()
library(tidyverse)
library(ggplot2)
library(ggcorrplot )
library(corrplot)
library(limma)

countTable <- read.table(file="data/matrix.featureCountsV2.tsv", header=TRUE, sep="\t", check.names=FALSE, row.names="Geneid")

sample_meta <- read.csv(file="data/AMP_RNAseq_metaData.csv", header=TRUE)
sample_meta <- sample_meta[sample_meta$PCT_CHIMERAS < 0.03,]
sample_meta <- sample_meta[sample_meta$study != 'BioFIND',]
sample_meta <- sample_meta[sample_meta$case_control_other_at_baseline %in% c('Case', 'Control'),]
sample_meta$ageSquared <- sample_meta$age_at_baseline**2
sample_meta$highestScore <- NULL
sample_meta$enrichedScore <- NULL
sample_meta$elevatedScore <- NULL

metaData <- read.csv("data/amp_meta_neutScos.csv", check.names = FALSE)
metaData <- metaData[metaData$sample_id %in% sample_meta$sample_id,] # apply additional filtering

predNeutPer <- read.csv('neutrophil_percent_prediction/knownPredNeutPer_split.csv')
rownames(predNeutPer) <- predNeutPer$X
predNeutPer <- predNeutPer[sample_meta$sample_id,]
sample_meta$neutPer <- predNeutPer$neutPer
rownames(sample_meta) <- sample_meta$sample_id

countSub <- countTable[, sample_meta$sample_id]
all(colnames(countSub) == sample_meta$sample_id) # TRUE

# create vst from featureCounts
dds <- DESeqDataSetFromMatrix(countData = countSub,
                              colData = sample_meta,
                              design = ~1)
vst <- vst(dds, blind = TRUE)

#add in additional PCs to exportated data frame
plotPCA.addPCs <- function (object, ...)
{
  .local <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE)
  {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    } else {
      colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], PC5 = pca$x[, 5], PC6 = pca$x[, 6], PC7 = pca$x[, 7], PC8 = pca$x[, 8], PC9 = pca$x[, 9], group = group, intgroup.df, name = colnames(object))
    if (returnData) {
      attr(d, "percentVar") <- percentVar[1:10]
      return(d)
    }
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
      geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed()
  }
  .local(object, ...)
}


### Without Corrections
pcaData <- plotPCA.addPCs(vst, intgroup='case_control_other_at_baseline', returnData=TRUE, ntop=500)
write.csv(pcaData, 'variance_analysis/pcaData.csv')
attributes(pcaData)$percentVar 
# 0.22860275 0.09146472 0.07876938 0.06485790 0.05060349 0.02625920 0.02416999 0.01725661 0.01288335 0.01111614

var_explained <- data.frame(
  PC=c("PC1","PC2","PC3","PC4","PC5", "PC6", "PC7", "PC8", "PC9") ,  
  perc=c(0.22860275, 0.09146472, 0.07876938, 0.06485790, 0.05060349, 0.02625920, 0.02416999, 0.01725661, 0.01288335)
)

p <- ggplot(var_explained, aes(x=PC, y=perc)) + 
  geom_bar(stat = "identity") +
  ggtitle("Percent Variance Explain - No Correction") +
  ylab("Percent Variance Explained") +
  scale_y_continuous(breaks = seq(0, 0.3, 0.05)) +
  theme_bw()

ggsave('variance_analysis/varianceExplained.pdf', p, width = 7, height = 7)


varianceCor <- function(pca, meta) {
  con_meta <- select_if(meta, is.numeric)
  #con_meta$visit_month <- NULL
  cat_meta <- meta[, !(colnames(meta) %in% colnames(con_meta))]
  con_meta$PC1 <- NULL
  con_meta$PC2 <- NULL
  con_meta$PC3 <- NULL
  con_meta$PC4 <- NULL
  con_meta$PC5 <- NULL
  con_meta$PC6 <- NULL
  con_meta$PC7 <- NULL
  con_meta$PC8 <- NULL
  con_meta$PC9 <- NULL
  
  all(rownames(cat_meta) == rownames(pca))
  icc_values <- data.frame(matrix(nrow=9, ncol=ncol(cat_meta)), row.names = colnames(pca))
  colnames(icc_values) <- colnames(cat_meta)
  for (i in 1:ncol(cat_meta)) {
    x <- numeric(9)
    for (j in 1:9) {
      df <- data.frame(cat_meta[,i], pca[,j])
      x[j]<- ICCbare(colnames(df)[1], colnames(df)[2], df)
    }
    icc_values[,i] <- x
  }
  
  pca_cors_spearman <- cbind(icc_values, cor(pca, con_meta, method = 'spearman'))
  pca_cors_spearman <- as.matrix(pca_cors_spearman^2)
  return(pca_cors_spearman)
}

test <- varianceCor(pcaData[,1:9], sample_meta)
p <- ggcorrplot(t(test), tl.cex = 8, lab=TRUE, lab_size = 2, show.legend = FALSE)
ggsave('variance_analysis/varCorrelation.pdf', p, width = 13, height=4)



variancePvals <- function(pca, meta) {
  con_meta <- select_if(meta, is.numeric)
  #con_meta$visit_month <- NULL
  cat_meta <- meta[, !(colnames(meta) %in% colnames(con_meta))]
  con_meta$PC1 <- NULL
  con_meta$PC2 <- NULL
  con_meta$PC3 <- NULL
  con_meta$PC4 <- NULL
  con_meta$PC5 <- NULL
  con_meta$PC6 <- NULL
  con_meta$PC7 <- NULL
  con_meta$PC8 <- NULL
  con_meta$PC9 <- NULL
  
  all(rownames(cat_meta) == rownames(pca))
  
  anova_p <- data.frame(matrix(nrow=9, ncol=ncol(cat_meta)), row.names = colnames(pca))
  colnames(anova_p) <- colnames(cat_meta)
  for (i in 1:ncol(cat_meta)) {
    if (colnames(cat_meta)[i] %in% c('sample_id', 'participant_id', 'diagnosis_at_baseline', 'diagnosis_latest')) {
      anova_p[,i] <- 1
    } else {
      x <- numeric(9)
      for (j in 1:9) {
        anova <- oneway.test(pca[,j] ~ cat_meta[,i])
        if(anova$p.value == 0) {
          x[j] <- 1e-300
        } else {
          x[j]<- anova$p.value
        }
        
      }
      anova_p[,i] <- x
    }
    
  }
  
  spearman_p <- data.frame(matrix(nrow=9, ncol=ncol(con_meta)), row.names = colnames(pca))
  colnames(spearman_p) <- colnames(con_meta)
  for (i in 1:ncol(con_meta)) {
    x <- numeric(9)
    for (j in 1:9) {
      spear <- cor.test(pca[,j], con_meta[,i], method = 'spearman')
      if(spear$p.value == 0) {
        x[j] <- 1e-300
      } else {
        x[j]<- spear$p.value
      }
      
    }
    spearman_p[,i] <- x
  }
  
  p_values <- cbind(anova_p, spearman_p)
  p_values <- -log10(p_values)
  
  # following columns are null because some values only appear once
  p_values$sample_id <- 0
  p_values$participant_id <- 0
  p_values$diagnosis_at_baseline <- 0
  p_values$diagnosis_latest <- 0
  
  return(p_values)
}

testPvals <- variancePvals(pcaData[,1:9], sample_meta)
p <- ggcorrplot(t(testPvals), tl.cex = 8, lab=TRUE, lab_size = 1.5, show.legend = TRUE)  +  scale_fill_gradient2(limit = c(0,325), low = "white", high =  "red", midpoint = 175)
ggsave('variance_analysis/varCorrelation_pvals.pdf', p, width = 13, height=4)


# With design corrections
vst_corrected <- vst
assay(vst_corrected) <- removeBatchEffect(assay(vst_corrected), 
                                          batch = sample_meta$sex, 
                                          batch2 = sample_meta$case_control_other_at_baseline,
                                          covariates = cbind(sample_meta$neutPer, sample_meta$PCT_MRNA_BASES, sample_meta$ageSquared))
pcaData_corrected <- plotPCA.addPCs(vst_corrected, intgroup='case_control_other_at_baseline', returnData=TRUE, ntop=500)
write.csv(pcaData_corrected, 'variance_analysis/pcaData_corrected.csv')
attributes(pcaData_corrected)$percentVar 
# 0.13068145 0.10262326 0.07091907 0.04130929 0.03175014 0.02692456 0.01650982 0.01500948 0.01455790 0.01264541

var_explained <- data.frame(
  PC=c("PC1","PC2","PC3","PC4","PC5", "PC6", "PC7", "PC8", "PC9") ,  
  perc=c(0.13068145, 0.10262326, 0.07091907, 0.04130929, 0.03175014, 0.02692456, 0.01650982, 0.01500948, 0.01455790)
)

p <- ggplot(var_explained, aes(x=PC, y=perc)) + 
  geom_bar(stat = "identity") +
  ggtitle("Percent Variance Explain - Corrected") +
  ylab("Percent Variance Explained") +
  scale_y_continuous(breaks = seq(0, 0.3, 0.05)) +
  theme_bw()

ggsave('variance_analysis/varianceExplained_corrected.pdf', p, width = 7, height = 7)

varCor <- varianceCor(pcaData_corrected[,1:9], sample_meta)
p <- ggcorrplot(t(varCor), tl.cex = 8, lab=TRUE, lab_size = 2, show.legend = FALSE)
ggsave('variance_analysis/varCorrelation_corrected.pdf', p, width = 13, height=4)


testPvals <- variancePvals(pcaData_corrected[,1:9], sample_meta)
p <- ggcorrplot(t(testPvals), tl.cex = 8, lab=TRUE, lab_size = 1.5, show.legend = TRUE)  +  scale_fill_gradient2(limit = c(0,325), low = "white", high =  "red", midpoint = 175)
ggsave('variance_analysis/varCorrelation_corrected_pvals.pdf', p, width = 13, height=4)


