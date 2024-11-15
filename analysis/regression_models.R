##### Developing Linear Models for Neutrophil % Prediction ######

# load libraries
library(DESeq2)
library(data.table)
library(Metrics)

# load data
countTable <- read.table(file="data/matrix.featureCountsV2.tsv", header=TRUE, sep="\t", check.names=FALSE, row.names="Geneid")
metaData <- read.csv("data/amp_meta_neutScos.csv", check.names = FALSE)
sample_meta <- read.csv(file="data/AMP_RNAseq_metaData.csv", header=TRUE)
sample_meta <- sample_meta[sample_meta$PCT_CHIMERAS < 0.03,]
sample_meta <- sample_meta[sample_meta$study != 'BioFIND',]
metaData <- metaData[metaData$sample_id %in% sample_meta$sample_id,] # apply additional filtering
write.csv(metaData, 'data/amp_meta_neutScos_filtered.csv')

countSub <- countTable[, metaData$sample_id]
all(colnames(countSub) == metaData$sample_id) # TRUE


# load in blood cell enriched genes
# import 355 neutrophil enriched genes
neut <- fread("data/neutrophil_lineageEnriched.tsv", select = c("Ensembl", "Gene"))

# import 225 basophil enriched genes
baso <- fread("data/basophil_lineageEnriched.tsv", select = c("Ensembl", "Gene"))

# import 106 eosinophil enriched genes
eosino <- fread("data/eosinophil_lineageEnriched.tsv", select = c("Ensembl", "Gene"))

# import 189 monocyte lineage enriched genes
mono <- fread("data/monocytes_lineageEnriched.tsv", select = c("Ensembl", "Gene"))

# import 609 t-cell lineage enriched genes (lymphocytes)
t_cells <- fread("data/t-cells_lineageEnriched.tsv", select = c("Ensembl", "Gene"))

# import 242 b-cell lineage enriched genes (lymphocytes)
b_cells <- fread("data/b-cells_lineageEnriched.tsv", select = c("Ensembl", "Gene"))

# import 347 dendritic lineage enriched genes
den <- fread("data/dendritic_lineageEnriched.tsv", select = c("Ensembl", "Gene"))

set.seed(777)
patno <- unique(metaData$participant_id)
sample_size <- floor(0.80*length(patno))

patID <- sample(patno,size = sample_size)
picked <- metaData[metaData$participant_id %in% patID,]$sample_id
development <- countSub[,colnames(countSub) %in% picked]
holdout <- countSub[, !(colnames(countSub) %in% picked)]
metaDev <- metaData[colnames(countSub) %in% picked,]
rownames(metaDev) <- metaDev$sample_id
metaHold <- metaData[!(colnames(countSub) %in% picked),]
rownames(metaHold) <- metaHold$sample_id
write.csv(metaDev, 'data/metaDev')

# create vst  for development
dds <- DESeqDataSetFromMatrix(countData = development,
                              colData = metaDev,
                              design = ~1)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
vstDev <- as.data.frame(t(assay(vsd)))
write.csv(vstDev, 'data/vstCounts_dev.csv')


# create vst for holdout
dds <- DESeqDataSetFromMatrix(countData = holdout,
                              colData = metaHold,
                              design = ~1)
dds <- estimateSizeFactors(dds)
normCounts <- as.data.frame(counts(dds, normalize = TRUE))
vsd <- vst(dds, blind = TRUE)
vstHoldout <- as.data.frame(t(assay(vsd)))
write.csv(vstHoldout, 'data/vstCounts_holdout.csv')


colnames(vstDev) <- sub('\\.[0-9]*$', '', colnames(vstDev))
vstDev <- vstDev[, colnames(vstDev) %in% c(neut$Ensembl, baso$Ensembl, eosino$Ensembl, mono$Ensembl, t_cells$Ensembl, b_cells$Ensembl, den$Ensembl)] 
all(rownames(vstDev) == metaDev$sample_id) #TRUE

colnames(vstHoldout) <- sub('\\.[0-9]*$', '', colnames(vstHoldout))
vstHoldout <- vstHoldout[,colnames(vstHoldout) %in% c(neut$Ensembl, baso$Ensembl, eosino$Ensembl, mono$Ensembl, t_cells$Ensembl, b_cells$Ensembl, den$Ensembl)] 
all(rownames(vstHoldout) == metaHold$sample_id) #TRUE

##### Blood Cell Based Model

# Neutrophil only
vstSub <- vstDev[, colnames(vstDev) %in% neut$Ensembl]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.8483
neut_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8741217
sigGenes <- neut_coef[neut_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.8039
neut_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8911026
sigGenes <- neut_coef[neut_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.8032
neut_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8906584
sigGenes <- neut_coef[neut_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.8029
neut_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8902983


# Basophil only
vstSub <- vstDev[, colnames(vstDev) %in% baso$Ensembl]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.7769
baso_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8239834
sigGenes <- baso_coef[baso_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.7428
baso_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8261513
sigGenes <- baso_coef[baso_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.7415
baso_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8286698
sigGenes <- baso_coef[baso_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.7409
baso_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8272422


# Eosinophil only
vstSub <- vstDev[, colnames(vstDev) %in% eosino$Ensembl]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.6969
eosino_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8030545
sigGenes <- eosino_coef[eosino_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.6909
eosino_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.7871803


# Monoctyes only
vstSub <- vstDev[, colnames(vstDev) %in% mono$Ensembl]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.74
mono_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8009005
sigGenes <- mono_coef[mono_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.7233
mono_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8199335
sigGenes <- mono_coef[mono_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.7223
mono_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.824762


# Lymphocytes only
vstSub <- vstDev[, colnames(vstDev) %in% c(t_cells$Ensembl, b_cells$Ensembl)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.8232
lymph_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.5815887
sigGenes <- lymph_coef[lymph_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.6562
lymph_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8107598
sigGenes <- lymph_coef[lymph_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.6559
lymph_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.81855
sigGenes <- lymph_coef[lymph_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.655
lymph_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8209289



# Dendritic cells only
vstSub <- vstDev[, colnames(vstDev) %in% den$Ensembl]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.7694
den_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8073885
sigGenes <- den_coef[den_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.6709
den_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.7645713
sigGenes <- den_coef[den_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.6687
den_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.7689629


# with all significant genes
allSigGenes <- c(rownames(neut_coef), rownames(baso_coef), rownames(eosino_coef), rownames(mono_coef), rownames(lymph_coef), rownames(den_coef))

vstSub <- vstDev[, colnames(vstDev) %in% allSigGenes]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.8466
all_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8971746
sigGenes <- all_coef[all_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.8263
all_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8935355
sigGenes <- all_coef[all_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.8258
all_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8907272

write.csv(all_coef, 'neutrophil_percent_prediction/blood_cell_model_genes_split.csv')


##### MI Based Model
mi_genes <- colnames(read.delim("neutrophil_percent_prediction/mi_genes.txt", sep = ","))
mi_genes <- sub('\\.[0-9]*$', '', mi_genes)

vstSub <- vstDev[, colnames(vstDev) %in% mi_genes]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.7976
mi_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8985171
sigGenes <- mi_coef[mi_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.7903
mi_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8923074
sigGenes <- mi_coef[mi_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.7891
mi_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.8927676

write.csv(mi_coef, 'neutrophil_percent_prediction/mi_model_genes_split.csv')



##### Combined Model

vstSub <- vstDev[, colnames(vstDev) %in% c(rownames(mi_coef), rownames(all_coef))]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.8382
comb_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.9041311
sigGenes <- comb_coef[comb_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.8358
comb_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.9033218
sigGenes <- comb_coef[comb_coef$coefficients.Pr...t.. < 0.05,]

vstSub <- vstDev[, colnames(vstDev) %in% rownames(sigGenes)]
vstSub$neutPer <- metaDev$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstSub)
summary(fit) # 0.8349
comb_coef <- as.data.frame(summary(fit)[4])
cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`) #0.9034231

write.csv(comb_coef, 'neutrophil_percent_prediction/combined_model_genes_split.csv')

##### XGBoost Model
library(xgboost)

test <- xgboost(data = as.matrix(vstDev), label = metaDev$`Neutrophils (%) (%)`, max.depth = 3, eta = 0.3, nrounds =10)

pred <- predict(test, as.matrix(vstHoldout))
cor(pred, metaHold$`Neutrophils (%) (%)`)**2


xgboost_eval <- function(genes) {
  results <- data.frame(matrix(nrow=100, ncol=3))
  colnames(results) <- c('R2', 'RMSE', 'MAE')
  sample_size <- floor(0.80*length(patno))
  for (i in 1:100) {
    patID <- sample(patno,size = sample_size)
    picked <- metaData[metaData$participant_id %in% patID,]$sample_id
    development <- countSub[,colnames(countSub) %in% picked]
    holdout <- countSub[, !(colnames(countSub) %in% picked)]
    metaDev <- metaData[colnames(countSub) %in% picked,]
    rownames(metaDev) <- metaDev$sample_id
    metaHold <- metaData[!(colnames(countSub) %in% picked),]
    rownames(metaHold) <- metaHold$sample_id
    
    # create vst  for development
    dds <- DESeqDataSetFromMatrix(countData = development,
                                  colData = metaDev,
                                  design = ~1)
    dds <- estimateSizeFactors(dds)
    normCounts <- as.data.frame(counts(dds, normalize = TRUE))
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    vstDev <-t(assay(vsd))
    #write.csv(vstCounts, 'data/vstCounts_updated.csv')
    
    
    # create vst for holdout
    dds <- DESeqDataSetFromMatrix(countData = holdout,
                                  colData = metaHold,
                                  design = ~1)
    dds <- estimateSizeFactors(dds)
    normCounts <- as.data.frame(counts(dds, normalize = TRUE))
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    vstHoldout <- t(assay(vsd))
    
    test <- xgboost(data = vstDev, label = metaDev$`Neutrophils (%) (%)`, max.depth = 3, eta = 0.3, nrounds = 10)
    
    results$R2[i] <- cor(predict(test, vstHoldout), metaHold$`Neutrophils (%) (%)`)**2 
    results$RMSE[i] <- rmse(metaHold$`Neutrophils (%) (%)`, predict(test, vstHoldout))
    results$MAE[i] <- mae(metaHold$`Neutrophils (%) (%)`, predict(test, vstHoldout))
  }
  return(results)
}

rownames(countSub) <- sub('\\.[0-9]*$', '', rownames(countSub))
rownames(metaData) <- metaData$sample_id
all(rownames(metaData) ==colnames(countSub))

xgb_eval <- xgboost_eval(countSub)

write.csv(xgb_eval, 'neutrophil_percent_prediction/xgb_eval_split.csv')

xgb_eval <- read.csv('neutrophil_percent_prediction/xgb_eval_split.csv')
xgb_eval$X <- NULL

##### Model Evaluation

all_coef <- read.csv('neutrophil_percent_prediction/blood_cell_model_genes_split.csv')
mi_coef <- read.csv('neutrophil_percent_prediction/mi_model_genes_split.csv')
comb_coef <- read.csv('neutrophil_percent_prediction/combined_model_genes_split.csv')

rownames(countSub) <- sub('\\.[0-9]*$', '', rownames(countSub))
countSub_bio <- countSub[rownames(countSub) %in% all_coef$X,]
countSub_mi <- countSub[rownames(countSub) %in% mi_coef$X,]
countSub_comb <- countSub[rownames(countSub) %in% comb_coef$X,]

sample_size <- floor(0.80*length(patno))
model_eval <- function(genes) {
  results <- data.frame(matrix(nrow=100, ncol=3))
  colnames(results) <- c('R2', 'RMSE', 'MAE')
  for (i in 1:100) {
    patID <- sample(patno,size = sample_size)
    picked <- metaData[metaData$participant_id %in% patID,]$sample_id
    development <- genes[,colnames(genes) %in% picked]
    holdout <- genes[, !(colnames(genes) %in% picked)]
    metaDev <- metaData[colnames(genes) %in% picked,]
    rownames(metaDev) <- metaDev$sample_id
    metaHold <- metaData[!(colnames(genes) %in% picked),]
    rownames(metaHold) <- metaHold$sample_id
    
    # create vst  for development
    dds <- DESeqDataSetFromMatrix(countData = development,
                                  colData = metaDev,
                                  design = ~1)
    dds <- estimateSizeFactors(dds)
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    vstDev <- as.data.frame(t(assay(vsd)))
    
    
    # create vst for holdout
    dds <- DESeqDataSetFromMatrix(countData = holdout,
                                  colData = metaHold,
                                  design = ~1)
    dds <- estimateSizeFactors(dds)
    vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    vstHoldout <- as.data.frame(t(assay(vsd)))
    
    vstDev$neutPer <- metaDev$`Neutrophils (%) (%)`
    fit <- lm(neutPer~., data=vstDev)
    results$R2[i] <- cor(predict(fit, vstHoldout), metaHold$`Neutrophils (%) (%)`)**2 
    results$RMSE[i] <- rmse(metaHold$`Neutrophils (%) (%)`, predict(fit, vstHoldout))
    results$MAE[i] <- mae(metaHold$`Neutrophils (%) (%)`, predict(fit, vstHoldout))
  }
  return(results)
}

mi_eval <- model_eval(countSub_mi)
bio_eval <- model_eval(countSub_bio)
combined_eval <- model_eval(countSub_comb)
write.csv(mi_eval, "neutrophil_percent_prediction/mi_eval_split.csv")
write.csv(bio_eval, 'neutrophil_percent_prediction/bio_eval_split.csv')
write.csv(combined_eval, 'neutrophil_percent_prediction/combined_eval_split.csv')


colMeans(combined_eval) 
#  R2      RMSE       MAE 
#0.8152391 3.9585297 3.0373661  
colMeans(bio_eval)
#  R2      RMSE       MAE 
#0.7973307 3.8049697 2.8381028 
colMeans(mi_eval)
#  R2      RMSE       MAE 
#0.4828976 6.0779484 4.8006637
colMeans(xgb_eval)
#  R2      RMSE       MAE 
#0.7536884 4.5898275 3.4590767 


###### Predicted Neutrophil Percentage for all samples

rownames(countTable) <- sub('\\.[0-9]*$', '', rownames(countTable))
countTable <- countTable[, sample_meta$sample_id]
known_sub <- countTable[rownames(countTable) %in% comb_coef$X,]
known_sub <- known_sub[, metaData$sample_id]
rownames(metaData) <- metaData$sample_id

dds <- DESeqDataSetFromMatrix(countData = known_sub,
                              colData = metaData,
                              design = ~1)
dds <- estimateSizeFactors(dds)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vstCounts <- as.data.frame(t(assay(vsd)))
all(rownames(vstCounts) == metaData$sample_id)

vstCounts$neutPer <- metaData$`Neutrophils (%) (%)`
fit <- lm(neutPer~. , data=vstCounts) 
summary(fit) #0.8217

unknown_sub <- countTable[rownames(countTable) %in% comb_coef$X, !(colnames(countTable) %in% metaData$sample_id)]
rownames(sample_meta) <- sample_meta$sample_id
unknown_meta <- sample_meta[colnames(unknown_sub),]
all(unknown_meta$sample_id == colnames(unknown_sub))

dds <- DESeqDataSetFromMatrix(countData = unknown_sub,
                              colData = unknown_meta,
                              design = ~1)
dds <- estimateSizeFactors(dds)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vstCounts_unknown <- as.data.frame(t(assay(vsd)))
all(unknown_meta$sample_id == colnames(unknown_sub))

predNeut <- as.data.frame(predict(fit, vstCounts_unknown))
colnames(predNeut) <- c('neutPer')

knownNeut <- as.data.frame(metaData$`Neutrophils (%) (%)`, row.names = metaData$sample_id)
colnames(knownNeut) <- c('neutPer')

allNeut <- rbind(predNeut, knownNeut)

write.csv(allNeut, 'neutrophil_percent_prediction/knownPredNeutPer_split.csv')


### Comparing significant genes

all_coef
comb_coef$in_bio <- comb_coef$X %in% all_coef$X
comb_coef$in_mi <- comb_coef$X %in% mi_coef$X
mi_coef$in_bio <- mi_coef$X %in% all_coef$X

all_coef$cell <- ifelse(all_coef$X %in% neut$Ensembl, 'Neutrophil', 
                        ifelse(all_coef$X %in% baso$Ensembl, 'Basophil',
                               ifelse(all_coef$X %in% eosino$Ensembl, 'Eosinophil',
                                      ifelse(all_coef$X %in% mono$Ensembl, 'Monocyte',
                                             ifelse(all_coef$X %in% t_cells, 'Lymphocyte',
                                                    ifelse(all_coef$X %in% b_cells, 'Lymphocyte',
                                                           ifelse(all_coef$X %in% den$Ensembl, 'Dendritic Cell', 'None')))))))
comb_coef$cell <- ifelse(comb_coef$X %in% neut$Ensembl, 'Neutrophil', 
                        ifelse(comb_coef$X %in% baso$Ensembl, 'Basophil',
                               ifelse(comb_coef$X %in% eosino$Ensembl, 'Eosinophil',
                                      ifelse(comb_coef$X %in% mono$Ensembl, 'Monocyte',
                                             ifelse(comb_coef$X %in% t_cells, 'Lymphocyte',
                                                    ifelse(comb_coef$X %in% b_cells, 'Lymphocyte',
                                                           ifelse(comb_coef$X %in% den$Ensembl, 'Dendritic Cell', 'None')))))))
mi_coef$cell <- ifelse(mi_coef$X %in% neut$Ensembl, 'Neutrophil', 
                         ifelse(mi_coef$X %in% baso$Ensembl, 'Basophil',
                                ifelse(mi_coef$X %in% eosino$Ensembl, 'Eosinophil',
                                       ifelse(mi_coef$X %in% mono$Ensembl, 'Monocyte',
                                              ifelse(mi_coef$X %in% t_cells, 'Lymphocyte',
                                                     ifelse(mi_coef$X %in% b_cells, 'Lymphocyte',
                                                            ifelse(mi_coef$X %in% den$Ensembl, 'Dendritic Cell', 'None')))))))

