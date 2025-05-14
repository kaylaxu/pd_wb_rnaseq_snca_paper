library(data.table)
library(EnhancedVolcano)
library(readr)
library(ggpubr)
library(patchwork)

setwd('/Users/kxu/Library/CloudStorage/OneDrive-CityofHopeNationalMedicalCenter/CraigLab/Xu/PD Whole Blood DE/PD_FINAL.111524')

mito_genes <- c('MT-LIPCAR',	'MT-ATP6',	'MT-ATP8', 'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-CYB',	'MT-ND1', 'MT-ND2', 'MT-ND3',	'MT-ND4',	'MT-ND4L',	'MT-ND5',	'MT-ND6',	'MT-7SDNA',
                'MT-ATT','MT-CSB1',	'MT-CSB2',	'MT-CSB3',	'MT-HPR',	'MT-HSP1',	'MT-HSP2',	'MT-LSP',	'MT-OHR',	'MT-OLR',	'MT-RNR3',	'MT-TAS', 'MT-TER',	'MT-TFH',	'MT-TFL',
                'MT-TFX',	'MT-TFY',	'MT-RNR1',	'MT-RNR2',	'MT-TA',	'MT-TC',	'MT-TD',	'MT-TE',	'MT-TF',	'MT-TG',	'MT-TH',	'MT-TI',	'MT-TK',	'MT-TL1',	'MT-TL2',	'MT-TM',
                'MT-TN',	'MT-TP',	'MT-TQ',	'MT-TR',	'MT-TS1',	'MT-TS2', 'MT-TT', 'MT-TV', 'MT-TW', 'MT-TY')
pd_genes <- read_tsv('data/PD_genes.tsv')

neut <- fread("data/neutrophil_lineageEnriched.tsv", select = c("Ensembl", "Gene"))
baso <- fread("data/basophil_lineageEnriched.tsv", select = c("Ensembl", "Gene"))
eosino <- fread("data/eosinophil_lineageEnriched.tsv", select = c("Ensembl", "Gene"))
mono <- fread("data/monocytes_lineageEnriched.tsv", select = c("Ensembl", "Gene"))
t_cells <- fread("data/t-cells_lineageEnriched.tsv", select = c("Ensembl", "Gene"))
b_cells <- fread("data/b-cells_lineageEnriched.tsv", select = c("Ensembl", "Gene"))
den <- fread("data/dendritic_lineageEnriched.tsv", select = c("Ensembl", "Gene"))


logfc <- 0.1
pval <- 0.05
annotate_df <- function(df) {
  df$significance <- ifelse(df$logFC > logfc & df$adj.P.Val < pval, 'upregulated', ifelse(df$logFC < -logfc & df$adj.P.Val < pval, 'downregulated', 'insignificant'))
  df$gene_id <- sub('\\.[0-9]*$', '', df$gene_id)
  df$category <- ifelse(df$gene_name %in% pd_genes$gene_name, 'PD causal variant', 
                        ifelse(df$gene_name %in% neut$Gene, 'Neutrophil gene', 
                               ifelse(df$gene_name %in% baso$Gene, 'Basophil gene', 
                                      ifelse(df$gene_name %in% eosino$Gene, 'Eosinophil gene', 
                                             ifelse(df$gene_name %in% t_cells$Gene, 'Lymphocyte gene',
                                                    ifelse(df$gene_name %in% b_cells$Gene, 'Lymphocyte gene',
                                                           ifelse(df$gene_name %in% den$Gene, 'Dendritic cell gene',
                                                                  ifelse(df$gene_name %in% mito_genes, 'Mitochondria gene', 'None'))))))))
  
  df <- df[order(factor(df$category, levels = c('None', 'Basophil gene', 'Eosinophil gene', 'Lymphocyte gene', 'Dendritic cell gene', 'Neutrophil gene', 'Mitochondria gene',  'PD casual variant'))),]
  return(df)
}

createVolcanoCustom <- function(df, plotTitle, keyvals.colour, geneLabels, xlims, ylims, legend='none') {
  EnhancedVolcano(df,
                  lab = df$gene_name,
                  selectLab = geneLabels,
                  boxedLabels = TRUE,
                  drawConnectors = TRUE,
                  x = 'logFC',
                  y = 'P.Value',
                  title = plotTitle,
                  ylab = bquote(~-Log[10]~italic(P-value)),
                  colCustom = keyvals.colour,
                  colAlpha = 0.6,
                  pCutoffCol = 'adj.P.Val',
                  pCutoff = 0.05,
                  FCcutoff = 0.1,
                  pointSize = 2.5,
                  labSize = 3,
                  legendPosition = legend,
                  xlim = xlims,
                  ylim = ylims) +
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_blank())
}




########### IPD v Control #########
df <- read.table(file="DE_allVisits/limmaResults/ipd_v_control_All_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)

df <- annotate_df(df)

keyvals <- rep('gray60', nrow(df))
names(keyvals) <- rep('other', nrow(df))
keyvals[which(df$category == "PD causal variant")] <- '#440154FF'
names(keyvals)[which(df$category == "PD causal variant")] <- 'PD causal variant'
keyvals[which(df$category == "Mitochondria gene")] <- 'darkorange'
names(keyvals)[which(df$category == "Mitochondria gene")] <- 'Mitochondria Genes'
geneLabels <- df[which((df$category == "PD causal variant")& df$significance != 'insignificant'),]$gene_name
ipd <- createVolcanoCustom(df, 'IPD vs Control - All', keyvals, geneLabels, c(-0.56,0.35), c(0,25))
ipd

########### GBA+ v Control #########
df <- read.table(file="DE_allVisits/limmaResults/GBA+_v_control_All_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)

df <- annotate_df(df)

keyvals <- rep('gray60', nrow(df))
names(keyvals) <- rep('other', nrow(df))
keyvals[which(df$category == "PD causal variant")] <- '#440154FF'
names(keyvals)[which(df$category == "PD causal variant")] <- 'PD causal variant'
keyvals[which(df$category == "Mitochondria gene")] <- 'darkorange'
names(keyvals)[which(df$category == "Mitochondria gene")] <- 'Mitochondria Genes'
geneLabels <- df[which((df$category == "PD causal variant")& df$significance != 'insignificant'),]$gene_name
gba <- createVolcanoCustom(df, 'GBA1+ vs Control - All'), keyvals, geneLabels,c(-0.81, 0.61), c(0,18))
gba

########### LRRK2+ v Control #########
df <- read.table(file="DE_allVisits/limmaResults/LRRK2+_v_control_All_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)

df <- annotate_df(df)

keyvals <- rep('gray60', nrow(df))
names(keyvals) <- rep('other', nrow(df))
keyvals[which(df$category == "PD causal variant")] <- '#440154FF'
names(keyvals)[which(df$category == "PD causal variant")] <- 'PD causal variant'
keyvals[which(df$category == "Mitochondria gene")] <- 'darkorange'
names(keyvals)[which(df$category == "Mitochondria gene")] <- 'Mitochondria Genes'
geneLabels <- df[which((df$category == "PD causal variant" )& df$significance != 'insignificant'),]$gene_name
lrrk2 <- createVolcanoCustom(df, 'LRRK2+ vs Control - All', keyvals, geneLabels, c(-1.13,1.1), c(0,75))
lrrk2



########### SNCA+ v Control #########
df <- read.table(file="DE_allVisits/limmaResults/SNCA+_v_control_PPMI_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)

df <- annotate_df(df)

keyvals <- rep('gray60', nrow(df))
names(keyvals) <- rep('other', nrow(df))
keyvals[which(df$category == "PD causal variant")] <- '#440154FF'
names(keyvals)[which(df$category == "PD causal variant")] <- 'PD causal variant'
keyvals[which(df$category == "Mitochondria gene")] <- 'darkorange'
names(keyvals)[which(df$category == "Mitochondria gene")] <- 'Mitochondria Genes'
geneLabels <- df[which((df$category == "PD causal variant" )& df$significance != 'insignificant'),]$gene_name
snca <- createVolcanoCustom(df, 'SNCA+ vs Control - All', keyvals, geneLabels, c(-3.1,3), c(0,25))
snca

########### LRRK2+ v LRRK2- #########
df <- read.table(file="DE_allVisits/limmaResults/LRRK2+_v_LRRK2-_All_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)

df <- annotate_df(df)

keyvals <- rep('gray60', nrow(df))
names(keyvals) <- rep('other', nrow(df))
keyvals[which(df$category == "PD causal variant")] <- '#440154FF'
names(keyvals)[which(df$category == "PD causal variant")] <- 'PD causal variant'
keyvals[which(df$category == "Mitochondria gene")] <- 'darkorange'
names(keyvals)[which(df$category == "Mitochondria gene")] <- 'Mitochondria Genes'
geneLabels <- df[which((df$category == "PD causal variant" )& df$significance != 'insignificant'),]$gene_name
lrrk2_2 <- createVolcanoCustom(df, 'LRRK2+ vs LRRK2- - All', keyvals, geneLabels, c(-1.36, 0.76), c(0,60))
lrrk2_2

########### GBA1+ v GBA1- #########
df <- read.table(file="DE_allVisits/limmaResults/GBA+_v_GBA-_All_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)

df <- annotate_df(df)

keyvals <- rep('gray60', nrow(df))
names(keyvals) <- rep('other', nrow(df))
keyvals[which(df$category == "PD causal variant")] <- '#440154FF'
names(keyvals)[which(df$category == "PD causal variant")] <- 'PD causal variant'
keyvals[which(df$category == "Mitochondria gene")] <- 'darkorange'
names(keyvals)[which(df$category == "Mitochondria gene")] <- 'Mitochondria Genes'
geneLabels <- df[which((df$category == "PD causal variant" )& df$significance != 'insignificant'),]$gene_name
gba_2 <- createVolcanoCustom(df, 'GBA+ vs GBA- - All', keyvals, geneLabels, c(-0.6, 0.6), c(0,12))
gba_2



### SNCA by cohort ###
snca_data <- read.csv('snca/snca_plot_data.csv')
limma_pval <- data.frame(group1=rep('HC', 4), group2=c('IPD', 'SNCA+ PD', 'GBA+ PD', 'LRRK2+ PD'), adj.p=c(1.1e-05, 1.7e-08, 6.4e-05, 1.9e-08), sig=c('***', '***', '***', '***'))
snca_data$cohort <- ifelse(snca_data$cohort %in% c('SNCA+', 'GBA+', 'LRRK2+'),
                           ifelse(snca_data$case_control_other %in% c('Case'), snca_data$cohort, NA), snca_data$cohort)
snca_data$cohort <- ifelse(snca_data$cohort %in% c('SNCA+', 'GBA+', 'LRRK2+'), paste0(snca_data$cohort, ' PD'), snca_data$cohort)

p <- ggplot(snca_data[snca_data$cohort %in% c('HC', 'IPD', 'SNCA+ PD', 'GBA+ PD', 'LRRK2+ PD'),], aes(x=factor(cohort, levels= c('HC', 'IPD', 'SNCA+ PD', 'GBA+ PD', 'LRRK2+ PD')), y=SNCA)) +
  geom_violin(alpha=0.5, aes(fill = cohort)) +
  geom_boxplot(width=0.25, alpha=0.5,  aes(fill = cohort)) +
  theme_classic()+
  theme(text = element_text(size=15, colour = 'black')) +
  stat_pvalue_manual(limma_pval, label = 'sig', y.position = 13.5, step.increase = 0.05, label.size = 6) +
  geom_hline(yintercept = 9.948, linetype='dotted', col = 'red', linewidth=0.8) +
  scale_fill_brewer(palette="Set2") +
  xlab('') +
  ylab(substitute(paste(italic('SNCA')))) +
  ggtitle('Expression by Genetic Status') +
  guides(fill='none')




ggsave('figure4.jpg', plot = ipd  + gba + lrrk2 + p + snca + lrrk2_2 + gba_2 + p + plot_layout(nrow = 2), dpi=300, height = 16, width=25)






