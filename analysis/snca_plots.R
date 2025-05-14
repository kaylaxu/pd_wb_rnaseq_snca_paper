##### SNCA Over Demongraphics  #####

library(readxl)
library(tidyverse)
library(limma)
library(edgeR)
library(umap)
library(ggbeeswarm)
library(ggpubr)
library(DESeq2)

countTable <- read.table(file="data/matrix.featureCountsV2.tsv", header=TRUE, sep="\t", check.names=FALSE, row.names="Geneid")
sample_meta <- read.csv(file="data/AMP_RNAseq_metaData.csv", header=TRUE)
sample_meta <- sample_meta[sample_meta$PCT_CHIMERAS < 0.03,]
sample_meta <- sample_meta[sample_meta$study != 'BioFIND',]
sample_meta <- sample_meta[sample_meta$case_control_other_at_baseline %in% c('Case', 'Control'),]
sample_meta$ageSquared <- sample_meta$age_at_baseline**2

predNeutPer <- read.csv('neutrophil_percent_prediction/knownPredNeutPer_split.csv')
rownames(predNeutPer) <- predNeutPer$X
predNeutPer <- predNeutPer[sample_meta$sample_id,]
sample_meta$neutPer <- predNeutPer$neutPer
rownames(sample_meta) <- sample_meta$sample_id

countSub <- countTable[, sample_meta$sample_id]
all(colnames(countSub) == sample_meta$sample_id) # TRUE

tmpG <- data.table::fread("../AMP_RNAseq/GRCh38_GENCODE29_geneInfo.txt", skip = 1)
genes.anno <- tmpG[,c("gene_id", "gene_name", "gene_type", "seqnames", "start", "end", "width", "strand")]
genes.anno <- genes.anno[with(genes.anno, order(gene_id, -width)), ]
head(genes.anno, 12)
genes.anno <- distinct(genes.anno, gene_id, .keep_all = TRUE)
countSub <- countSub[genes.anno$gene_id, ] 

dge <- DGEList(countSub) 
dge  <- calcNormFactors(dge)

#order gene list

genes.anno.Match <- genes.anno[genes.anno$gene_id %in% sort(row.names(dge$counts)),]
dim(genes.anno.Match)

i <- match(row.names(dge$counts), genes.anno.Match$gene_id, nomatch = 0)
genes.anno.Match <- genes.anno.Match[i,]

identical(genes.anno.Match$gene_id, row.names(dge$counts)) #double checking


dge$genes <- genes.anno.Match

CPM <- cpm(dge)

logCPM <- as.data.frame(cpm(dge, log = TRUE, prior.count = 3))
all(colnames(logCPM) == sample_meta$sample_id)

logCPM_neutPer <- as.data.frame(removeBatchEffect(logCPM, covariates = sample_meta$neutPer))


snca_data <- as.data.frame(t(logCPM[c('ENSG00000145335'),]))
colnames(snca_data)[1] <- 'SNCA'
all(rownames(snca_data) == sample_meta$sample_id)

snca_data$GBA_mut <- sample_meta$has_known_GBA_mutation_in_WGS
snca_data$SNCA_mut <- sample_meta$has_known_SNCA_mutation_in_WGS
snca_data$LRRK2_mut <- sample_meta$has_known_LRRK2_mutation_in_WGS
snca_data$PD_mut <- sample_meta$has_known_PD_mutation_in_WGS
snca_data$case_control_other <- factor(sample_meta$case_control_other_at_baseline, levels = c('Control', 'Case', 'Other'))

snca_data$genetic_status <- ifelse(is.na(snca_data$PD_mut), 'Unknown', 
                                   ifelse(snca_data$SNCA_mut == 'Yes', 'SNCA+', 
                                          ifelse(snca_data$LRRK2_mut == 'Yes' & snca_data$GBA_mut == 'No', 'LRRK2+',
                                                 ifelse(snca_data$LRRK2_mut == 'No' & snca_data$GBA_mut == 'Yes', 'GBA+',
                                                        ifelse(snca_data$LRRK2_mut == 'Yes' & snca_data$GBA_mut == 'Yes', 'GBA+/LRRK2+', 'SNCA-/GBA-/LRRK2-')))))
snca_data$age_at_baseline <- sample_meta$age_at_baseline
snca_data$ageBin <- factor(ifelse(snca_data$age_at_baseline < 30, NA,
                                  ifelse(snca_data$age_at_baseline < 40, '30s',
                                         ifelse(snca_data$age_at_baseline < 50, '40s', 
                                                ifelse(snca_data$age_at_baseline < 60, '50s',
                                                       ifelse(snca_data$age_at_baseline < 70, '60s', 
                                                              ifelse(snca_data$age_at_baseline < 80, '70s', '>80')))))),
                           levels = c('<30', '30s', '40s', '50s', '60s', '70s', '>80'))
  

snca_data$visit <- factor(sample_meta$visit_month,levels = c('0', '6', '12', '18', '24', '36'))


### Corrected for predicted neutrophil percentage
snca_data_neutPer <- snca_data
snca_neutPer <- as.data.frame(t(logCPM_neutPer[c('ENSG00000145335'),]))
snca_data_neutPer$SNCA <- snca_neutPer$ENSG00000145335.15


### SNCA by genetic cohort
# x = hc, ipd, snca, gba, lrrk
# y = snca
ipd_df <- read.table(file="DE_allVisits/limmaResults/ipd_v_control_All_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)
ipd_df <- ipd_df[ipd_df$gene_name %in% c('SNCA'),] #1.103721e-05

snca_df <- read.table(file="DE_allVisits/limmaResults/SNCA+_v_control_PPMI_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)
snca_df <- snca_df[snca_df$gene_name %in% c('SNCA'),] #1.704208e-08

gba_df <- read.table(file="DE_allVisits/limmaResults/GBA+_v_control_All_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)
gba_df <- gba_df[gba_df$gene_name %in% c('SNCA'),] # 6.419258e-05

lrrk2_df <- read.table(file="DE_allVisits/limmaResults/LRRK2+_v_control_All_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)
lrrk2_df <- lrrk2_df[lrrk2_df$gene_name %in% c('SNCA'),] #1.902129e-08

limma_pval <- data.frame(group1=rep('HC', 4), group2=c('IPD', 'SNCA+', 'GBA+', 'LRRK2+'), adj.p=c(1.1e-05, 1.7e-08, 6.4e-05, 1.9e-08))


snca_data$cohort <- factor(ifelse(snca_data$genetic_status == 'SNCA-/GBA-/LRRK2-', 
                                          ifelse(snca_data$case_control_other == 'Case', 'IPD', 'HC'), snca_data$genetic_status), levels = c('HC', 'IPD', 'SNCA+', 'GBA+', 'LRRK2+', 'GBA+/LRRK2+', 'Unknown'))
snca_data_neutPer$cohort <- factor(ifelse(snca_data_neutPer$genetic_status == 'SNCA-/GBA-/LRRK2-', 
                                   ifelse(snca_data_neutPer$case_control_other == 'Case', 'IPD', 'HC'), snca_data_neutPer$genetic_status), levels = c('HC', 'IPD', 'SNCA+', 'GBA+', 'LRRK2+', 'GBA+/LRRK2+', 'Unknown'))

p <- ggplot(snca_data_neutPer[snca_data_neutPer$cohort %in% c('HC', 'IPD', 'SNCA+', 'GBA+', 'LRRK2+'),], aes(x=cohort, y=SNCA)) +
  geom_violin(alpha=0.5, aes(fill = cohort)) +
  geom_boxplot(width=0.25, alpha=0.5,  aes(fill = cohort)) +
  theme_linedraw()+
  theme(text = element_text(size=20)) +
  stat_pvalue_manual(limma_pval, label = 'adj.p', y.position = 13.5, step.increase = 0.05, label.size = 6) +
  xlab('') +
  ylab('log(CPM)') +
  ggtitle('SNCA by Genetic Cohort')

ggsave('snca/snca_by_cohort.jpeg', p, dpi=300, width = 10, height = 8)


### SNCA by predicted neutrophil percentage
snca_data$predNeut <- sample_meta$neutPer

p1 <- ggplot(snca_data[snca_data$cohort %in% c('HC', 'IPD', 'SNCA+', 'GBA+', 'LRRK2+'),], aes(x=predNeut, y=SNCA)) +
  geom_point(size = 0.5, aes(color = case_control_other), alpha=0.5) +
  geom_smooth(method = 'lm', aes(color = case_control_other)) +
  scale_color_manual(values = c("red", "blue")) +
  theme_bw() +
  theme(text = element_text(size=20)) +
  xlab('Predicted Neutrophil Percentage')+
  ylab(substitute(paste(italic('SNCA'))))+
  guides(color = guide_legend(title = "Status"))  +
  stat_cor(aes(color = case_control_other), label.x = 3) 
  
p2 <- ggplot(snca_data[snca_data$cohort %in% c('HC', 'IPD', 'SNCA+', 'GBA+', 'LRRK2+'),], aes(x=predNeut, y=SNCA)) +
  geom_point(size = 0.5, aes(color = GBA_mut), alpha=0.5) +
  geom_smooth(method = 'lm', aes(color = GBA_mut)) +
  scale_color_manual(values = c("darkgreen", "darkorange")) +
  theme_bw() +
  theme(text = element_text(size=20)) +
  xlab('Predicted Neutrophil Percentage')+
  ylab(substitute(paste(italic('SNCA'))))+
  guides(color = guide_legend(title = substitute(paste(italic('GBA1+')))))  +
  stat_cor(aes(color = GBA_mut), label.x = 3) 

p3 <- ggplot(snca_data[snca_data$cohort %in% c('HC', 'IPD', 'SNCA+', 'GBA+', 'LRRK2+'),], aes(x=predNeut, y=SNCA, color = SNCA_mut)) +
  geom_point(size = 0.5, aes(color = SNCA_mut), alpha=0.5) +
  geom_smooth(method = 'lm', aes(color = SNCA_mut)) +
  scale_color_manual(values = c("darkgreen", "darkorange")) +
  theme_bw() +
  theme(text = element_text(size=20)) +
  xlab('Predicted Neutrophil Percentage')+
  ylab(substitute(paste(italic('SNCA'))))+
  guides(color = guide_legend(title = substitute(paste(italic('SNCA+')))))  +
  stat_cor(aes(color = SNCA_mut), label.x = 3) 

p4 <- ggplot(snca_data[snca_data$cohort %in% c('HC', 'IPD', 'SNCA+', 'GBA+', 'LRRK2+'),], aes(x=predNeut, y=SNCA, color = LRRK2_mut)) +
  geom_point(size = 0.5, aes(color = LRRK2_mut), alpha=0.5) +
  geom_smooth(method = 'lm', aes(color = LRRK2_mut)) +
  scale_color_manual(values = c("darkgreen", "darkorange")) +
  theme_bw() +
  theme(text = element_text(size=20)) +
  xlab('Predicted Neutrophil Percentage') +
  ylab(substitute(paste(italic('SNCA'))))+
  guides(color = guide_legend(title = substitute(paste(italic('LRR2+')))))  +
  stat_cor(aes(color = LRRK2_mut), label.x = 3) 

ggsave('snca/snca_neutPer.jpg', plot=p1+p2+p3+p4, dpi=300)


### Figure 6 
snca_data_neutPer$genetic_status <- factor(snca_data_neutPer$genetic_status, levels = c('GBA+', 'LRRK2+', 'SNCA+', 'SNCA-/GBA-/LRRK2-'))
#write.csv(snca_data_neutPer, 'snca/snca_plot_data.csv')

p1 <- ggplot(snca_data_neutPer[snca_data_neutPer$cohort %in% c('HC', 'IPD', 'SNCA+', 'GBA+', 'LRRK2+') & !is.na(snca_data_neutPer$ageBin),], aes(x=ageBin, y=SNCA, fill=case_control_other)) +
  geom_boxplot(width=0.4) +
  facet_grid(genetic_status~case_control_other) +
  geom_hline(yintercept=9.948, linetype='dotted', col = 'red', linewidth=0.8)+
  theme_classic() +
  theme(text = element_text(size=20)) +
  xlab('Age at Baseline') 

snca_data_neutPer$cohort2 <- ifelse(snca_data_neutPer$cohort %in% c('HC', 'IPD'), as.character(snca_data_neutPer$cohort),
                                    ifelse(snca_data_neutPer$case_control_other == 'Control', NA, 
                                           ifelse(snca_data_neutPer$cohort %in% c('SNCA+'), 'SNCA+ PD',
                                                  ifelse(snca_data_neutPer$cohort %in% c('GBA+'), 'GBA+ PD',
                                                         ifelse(snca_data_neutPer$cohort %in% c('LRRK2+'), 'LRRK2+ PD', NA)))))

p2 <- ggplot(snca_data_neutPer[snca_data_neutPer$visit %in% c(0) & snca_data_neutPer$cohort2 %in% c('HC', 'IPD', 'SNCA+ PD', 'GBA+ PD', 'LRRK2+ PD'),], aes(x=factor(cohort2, levels = c('HC', 'IPD', 'SNCA+ PD', 'GBA+ PD', 'LRRK2+ PD')), y=SNCA, fill=cohort2)) +
  geom_violin(alpha=0.5) +
  geom_boxplot(width=0.15) +
  theme_classic() +
  scale_fill_brewer(palette="Set2") +
  theme(text = element_text(size=20)) +
  stat_compare_means(comparisons = list(c('HC', 'IPD'),
                                 c('HC', 'SNCA+ PD'),
                                 c('HC', 'GBA+ PD'),
                                 c('HC', 'LRRK2+ PD')),
              test='wilcox.test', label = 'p.signif') +
  xlab('Cohort') +
  ggtitle('X')


p3 <- ggplot(snca_data_neutPer[snca_data_neutPer$cohort %in% c('IPD') & !is.na(snca_data_neutPer$ageBin),], aes(x=ageBin, y=SNCA, fill=ageBin)) +
  geom_violin(alpha=0.5)+
  geom_boxplot(width=0.15) +
  theme_classic() +
  theme(text = element_text(size=20)) +
  scale_fill_brewer(palette="Oranges") +
  stat_compare_means(comparisons = list(c('30s', '40s'),
                                 c('30s', '50s'),
                                 c('30s', '60s'),
                                 c('30s', '70s'),
                                 c('30s', '>80')),
              test='wilcox.test', label = 'p.signif') +
  xlab('Age at Baseline')+
  ggtitle('X')

p <- ggarrange(p1 + guides(fill="none"), ggarrange(p2 + guides(fill="none"), p3 + guides(fill="none"), nrow = 2, labels = c("B", "C")), 
          labels = c("a", "b", "c"),
          ncol = 2)

ggsave('snca/test.jpeg', plot=p, dpi=300, height = 15, width=15)

### All case vs all control snca expression by age

p <- ggplot(snca_data_neutPer[snca_data_neutPer$visit == 0,], aes(x=case_control_other, y=SNCA, fill=case_control_other)) +
  geom_violin(alpha=0.5)+
  geom_boxplot(alpha=0.5, width=0.25) +
  facet_wrap(~ageBin, ncol = 5)+
  stat_compare_means(comparisons = list(c('Control', 'Case')))+
  theme_bw() +
  theme(text = element_text(size=20)) +
  ylab('logCPM') +
  xlab('')+
  ggtitle('SNCA Expression by Age at Baseline - Neutrophil % Correction')


ggsave('snca/case_v_control_SNCA_exp_neutPer.jpeg', p, dpi=300)



### All case vs all control snca expression by visit

p <- ggplot(snca_data_neutPer, aes(x=case_control_other, y=SNCA, fill = case_control_other)) +
  geom_boxplot(alpha=0.5, width=0.5) +
  geom_violin(alpha=0.5)+
  facet_grid(~visit) +
  geom_signif(comparisons = list(c('Case', 'Control')),
              test = 'wilcox.test', step_increase = 0.04, textsize = 5) +
  theme_bw() +
  xlab('')+
  ylab('logCPM') +
  theme(text = element_text(size=20))

ggsave('snca/snca_by_visit_neutPer.jpeg', p, dpi=300)



















