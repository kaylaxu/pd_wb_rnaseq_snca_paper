###### UMAP Plotting and Analysis ######

library(readxl)
library(tidyverse)
library(limma)
library(DESeq2)
library(umap)

countTable <- read.table(file="data/matrix.featureCountsV2.tsv", header=TRUE, sep="\t", check.names=FALSE, row.names="Geneid")
sample_meta <- read.csv(file="data/AMP_RNAseq_metaData.csv", header=TRUE)
sample_meta <- sample_meta[sample_meta$PCT_CHIMERAS < 0.03,]
sample_meta <- sample_meta[sample_meta$study != 'BioFIND',]
sample_meta <- sample_meta[sample_meta$case_control_other_at_baseline %in% c('Case', 'Control'),]
sample_meta$ageSquared <- sample_meta$age_at_baseline**2

sample_meta$heatmap_clusters <- ifelse(sample_meta$case_control_other_at_baseline == 'Control' & sample_meta$has_known_PD_mutation_in_WGS == 'No', "Control",
                                       ifelse(sample_meta$case_control_other_at_baseline == 'Case',
                                              ifelse(sample_meta$has_known_GBA_mutation_in_WGS == 'Yes' & sample_meta$has_known_LRRK2_mutation_in_WGS == 'No', "GBA+",
                                                     ifelse(sample_meta$has_known_LRRK2_mutation_in_WGS == 'Yes' & sample_meta$has_known_GBA_mutation_in_WGS == 'No', "LRRK2+",
                                                            ifelse(sample_meta$has_known_SNCA_mutation_in_WGS == 'Yes', "SNCA+", 
                                                                   ifelse(sample_meta$has_known_GBA_mutation_in_WGS == 'Yes' & sample_meta$has_known_LRRK2_mutation_in_WGS == 'Yes', NA, 'Idiopathic PD')))), NA))

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
dds <- estimateSizeFactors(dds)
normCounts <- as.data.frame(counts(dds, normalize = TRUE))
vsd <- vst(dds, blind = TRUE)
vstCounts_corrected <- removeBatchEffect(assay(vsd), 
                                         batch = sample_meta$sex, 
                                         covariates = cbind(sample_meta$neutPer, sample_meta$PCT_MRNA_BASES, sample_meta$ageSquared))
vstCounts_corrected <- as.data.frame(vstCounts_corrected)
vstCounts <- as.data.frame(assay(vsd))


mito_genes <- c('MT-LIPCAR',	'MT-ATP6',	'MT-ATP8', 'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-CYB',	'MT-ND1', 'MT-ND2', 'MT-ND3',	'MT-ND4',	'MT-ND4L',	'MT-ND5',	'MT-ND6',	'MT-7SDNA',
                'MT-ATT','MT-CSB1',	'MT-CSB2',	'MT-CSB3',	'MT-HPR',	'MT-HSP1',	'MT-HSP2',	'MT-LSP',	'MT-OHR',	'MT-OLR',	'MT-RNR3',	'MT-TAS', 'MT-TER',	'MT-TFH',	'MT-TFL',
                'MT-TFX',	'MT-TFY',	'MT-RNR1',	'MT-RNR2',	'MT-TA',	'MT-TC',	'MT-TD',	'MT-TE',	'MT-TF',	'MT-TG',	'MT-TH',	'MT-TI',	'MT-TK',	'MT-TL1',	'MT-TL2',	'MT-TM',
                'MT-TN',	'MT-TP',	'MT-TQ',	'MT-TR',	'MT-TS1',	'MT-TS2', 'MT-TT', 'MT-TV', 'MT-TW', 'MT-TY')

tmpG <- data.table::fread("../AMP_RNAseq/GRCh38_GENCODE29_geneInfo.txt", skip = 1)
genes.anno <- tmpG[,c("gene_id", "gene_name", "gene_type", "seqnames", "start", "end", "width", "strand")]
genes.anno <- genes.anno[with(genes.anno, order(gene_id, -width)), ]
head(genes.anno, 12)
genes.anno <- distinct(genes.anno, gene_id, .keep_all = TRUE)


# Genes from significant IPA pathways
mito_dys_genes <- read_excel('pathway_analysis/mito_dys_genes.xls')
mito_dys_genes <- mito_dys_genes[!is.na(mito_dys_genes$`Entrez Gene ID for Human`),]
pd_genes <- read_excel('pathway_analysis/pd_genes.xls')
pd_genes <- pd_genes[!is.na(pd_genes$`Entrez Gene ID for Human`),]
leuk_genes <- read_excel('pathway_analysis/leukocyte_extravasation_genes.xls')
leuk_genes <- leuk_genes[!is.na(leuk_genes$`Entrez Gene ID for Human`),]
bbsome_genes <- read_excel('pathway_analysis/bbsome_genes.xls')
bbsome_genes <- bbsome_genes[!is.na(bbsome_genes$`Entrez Gene ID for Human`),]


mito.anno <- genes.anno[genes.anno$gene_name %in% mito_genes,]
mito.dys.anno <- genes.anno[genes.anno$gene_name %in% mito_dys_genes$Symbol,]
pd.anno <- genes.anno[genes.anno$gene_name %in% pd_genes$Symbol,]
leuk.anno <- genes.anno[genes.anno$gene_name %in% leuk_genes$Symbol,]
bbsome.anno <- genes.anno[genes.anno$gene_name %in% bbsome_genes$Symbol,]


####### Hierarchical Clustering #################
library(pheatmap)
x <- vstCounts_corrected[pd.anno$gene_id,]
dfh <- sample_meta[!is.na(sample_meta$heatmap_clusters),c('sample_id', 'heatmap_clusters')]
rownames(dfh) <- dfh$sample_id
dfh$sample_id <- NULL
x <- x[, rownames(dfh)]

jpeg(filename = 'hierarchical_clustering/pd_clustering.jpg', quality = 300, height = 1000, width=3000)
pheatmap(x, scale = 'row', annotation_col= dfh, color= colorRampPalette(c("navy", "white", "red"))(50), clustering_method = "ward.D")
dev.off()


#### UMAP from PD pathway genes
pd_counts <- vstCounts_corrected[pd.anno$gene_id,]
pd_umap <- umap(t(pd_counts))
pd_layout <- data.frame(pd_umap$layout)
pd_layout$cluster <- sample_meta[rownames(pd_layout),]$heatmap_clusters
pd_layout <- as.data.frame(pd_layout[!is.na(pd_layout$cluster),])

p <- ggplot(pd_layout, aes(x = X1, y = X2, col = cluster))  +
  geom_point(pch=19, size=1.5, alpha=0.7) +
  theme_void()

ggsave('umap/umap_pd_corrected.jpeg', p, dpi = 600)

density_plot <- ggplot(pd_layout[!is.na(pd_layout$cluster),], aes(x=X1, y=X2)) +
  stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=F) + 
  geom_point(data=pd_layout[, c('X1', 'X2')], shape=16, size=0.1, alpha=0.2, color="white") +
  scale_fill_gradientn(colours=viridisLite::mako(100), name="Density") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(~factor(cluster, levels = c('Control', 'Idiopathic PD', 'SNCA+', 'LRRK2+', 'GBA+')), ncol=6) +
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12, color="black"),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"))

ggsave('umap/density_pd_corrected.jpeg', density_plot, dpi = 600, height = 6, width = 32)


#### UMAP from mitochondria genes
mito_counts <- vstCounts_corrected[mito.anno$gene_id,]
mito_umap <- umap(t(mito_counts))
mito_layout <- data.frame(mito_umap$layout)
mito_layout$cluster <- sample_meta[rownames(mito_layout),]$heatmap_clusters
mito_layout <- as.data.frame(mito_layout[!is.na(mito_layout$cluster),])

p <- ggplot(mito_layout, aes(x = X1, y = X2, col = cluster))  +
  geom_point(pch=19, size=2, alpha=0.7) +
  theme_void()

ggsave('umap/umap_mito_corrected.jpeg', p, dpi = 600)

density_plot <- ggplot(mito_layout[!is.na(mito_layout$cluster),], aes(x=X1, y=X2)) +
  stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=F) + 
  geom_point(data=mito_layout[, c('X1', 'X2')], shape=16, size=0.2, alpha=0.2, color="white") +
  scale_fill_gradientn(colours=viridisLite::mako(100), name="Density") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(~factor(cluster, levels = c('Control', 'Idiopathic PD', 'SNCA+', 'LRRK2+', 'GBA+')), ncol=6) +
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12, color="black"),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"))

ggsave('umap/density_mito_corrected.jpeg', density_plot, dpi = 600, height = 6, width = 32)



#### UMAP from mitochondria dysfunction pathway genes
mito_counts <- vstCounts_corrected[mito.dys.anno$gene_id,]
mito_umap <- umap(t(mito_counts))
mito_layout <- data.frame(mito_umap$layout)
mito_layout$cluster <- sample_meta[rownames(mito_layout),]$heatmap_clusters
mito_layout <- as.data.frame(mito_layout[!is.na(mito_layout$cluster),]) #5470 
mito_layout <- mito_layout[mito_layout$X2 > -10,] #5339

p <- ggplot(mito_layout, aes(x = X1, y = X2, col = cluster))  +
  geom_point(pch=19, size=2, alpha=0.7) +
  theme_void()

ggsave('umap/umap_mitoDys_corrected.jpeg', p, dpi = 600)


density_plot <- ggplot(mito_layout[!is.na(mito_layout$cluster),], aes(x=X1, y=X2)) +
  stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=F) + 
  geom_point(data=mito_layout[, c('X1', 'X2')], shape=16, size=0.2, alpha=0.2, color="white") +
  scale_fill_gradientn(colours=viridisLite::mako(100), name="Density") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(~factor(cluster, levels = c('Control', 'Idiopathic PD', 'SNCA+', 'LRRK2+', 'GBA+')), ncol=6) +
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12, color="black"),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"))

ggsave('umap/density_mitoDys_corrected.jpeg',density_plot, dpi = 600, height = 6, width = 32)


#### UMAP from Leukocyte Extravasation Signaling   BEST 
leuk_counts <- vstCounts_corrected[leuk.anno$gene_id,]
leuk_umap <- umap(t(leuk_counts))
leuk_layout <- data.frame(leuk_umap$layout)
leuk_layout$cluster <- sample_meta[rownames(leuk_layout),]$heatmap_clusters
leuk_layout <- as.data.frame(leuk_layout[!is.na(leuk_layout$cluster),])
leuk_layout <- leuk_layout[leuk_layout$X2 >-8,] #5339

p <- ggplot(leuk_layout, aes(x = X1, y = X2, col = cluster))  +
  geom_point(pch=19, size=2, alpha=0.7) +
  theme_void()
ggsave('umap/umap_leuk_corrected.jpeg', p, dpi = 600)


density_plot <- ggplot(leuk_layout[!is.na(leuk_layout$cluster),], aes(x=X1, y=X2)) +
  stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=F) + 
  geom_point(data=leuk_layout[, c('X1', 'X2')], shape=16, size=0.2, alpha=0.2, color="white") +
  scale_fill_gradientn(colours=viridisLite::mako(100), name="Density") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(~factor(cluster, levels = c('Control', 'Idiopathic PD', 'SNCA+', 'LRRK2+', 'GBA+')), ncol=6) +
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12, color="black"),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"))

ggsave('umap/density_leuk_corrected.jpeg', density_plot, dpi = 600, height = 6, width = 32)


#### UMAP from BBSome Signaling Pathway
bbsome_counts <- vstCounts_corrected[bbsome.anno$gene_id,]
bbsome_umap <- umap(t(bbsome_counts))
bbsome_layout <- data.frame(bbsome_umap$layout)
bbsome_layout$cluster <- sample_meta[rownames(bbsome_layout),]$heatmap_clusters
bbsome_layout <- as.data.frame(bbsome_layout[!is.na(bbsome_layout$cluster),])
bbsome_layout <- bbsome_layout[bbsome_layout$X1 < 5,] #5358

p <- ggplot(bbsome_layout, aes(x = X1, y = X2, col = cluster))  +
  geom_point(pch=19, size=2, alpha=0.7) +
  theme_void()
ggsave('umap/umap_bbsome_corrected.jpeg', p, dpi = 600)


density_plot <- ggplot(bbsome_layout[!is.na(bbsome_layout$cluster),], aes(x=X1, y=X2)) +
  stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=F) + 
  geom_point(data=bbsome_layout[, c('X1', 'X2')], shape=16, size=0.2, alpha=0.2, color="white") +
  scale_fill_gradientn(colours=viridisLite::mako(100), name="Density") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(~factor(cluster, levels = c('Control', 'Idiopathic PD', 'SNCA+', 'LRRK2+', 'GBA+')), ncol=6) +
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12, color="black"),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"))

ggsave('umap/density_bbsome_corrected.jpeg', density_plot, dpi = 600, height = 6, width = 32)























