###### Volcano Plots for All Visits w/Neutrophil Correction at Baseline ########

library(data.table)
library(EnhancedVolcano)
library(readr)

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
  
  df <- df[order(factor(df$category, levels = c('None', 'Basophil gene', 'Eosinophil gene', 'Lymphocyte gene', 'Dendritic cell gene', 'Neutrophil gene', 'Mitochondria gene', 'Mitochondria-related pathway gene', 'PD casual variant'))),]
  return(df)
}

createVolcanoCustom <- function(df, plotTitle, keyvals.colour, geneLabels, xlims, ylims) {
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
                  legendPosition = "bottom",
                  xlim = xlims,
                  ylim = ylims, max.overlaps = 100) +
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_blank())
}



########### Case v Control #########
df <- read.table(file="DE_BL/limmaResults/case_v_control_All_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)

df <- annotate_df(df)
volPlot <- EnhancedVolcano(df,
                           lab = df$gene_name,
                           x = 'logFC',
                           y = 'P.Value',
                           title = 'Case vs Control - All at Baseline',
                           ylab = bquote(~-Log[10]~italic(P)),
                           colAlpha = 0.6,
                           pCutoffCol = 'adj.P.Val',
                           pCutoff = 0.05,
                           FCcutoff = 0.1,
                           pointSize = 2.5,
                           labSize = 0,
                           legendPosition = "bottom",
                           xlim = c(-0.7, 0.35),
                           ylim = c(0,13)) +
  theme(plot.title = element_text(size = 20),
        plot.subtitle = element_blank())

volPlot

ggsave('DE_BL/volcanoPlots/case_v_control_All_EVol.jpeg',
       volPlot,
       dpi = 300, width = 6, height = 8)

# Mito volcano plot
keyvals <- rep('gray60', nrow(df))
names(keyvals) <- rep('other', nrow(df))
keyvals[which(df$category == "PD causal variant")] <- '#440154FF'
names(keyvals)[which(df$category == "PD causal variant")] <- 'PD causal variant'
keyvals[which(df$category == "Mitochondria gene")] <- 'darkorange'
names(keyvals)[which(df$category == "Mitochondria gene")] <- 'Pathway Genes'
geneLabels <- df[which((df$category == "PD causal variant" | df$category == 'Mitochondria gene')& df$significance != 'insignificant'),]$gene_name
volPlot <- createVolcanoCustom(df, 'Case vs Control -  All at Baseline', keyvals, geneLabels, c(-0.7,0.35), c(0,13))
volPlot
ggsave('DE_BL/volcanoPlots/case_v_control_All_EVol_mito.jpeg',
       volPlot,
       dpi = 300, width = 6, height = 8)

# blood cell plot
df_sig <- df[df$adj.P.Val < pval & abs(df$logFC) > logfc ,] 
df_sig$category <- ifelse(df_sig$category %in% c('Lymphocyte gene', 'Neutrophil gene'),  df_sig$category, NA)
plot <- ggplot(df_sig[!is.na(df_sig$category),], aes(x=logFC, fill=category)) +
  geom_histogram(alpha=0.7, bins = 100) +
  xlim(-0.7,0.35) + 
  theme_light() +
  theme(text = element_text(size=15, color='black'),
        axis.text = element_text(size=15, color='black')) 
plot

ggsave('DE_BL/blood_cell_plots/case_v_control_All_blood_genes.jpeg',
       plot,
       dpi = 300, width = 6, height = 3)


########### IPD v Control #########
df <- read.table(file="DE_BL/limmaResults/ipd_v_control_All_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)

df <- annotate_df(df)
volPlot <- EnhancedVolcano(df,
                           lab = df$gene_name,
                           x = 'logFC',
                           y = 'P.Value',
                           title = 'IPD vs Control - All at Baseline',
                           ylab = bquote(~-Log[10]~italic(P)),
                           colAlpha = 0.6,
                           pCutoffCol = 'adj.P.Val',
                           pCutoff = 0.05,
                           FCcutoff = 0.1,
                           pointSize = 2.5,
                           labSize = 0,
                           legendPosition = "bottom",
                           xlim = c(-0.74, 0.26),
                           ylim = c(0,7)) +
  theme(plot.title = element_text(size = 20),
        plot.subtitle = element_blank())

volPlot

ggsave('DE_BL/volcanoPlots/ipd_v_control_All_EVol.jpeg',
       volPlot,
       dpi = 300, width = 6, height = 8)

# Mito volcano plot
keyvals <- rep('gray60', nrow(df))
names(keyvals) <- rep('other', nrow(df))
keyvals[which(df$category == "PD causal variant")] <- '#440154FF'
names(keyvals)[which(df$category == "PD causal variant")] <- 'PD causal variant'
keyvals[which(df$category == "Mitochondria gene")] <- 'darkorange'
names(keyvals)[which(df$category == "Mitochondria gene")] <- 'Pathway Genes'
geneLabels <- df[which((df$category == "PD causal variant" | df$category == 'Mitochondria gene')& df$significance != 'insignificant'),]$gene_name
volPlot <- createVolcanoCustom(df, 'IPD vs Control - All at Baseline', keyvals, geneLabels, c(-0.74,0.26), c(0,7))
volPlot
ggsave('DE_BL/volcanoPlots/ipd_v_control_All_EVol_mito.jpeg',
       volPlot,
       dpi = 300, width = 6, height = 8)

# blood cell plot
df_sig <- df[df$adj.P.Val < pval & abs(df$logFC) > logfc ,] 
df_sig$category <- ifelse(df_sig$category %in% c('Lymphocyte gene', 'Neutrophil gene'),  df_sig$category, NA)
plot <- ggplot(df_sig[!is.na(df_sig$category),], aes(x=logFC, fill=category)) +
  geom_histogram(alpha=0.7, bins = 100) +
  xlim(-0.74,0.26) + 
  theme_light() +
  theme(text = element_text(size=15, color='black'),
        axis.text = element_text(size=15, color='black')) 
plot

ggsave('DE_BL/blood_cell_plots/ipd_v_control_All_blood_genes.jpeg',
       plot,
       dpi = 300, width = 6, height = 3)




########### GBA+ v Control #########
df <- read.table(file="DE_BL/limmaResults/GBA+_v_control_All_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)

df <- annotate_df(df)
volPlot <- EnhancedVolcano(df,
                           lab = df$gene_name,
                           x = 'logFC',
                           y = 'P.Value',
                           title = 'GBA+ vs Control - All at Baseline',
                           ylab = bquote(~-Log[10]~italic(P)),
                           colAlpha = 0.6,
                           pCutoffCol = 'adj.P.Val',
                           pCutoff = 0.05,
                           FCcutoff = 0.1,
                           pointSize = 2.5,
                           labSize = 0,
                           legendPosition = "bottom",
                           xlim = c(-1.1, 0.8),
                           ylim = c(0,20)) +
  theme(plot.title = element_text(size = 20),
        plot.subtitle = element_blank())

volPlot

ggsave('DE_BL/volcanoPlots/GBA+_v_control_All_EVol.jpeg',
       volPlot,
       dpi = 300, width = 6, height = 8)

# Mito volcano plot
keyvals <- rep('gray60', nrow(df))
names(keyvals) <- rep('other', nrow(df))
keyvals[which(df$category == "PD causal variant")] <- '#440154FF'
names(keyvals)[which(df$category == "PD causal variant")] <- 'PD causal variant'
keyvals[which(df$category == "Mitochondria gene")] <- 'darkorange'
names(keyvals)[which(df$category == "Mitochondria gene")] <- 'Mitochondria Genes'
geneLabels <- df[which((df$category == "PD causal variant" | df$category == 'Mitochondria gene')& df$significance != 'insignificant'),]$gene_name
volPlot <- createVolcanoCustom(df, 'GBA+ vs Control - All at Baseline', keyvals, geneLabels, c(-1.1,0.8), c(0,20))
volPlot
ggsave('DE_BL/volcanoPlots/GBA+_v_control_All_EVol_mito.jpeg',
       volPlot,
       dpi = 300, width = 6, height = 8)

# blood cell plot
df_sig <- df[df$adj.P.Val < pval & abs(df$logFC) > logfc ,] 
df_sig$category <- ifelse(df_sig$category %in% c('Lymphocyte gene', 'Neutrophil gene'),  df_sig$category, NA)
plot <- ggplot(df_sig[!is.na(df_sig$category),], aes(x=logFC, fill=category)) +
  geom_histogram(alpha=0.7, bins = 100) +
  xlim(-1.1,0.8) + 
  theme_light() +
  theme(text = element_text(size=15, color='black'),
        axis.text = element_text(size=15, color='black')) 
plot

ggsave('DE_BL/blood_cell_plots/GBA+_v_control_All_blood_genes.jpeg',
       plot,
       dpi = 300, width = 6, height = 3)






########### LRRK2+ v Control #########
df <- read.table(file="DE_BL/limmaResults/LRRK2+_v_control_All_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)

df <- annotate_df(df)
volPlot <- EnhancedVolcano(df,
                           lab = df$gene_name,
                           x = 'logFC',
                           y = 'P.Value',
                           title = 'LRRK2+ vs Control - All at Baseline',
                           ylab = bquote(~-Log[10]~italic(P)),
                           colAlpha = 0.6,
                           pCutoffCol = 'adj.P.Val',
                           pCutoff = 0.05,
                           FCcutoff = 0.1,
                           pointSize = 2.5,
                           labSize = 0,
                           legendPosition = "bottom",
                           xlim = c(-2, 1.4),
                           ylim = c(0,35)) +
  theme(plot.title = element_text(size = 20),
        plot.subtitle = element_blank())

volPlot

ggsave('DE_BL/volcanoPlots/LRRK2+_v_control_All_EVol.jpeg',
       volPlot,
       dpi = 300, width = 6, height = 8)

# Mito volcano plot
keyvals <- rep('gray60', nrow(df))
names(keyvals) <- rep('other', nrow(df))
keyvals[which(df$category == "PD causal variant")] <- '#440154FF'
names(keyvals)[which(df$category == "PD causal variant")] <- 'PD causal variant'
keyvals[which(df$category == "Mitochondria gene")] <- 'darkorange'
names(keyvals)[which(df$category == "Mitochondria gene")] <- 'Mitochondria Genes'
geneLabels <- df[which((df$category == "PD causal variant" | df$category == 'Mitochondria gene')& df$significance != 'insignificant'),]$gene_name
volPlot <- createVolcanoCustom(df, 'LRRK2+ vs Control - All at Baseline', keyvals, geneLabels, c(-2,1.4), c(0,35))
volPlot
ggsave('DE_BL/volcanoPlots/LRRK2+_v_control_All_EVol_mito.jpeg',
       volPlot,
       dpi = 300, width = 6, height = 8)

# blood cell plot
df_sig <- df[df$adj.P.Val < pval & abs(df$logFC) > logfc ,] 
df_sig$category <- ifelse(df_sig$category %in% c('Lymphocyte gene', 'Neutrophil gene'),  df_sig$category, NA)
plot <- ggplot(df_sig[!is.na(df_sig$category),], aes(x=logFC, fill=category)) +
  geom_histogram(alpha=0.7, bins = 100) +
  xlim(-2,1.4) + 
  theme_light() +
  theme(text = element_text(size=15, color='black'),
        axis.text = element_text(size=15, color='black')) 
plot

ggsave('DE_BL/blood_cell_plots/LRRK2+_v_control_All_blood_genes.jpeg',
       plot,
       dpi = 300, width = 6, height = 3)

########### SNCA+ v Control #########
df <- read.table(file="DE_BL/limmaResults/SNCA+_v_control_PPMI_limmaResults.tsv", header=TRUE, sep="\t", check.names=FALSE)

df <- annotate_df(df)
volPlot <- EnhancedVolcano(df,
                           lab = df$gene_name,
                           x = 'logFC',
                           y = 'P.Value',
                           title = 'SNCA+ vs Control - PPMI',
                           ylab = bquote(~-Log[10]~italic(P)),
                           colAlpha = 0.6,
                           pCutoffCol = 'adj.P.Val',
                           pCutoff = 0.05,
                           FCcutoff = 0.1,
                           pointSize = 2.5,
                           labSize = 0,
                           legendPosition = "bottom",
                           xlim = c(-2.9, 2.2),
                           ylim = c(0,15)) +
  theme(plot.title = element_text(size = 20),
        plot.subtitle = element_blank())

volPlot

ggsave('DE_BL/volcanoPlots/SNCA+_v_control_All_EVol.jpeg',
       volPlot,
       dpi = 300, width = 6, height = 8)

# Mito volcano plot
keyvals <- rep('gray60', nrow(df))
names(keyvals) <- rep('other', nrow(df))
keyvals[which(df$category == "PD causal variant")] <- '#440154FF'
names(keyvals)[which(df$category == "PD causal variant")] <- 'PD causal variant'
keyvals[which(df$category == "Mitochondria gene")] <- 'darkorange'
names(keyvals)[which(df$category == "Mitochondria gene")] <- 'Mitochondria Genes'
geneLabels <- df[which((df$category == "PD causal variant" | df$category == 'Mitochondria gene')& df$significance != 'insignificant'),]$gene_name
volPlot <- createVolcanoCustom(df, 'SNCA+ vs Control - All at Baseline', keyvals, geneLabels, c(-2.9,2.2), c(0,15))
volPlot
ggsave('DE_BL/volcanoPlots/SNCA+_v_control_PPMI_EVol_mito.jpeg',
       volPlot,
       dpi = 300, width = 6, height = 8)

# blood cell plot
df_sig <- df[df$adj.P.Val < pval & abs(df$logFC) > logfc ,] 
df_sig$category <- ifelse(df_sig$category %in% c('Lymphocyte gene', 'Neutrophil gene'),  df_sig$category, NA)
plot <- ggplot(df_sig[!is.na(df_sig$category),], aes(x=logFC, fill=category)) +
  geom_histogram(alpha=0.7, bins = 100) +
  xlim(-2.9,2.2) + 
  theme_light() +
  theme(text = element_text(size=15, color='black'),
        axis.text = element_text(size=15, color='black')) 
plot

ggsave('DE_BL/blood_cell_plots/SNCA+_v_control_PPMI_blood_genes.jpeg',
       plot,
       dpi = 300, width = 6, height = 3)







