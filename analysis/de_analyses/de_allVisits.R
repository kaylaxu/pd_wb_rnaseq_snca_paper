### Differential Expression Analyses for All Visits w/ Neutrophil Correction ###


library(ggplot2) # for plotting
library(ggpubr) # calculate/plot stat significance in boxplots
library(tidyverse) #basis for data manipulation, best practices
library(edgeR) #expression normalization
library(limma) #differential expression

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

tmpG <- data.table::fread("../AMP_RNAseq/GRCh38_GENCODE29_geneInfo.txt", skip = 1)
genes.anno <- tmpG[,c("gene_id", "gene_name", "gene_type", "seqnames", "start", "end", "width", "strand")]
genes.anno <- genes.anno[with(genes.anno, order(gene_id, -width)), ]
head(genes.anno, 12)
genes.anno <- distinct(genes.anno, gene_id, .keep_all = TRUE)


### All case vs all control
# PPMI
setStudy <- c('PPMI')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" )
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Control")

de_analysis <- function(cases, controls, counts, genes, min = 10) {
  cases$group <- 'Case'
  controls$group <- 'Control'
  all_samples <- rbind(cases, controls)
  countSub <- counts[, all_samples$sample_id]
  all(all_samples$sample_id == colnames(countSub))
  
  neutPerc <- all_samples$neutPer
  case <- factor(gsub("\\+", "", all_samples$group), levels = c("Control", "Case"))
  sex <- factor(all_samples$sex,
                levels = c("Female", "Male"))
  mrnaBases <- all_samples$PCT_MRNA_BASES
  ageSquared <- all_samples$ageSquared
  
  design <- model.matrix(~ 0+case + sex  + neutPerc + ageSquared + mrnaBases)
  
  colnames(design) <- gsub("case", "", colnames(design))
  
  contr.matrix <- makeContrasts(
    Case_vs_Control = Case - Control,
    levels = colnames(design)
  )
  
  keep <- filterByExpr(countSub, group = case, min.count = min)
  sum(keep) #this is the number of genes we are testing for differential expression, aka the background set 
  
  dge <- DGEList(countSub[keep,])
  dge  <- calcNormFactors(dge)
  
  genes.anno.Match <- genes[genes$gene_id %in% sort(row.names(dge$counts)),]
  dim(genes.anno.Match)
  i <- match(row.names(dge$counts), genes.anno.Match$gene_id, nomatch = 0)
  genes.anno.Match <- genes.anno.Match[i,]
  identical(genes.anno.Match$gene_id, row.names(dge$counts)) #double checking
  dge$genes <- genes.anno.Match
  
  v <- voom(dge, design, plot=TRUE) 
  vfit <- lmFit(v, design)
  vfit.con1 <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit.con1)
  
  df <- topTable(efit, coef = 1, n = Inf)
  return(df)
}

df <- de_analysis(cases, controls, countTable, genes.anno)
write_tsv(df, "DE_allVisits/limmaResults/case_v_control_PPMI_limmaResults.tsv")

# PDBP
setStudy <- c('PDBP')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" )
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Control" )

df <- de_analysis(cases, controls, countTable, genes.anno)
write_tsv(df, "DE_allVisits/limmaResults/case_v_control_PDBP_limmaResults.tsv")

# All
setStudy <- c('PPMI', 'PDBP')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" )
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Control")

df <- de_analysis(cases, controls, countTable, genes.anno)
write_tsv(df, "DE_allVisits/limmaResults/case_v_control_All_limmaResults.tsv")


### IPD vs control (no mutation)
# PPMI
setStudy <- c('PPMI')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_PD_mutation_in_WGS == 'No')
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Control" & has_known_PD_mutation_in_WGS == 'No')

df <- de_analysis(cases, controls, countTable, genes.anno, min = 20)
write_tsv(df, "DE_allVisits/limmaResults/ipd_v_control_PPMI_limmaResults.tsv")

# PDBP
setStudy <- c('PDBP')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_PD_mutation_in_WGS == 'No')
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Control" & has_known_PD_mutation_in_WGS == 'No')

df <- de_analysis(cases, controls, countTable, genes.anno)
write_tsv(df, "DE_allVisits/limmaResults/ipd_v_control_PDBP_limmaResults.tsv")

# All
setStudy <- c('PPMI', 'PDBP')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_PD_mutation_in_WGS == 'No')
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Control" & has_known_PD_mutation_in_WGS == 'No')

df <- de_analysis(cases, controls, countTable, genes.anno)
write_tsv(df, "DE_allVisits/limmaResults/ipd_v_control_All_limmaResults.tsv")


### GBA+ vs control (no mutation)
# PPMI
setStudy <- c('PPMI')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_GBA_mutation_in_WGS == 'Yes')
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Control" & has_known_PD_mutation_in_WGS == 'No')

df <- de_analysis(cases, controls, countTable, genes.anno, min = 20)
write_tsv(df, "DE_allVisits/limmaResults/GBA+_v_control_PPMI_limmaResults.tsv")

# PDBP
setStudy <- c('PDBP')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_GBA_mutation_in_WGS == 'Yes')
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Control" & has_known_PD_mutation_in_WGS == 'No')

df <- de_analysis(cases, controls, countTable, genes.anno)
write_tsv(df, "DE_allVisits/limmaResults/GBA+_v_control_PDBP_limmaResults.tsv")

# All
setStudy <- c('PPMI', 'PDBP')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_GBA_mutation_in_WGS == 'Yes')
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Control" & has_known_PD_mutation_in_WGS == 'No')

df <- de_analysis(cases, controls, countTable, genes.anno, min = 20)
write_tsv(df, "DE_allVisits/limmaResults/GBA+_v_control_All_limmaResults.tsv")


### LRRK2+ vs control (no mutation)
# PPMI
setStudy <- c('PPMI')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_LRRK2_mutation_in_WGS == 'Yes')
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Control" & has_known_PD_mutation_in_WGS == 'No')

df <- de_analysis(cases, controls, countTable, genes.anno)
write_tsv(df, "DE_allVisits/limmaResults/LRRK2+_v_control_PPMI_limmaResults.tsv")

# PDBP
setStudy <- c('PDBP')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_LRRK2_mutation_in_WGS == 'Yes')
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Control" & has_known_PD_mutation_in_WGS == 'No')

df <- de_analysis(cases, controls, countTable, genes.anno, min = 80)
write_tsv(df, "DE_allVisits/limmaResults/LRRK2+_v_control_PDBP_limmaResults.tsv")

# All
setStudy <- c('PPMI', 'PDBP')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_LRRK2_mutation_in_WGS == 'Yes')
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Control" & has_known_PD_mutation_in_WGS == 'No')

df <- de_analysis(cases, controls, countTable, genes.anno, min = 20)
write_tsv(df, "DE_allVisits/limmaResults/LRRK2+_v_control_All_limmaResults.tsv")


### LRRK2+ vs LRRK2-
# PPMI
setStudy <- c('PPMI')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_LRRK2_mutation_in_WGS == 'Yes')
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_LRRK2_mutation_in_WGS == 'No')

df <- de_analysis(cases, controls, countTable, genes.anno, min = 40)
write_tsv(df, "DE_allVisits/limmaResults/LRRK2+_v_LRRK2-_PPMI_limmaResults.tsv")

# PDBP
setStudy <- c('PDBP')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_LRRK2_mutation_in_WGS == 'Yes')
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_LRRK2_mutation_in_WGS == 'No')

df <- de_analysis(cases, controls, countTable, genes.anno, min = 100)
write_tsv(df, "DE_allVisits/limmaResults/LRRK2+_v_LRRK2-_PDBP_limmaResults.tsv")

# All
setStudy <- c('PPMI', 'PDBP')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_LRRK2_mutation_in_WGS == 'Yes')
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_LRRK2_mutation_in_WGS == 'No')

df <- de_analysis(cases, controls, countTable, genes.anno, min = 45)
write_tsv(df, "DE_allVisits/limmaResults/LRRK2+_v_LRRK2-_All_limmaResults.tsv")


### SNCA+ vs control (no mutation)
# PPMI
setStudy <- c('PPMI')
cases <- sample_meta  %>% filter(study %in% setStudy & case_control_other_at_baseline == "Case" & has_known_SNCA_mutation_in_WGS == 'Yes')
controls <- sample_meta %>% filter(study %in% setStudy & case_control_other_at_baseline == "Control" & has_known_PD_mutation_in_WGS == 'No')

df <- de_analysis(cases, controls, countTable, genes.anno, min = 150)
write_tsv(df, "DE_allVisits/limmaResults/SNCA+_v_control_PPMI_limmaResults.tsv")




