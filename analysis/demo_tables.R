# Table 1, Supplemental Table 1, Supplemental Table 2

metaData <- read.csv("data/amp_meta_neutScos.csv", check.names = FALSE)
sample_meta <- read.csv(file="data/AMP_RNAseq_metaData.csv", header=TRUE)
sample_meta <- sample_meta[sample_meta$PCT_CHIMERAS < 0.03,]
sample_meta <- sample_meta[sample_meta$study != 'BioFIND',]
metaData <- metaData[metaData$sample_id %in% sample_meta$sample_id,] # apply additional filtering

metaData$ageBin <- factor(ifelse(metaData$age_at_baseline < 55, '< 55', ifelse(metaData$age_at_baseline > 65, '> 65', '55 to 65')), levels = c('< 55', '55 to 65', '> 65'))
sample_meta$ageBin <- factor(ifelse(sample_meta$age_at_baseline < 55, '< 55', ifelse(sample_meta$age_at_baseline > 65, '> 65', '55 to 65')), levels = c('< 55', '55 to 65', '> 65'))
sample_meta$genetic_status <- ifelse(is.na(sample_meta$has_known_PD_mutation_in_WGS), 'Unknown',
                                     ifelse(sample_meta$has_known_PD_mutation_in_WGS %in% c('No'), 'None',
                                            ifelse(sample_meta$has_known_SNCA_mutation_in_WGS == 'Yes' & sample_meta$has_known_GBA_mutation_in_WGS == 'No' & sample_meta$has_known_LRRK2_mutation_in_WGS == 'No', 'SNCA+',
                                                   ifelse(sample_meta$has_known_SNCA_mutation_in_WGS == 'No' & sample_meta$has_known_GBA_mutation_in_WGS == 'No' & sample_meta$has_known_LRRK2_mutation_in_WGS == 'Yes', 'LRRK2+',
                                                          ifelse(sample_meta$has_known_SNCA_mutation_in_WGS == 'No' & sample_meta$has_known_GBA_mutation_in_WGS == 'Yes' & sample_meta$has_known_LRRK2_mutation_in_WGS == 'No', 'GBA+', 
                                                                 'Multiple')))))

length(unique(sample_meta$participant_id)) #2776
length(sample_meta$sample_id) # 6897


cbc_merge <- merge(metaData, sample_meta, by='sample_id')
cbc_male_samples <- cbc_merge[cbc_merge$sex.x == 'Male',]
cbc_fem_samples <- cbc_merge[cbc_merge$sex.x == 'Female',]
cbc_male_part <- cbc_male_samples[!duplicated(cbc_male_samples$participant_id.x),]
cbc_fem_part <- cbc_fem_samples[!duplicated(cbc_fem_samples$participant_id.x),]


table(cbc_fem_part[cbc_fem_part$case_control_other_at_baseline.x == 'Control',]$ageBin.x, cbc_fem_part[cbc_fem_part$case_control_other_at_baseline.x == 'Control',]$genetic_status)
table(cbc_male_part[cbc_male_part$case_control_other_at_baseline.x == 'Control',]$ageBin.x, cbc_male_part[cbc_male_part$case_control_other_at_baseline.x == 'Control',]$genetic_status)
table(cbc_fem_samples[cbc_fem_samples$case_control_other_at_baseline.x == 'Control',]$ageBin.x, cbc_fem_samples[cbc_fem_samples$case_control_other_at_baseline.x == 'Control',]$genetic_status)
table(cbc_male_samples[cbc_male_samples$case_control_other_at_baseline.x == 'Control',]$ageBin.x, cbc_male_samples[cbc_male_samples$case_control_other_at_baseline.x == 'Control',]$genetic_status)


table(cbc_fem_part[cbc_fem_part$case_control_other_at_baseline.x == 'Case',]$ageBin.x, cbc_fem_part[cbc_fem_part$case_control_other_at_baseline.x == 'Case',]$genetic_status)
table(cbc_male_part[cbc_male_part$case_control_other_at_baseline.x == 'Case',]$ageBin.x, cbc_male_part[cbc_male_part$case_control_other_at_baseline.x == 'Case',]$genetic_status)
table(cbc_fem_samples[cbc_fem_samples$case_control_other_at_baseline.x == 'Case',]$ageBin.x, cbc_fem_samples[cbc_fem_samples$case_control_other_at_baseline.x == 'Case',]$genetic_status)
table(cbc_male_samples[cbc_male_samples$case_control_other_at_baseline.x == 'Case',]$ageBin.x, cbc_male_samples[cbc_male_samples$case_control_other_at_baseline.x == 'Case',]$genetic_status)


sub_samples <- sample_meta[!(sample_meta$sample_id %in% metaData$sample_id),]
sub_fem_sample <- sub_samples[sub_samples$sex == 'Female',]
sub_male_sample <- sub_samples[sub_samples$sex == 'Male',]
sub_fem_part <- sub_fem_sample[!duplicated(sub_fem_sample$participant_id),]
sub_male_part <- sub_male_sample[!duplicated(sub_male_sample$participant_id),]


table(sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Control',]$ageBin, sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Control',]$genetic_status)
table(sub_male_part[sub_male_part$case_control_other_at_baseline == 'Control',]$ageBin, sub_male_part[sub_male_part$case_control_other_at_baseline == 'Control',]$genetic_status)
table(sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Control',]$ageBin, sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Control',]$genetic_status)
table(sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Control',]$ageBin, sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Control',]$genetic_status)

table(sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Case',]$ageBin, sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Case',]$genetic_status)
table(sub_male_part[sub_male_part$case_control_other_at_baseline == 'Case',]$ageBin, sub_male_part[sub_male_part$case_control_other_at_baseline == 'Case',]$genetic_status)
table(sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Case',]$ageBin, sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Case',]$genetic_status)
table(sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Case',]$ageBin, sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Case',]$genetic_status)

table(sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Other',]$ageBin, sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Other',]$genetic_status)
table(sub_male_part[sub_male_part$case_control_other_at_baseline == 'Other',]$ageBin, sub_male_part[sub_male_part$case_control_other_at_baseline == 'Other',]$genetic_status)
table(sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Other',]$ageBin, sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Other',]$genetic_status)
table(sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Other',]$ageBin, sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Other',]$genetic_status)




#PPMI

table(cbc_fem_part[cbc_fem_part$case_control_other_at_baseline.x == 'Control' & cbc_fem_part$study %in% c('PPMI'),]$ageBin.x, cbc_fem_part[cbc_fem_part$case_control_other_at_baseline.x == 'Control' & cbc_fem_part$study %in% c('PPMI'),]$genetic_status)
table(cbc_male_part[cbc_male_part$case_control_other_at_baseline.x == 'Control' & cbc_male_part$study %in% c('PPMI'),]$ageBin.x, cbc_male_part[cbc_male_part$case_control_other_at_baseline.x == 'Control' & cbc_male_part$study %in% c('PPMI'),]$genetic_status)
table(cbc_fem_samples[cbc_fem_samples$case_control_other_at_baseline.x == 'Control' & cbc_fem_samples$study %in% c('PPMI'),]$ageBin.x, cbc_fem_samples[cbc_fem_samples$case_control_other_at_baseline.x == 'Control' & cbc_fem_samples$study %in% c('PPMI'),]$genetic_status)
table(cbc_male_samples[cbc_male_samples$case_control_other_at_baseline.x == 'Control' & cbc_male_samples$study %in% c('PPMI'),]$ageBin.x, cbc_male_samples[cbc_male_samples$case_control_other_at_baseline.x == 'Control' & cbc_male_samples$study %in% c('PPMI'),]$genetic_status)


table(cbc_fem_part[cbc_fem_part$case_control_other_at_baseline.x == 'Case' & cbc_fem_part$study %in% c('PPMI'),]$ageBin.x, cbc_fem_part[cbc_fem_part$case_control_other_at_baseline.x == 'Case' & cbc_fem_part$study %in% c('PPMI'),]$genetic_status)
table(cbc_male_part[cbc_male_part$case_control_other_at_baseline.x == 'Case'& cbc_male_part$study %in% c('PPMI'),]$ageBin.x, cbc_male_part[cbc_male_part$case_control_other_at_baseline.x == 'Case' & cbc_male_part$study %in% c('PPMI'),]$genetic_status)
table(cbc_fem_samples[cbc_fem_samples$case_control_other_at_baseline.x == 'Case'& cbc_fem_samples$study %in% c('PPMI'),]$ageBin.x, cbc_fem_samples[cbc_fem_samples$case_control_other_at_baseline.x == 'Case'& cbc_fem_samples$study %in% c('PPMI'),]$genetic_status)
table(cbc_male_samples[cbc_male_samples$case_control_other_at_baseline.x == 'Case'& cbc_male_samples$study %in% c('PPMI'),]$ageBin.x, cbc_male_samples[cbc_male_samples$case_control_other_at_baseline.x == 'Case'& cbc_male_samples$study %in% c('PPMI'),]$genetic_status)


table(sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Control' & sub_fem_part$study %in% c('PPMI'),]$ageBin, sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Control' & sub_fem_part$study %in% c('PPMI'),]$genetic_status)
table(sub_male_part[sub_male_part$case_control_other_at_baseline == 'Control' & sub_male_part$study %in% c('PPMI'),]$ageBin, sub_male_part[sub_male_part$case_control_other_at_baseline == 'Control' & sub_male_part$study %in% c('PPMI'),]$genetic_status)
table(sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Control' & sub_fem_sample$study %in% c('PPMI'),]$ageBin, sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Control' & sub_fem_sample$study %in% c('PPMI'),]$genetic_status)
table(sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Control' & sub_male_sample$study %in% c('PPMI'),]$ageBin, sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Control' & sub_male_sample$study %in% c('PPMI'),]$genetic_status)

table(sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Case' & sub_fem_part$study %in% c('PPMI'),]$ageBin, sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Case' & sub_fem_part$study %in% c('PPMI'),]$genetic_status)
table(sub_male_part[sub_male_part$case_control_other_at_baseline == 'Case' & sub_male_part$study %in% c('PPMI'),]$ageBin, sub_male_part[sub_male_part$case_control_other_at_baseline == 'Case' & sub_male_part$study %in% c('PPMI'),]$genetic_status)
table(sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Case' & sub_fem_sample$study %in% c('PPMI'),]$ageBin, sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Case' & sub_fem_sample$study %in% c('PPMI'),]$genetic_status)
table(sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Case' & sub_male_sample$study %in% c('PPMI'),]$ageBin, sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Case' & sub_male_sample$study %in% c('PPMI'),]$genetic_status)

table(sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Other' & sub_fem_part$study %in% c('PPMI'),]$ageBin, sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Other' & sub_fem_part$study %in% c('PPMI'),]$genetic_status)
table(sub_male_part[sub_male_part$case_control_other_at_baseline == 'Other' & sub_male_part$study %in% c('PPMI'),]$ageBin, sub_male_part[sub_male_part$case_control_other_at_baseline == 'Other' & sub_male_part$study %in% c('PPMI'),]$genetic_status)
table(sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Other' & sub_fem_sample$study %in% c('PPMI'),]$ageBin, sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Other' & sub_fem_sample$study %in% c('PPMI'),]$genetic_status)
table(sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Other' & sub_male_sample$study %in% c('PPMI'),]$ageBin, sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Other' & sub_male_sample$study %in% c('PPMI'),]$genetic_status)




# PDBP
table(sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Control' & sub_fem_part$study %in% c('PDBP'),]$ageBin, sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Control' & sub_fem_part$study %in% c('PDBP'),]$genetic_status)
table(sub_male_part[sub_male_part$case_control_other_at_baseline == 'Control' & sub_male_part$study %in% c('PDBP'),]$ageBin, sub_male_part[sub_male_part$case_control_other_at_baseline == 'Control' & sub_male_part$study %in% c('PDBP'),]$genetic_status)
table(sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Control' & sub_fem_sample$study %in% c('PDBP'),]$ageBin, sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Control' & sub_fem_sample$study %in% c('PDBP'),]$genetic_status)
table(sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Control' & sub_male_sample$study %in% c('PDBP'),]$ageBin, sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Control' & sub_male_sample$study %in% c('PDBP'),]$genetic_status)

table(sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Case' & sub_fem_part$study %in% c('PDBP'),]$ageBin, sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Case' & sub_fem_part$study %in% c('PDBP'),]$genetic_status)
table(sub_male_part[sub_male_part$case_control_other_at_baseline == 'Case' & sub_male_part$study %in% c('PDBP'),]$ageBin, sub_male_part[sub_male_part$case_control_other_at_baseline == 'Case' & sub_male_part$study %in% c('PDBP'),]$genetic_status)
table(sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Case' & sub_fem_sample$study %in% c('PDBP'),]$ageBin, sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Case' & sub_fem_sample$study %in% c('PDBP'),]$genetic_status)
table(sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Case' & sub_male_sample$study %in% c('PDBP'),]$ageBin, sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Case' & sub_male_sample$study %in% c('PDBP'),]$genetic_status)

table(sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Other' & sub_fem_part$study %in% c('PDBP'),]$ageBin, sub_fem_part[sub_fem_part$case_control_other_at_baseline == 'Other' & sub_fem_part$study %in% c('PDBP'),]$genetic_status)
table(sub_male_part[sub_male_part$case_control_other_at_baseline == 'Other' & sub_male_part$study %in% c('PDBP'),]$ageBin, sub_male_part[sub_male_part$case_control_other_at_baseline == 'Other' & sub_male_part$study %in% c('PDBP'),]$genetic_status)
table(sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Other' & sub_fem_sample$study %in% c('PDBP'),]$ageBin, sub_fem_sample[sub_fem_sample$case_control_other_at_baseline == 'Other' & sub_fem_sample$study %in% c('PDBP'),]$genetic_status)
table(sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Other' & sub_male_sample$study %in% c('PDBP'),]$ageBin, sub_male_sample[sub_male_sample$case_control_other_at_baseline == 'Other' & sub_male_sample$study %in% c('PDBP'),]$genetic_status)



## Double check duplicate participants in cbc and w/out cbc
sub_samples <- sample_meta[!(sample_meta$sample_id %in% metaData$sample_id) & sample_meta$study %in% c('PPMI'),] #2711 samples, correct
length(unique(sub_samples$participant_id)) #1475 participants that have samples w/out CBC
length(unique( metaData$participant_id)) #600 participants with samples that have CBC

sum(unique(sub_samples$participant_id) %in% unique( metaData$participant_id)) #578 participants with samples in both

# 600 + (1475 - 578) = 1497









