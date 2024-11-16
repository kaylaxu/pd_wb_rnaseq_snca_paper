# Decreased SNCA Expression in Whole-Blood RNA Analysis of Parkinson’s Disease Adjusting for Lymphocytes
Kayla Xu<sup>1</sup>, Ivo Violich<sup>2</sup>, Elizabeth Hutchins<sup>3</sup>, Eric Alsop<sup>3</sup>, Mike A. Nalls<sup>3</sup>, Cornelis Blauwendraat<sup>4</sup>, J. Raphael Gibbs<sup>4</sup>, Mark R. Cookson<sup>4</sup>, Anni Moore<sup>4</sup>, Kendall Van Keuren-Jensen<sup>5</sup>, David W. Craig<sup>1</sup>*
1.	Integrated Translational Sciences, Beckman Research Institute, City of Hope, Duarte, CA 
2.	Department of Translational Genomics, Keck School of Medicine, University of Southern California, Los Angeles, CA 
3.	Laboratory of Neurogenetics, National Institute on Aging, National Institutes of Health, Bethesda, MD, USA; DataTecnica International, Glen Echo, MD, USA. 
4.	Laboratory of Neurogenetics, National Institute on Aging, National Institutes of Health, Bethesda, MD, USA 
5.	Center for Alzheimer’s and Related dementias, National Institutes of Health, Bethesda, MD USA 

## analysis
* de_anlayses/
    * R scripts for differential expression analysis
    * R scripts for volcano plot generation 
* mi_feature_selection.ipynb
    * data-driven feature selection using mutual inforamation scores 
* regression_models.R
    * linear regression and XGBoost models for neutrophil percentage prediction
* snca_plots.R
    * CPM and log normalized SNCA gene counts, corrected for predicted neutrophil expression
    * plotted against PD demongraphic, clinical, and biological factors
* umap.R
    * UMAP analysis of whole blood RNAseq based on significantly expression pathways identified via Ingenuity Pathway Analysis
* variance_anlaysis.R
    * PCA of gene expression and correlation with sample metadata (e.g. clinical, QC, etc.)
 
## figures
Code for figure generation
* figure2.R
* figure4.R
