# TCGA Breast Cancer Gene Expression Analysis Pipeline

This repository contains an R-based analysis pipeline for investigating differential gene expression in HER2/ERBB2+ breast cancer using TCGA RNA-seq data.

## Prerequisites

### Required R Packages
```R
# Install BiocManager if not already installed
install.packages("BiocManager")

# Install required packages
BiocManager::install(c(
  "DESeq2",            # For differential expression analysis
  "clusterProfiler",   # For pathway analysis
  "pheatmap",          # For heatmap visualization
  "ComplexHeatmap",    # For advanced heatmap visualization
  "glmnet",           # For LASSO regression
  "org.Hs.eg.db"      # For gene annotation
))

# Additional CRAN packages
install.packages(c(
  "survival",        # For survival analysis
  "survminer",       # For survival visualization
  "ggplot2",         # For plotting
  "reshape2"         # For data manipulation
))
```

## Data Requirements

The pipeline expects three input files from the TCGA breast cancer dataset:
1. RNA-seq data (`data_mrna_seq_v2_rsem.txt`)
2. Clinical data (`data_clinical_patient.txt`)
3. Copy Number Alteration (CNA) data (`data_cna.txt`)

## Pipeline Components

### 1. Data Loading and Preprocessing
```R
# Load data files
rna_seq <- read.delim("data_mrna_seq_v2_rsem.txt", sep="\t", header=TRUE)
clinical <- read.delim("data_clinical_patient.txt", sep="\t", header=TRUE)
cna <- read.delim("data_cna.txt", sep="\t", header=TRUE)
```

Key preprocessing steps:
- ID standardization across datasets
- Removal of metadata rows from clinical data
- Handling of missing values
- Sample matching across datasets

### 2. ERBB2 Status Classification
The pipeline classifies samples based on ERBB2 amplification status:
- Amplified: CNA > 0
- Not Amplified: CNA ≤ 0

### 3. Differential Expression Analysis
Using DESeq2 for:
- Data normalization
- Differential expression testing
- Variance stabilizing transformation (VST)

### 4. Visualization Components

#### PCA Plot
```R
# Generate PCA plot
pca_data <- plotPCA(vst, intgroup="ERBB2_Status", returnData=TRUE)
ggplot(pca_data, aes(PC1, PC2, color=ERBB2_Status)) +
  geom_point(size=3) +
  geom_density2d()
```

#### Heatmap Generation
```R
# Generate heatmap of top DE genes
pheatmap(mat, 
         annotation_col=metadata_factors,
         scale="row",
         show_rownames=TRUE)
```

### 5. Pathway Analysis
Using clusterProfiler for:
- GO enrichment analysis
- GSEA analysis
- Pathway visualization

### 6. Survival Analysis
Implements LASSO-regularized Cox regression:
- Patient stratification
- Survival curve generation
- Risk score calculation

## Function Descriptions

### ID Standardization
```R
standardize_ids_strict <- function(ids) {
  ids <- toupper(ids)                      
  ids <- gsub("[^A-Z0-9.]", ".", ids)     
  ids <- gsub("\\.\\d+$", "", ids)        
  ids <- gsub("\\.+", ".", ids)           
  ids <- sub("\\.$", "", ids)             
  return(ids)
}
```

### Data Filtering
- Removes low count genes (< 10 counts)
- Handles missing values in survival data
- Matches samples across datasets

## Output Files

The pipeline generates:
1. Differential expression results
2. PCA plots
3. Heatmaps
4. Pathway enrichment results
5. Survival analysis plots
6. Risk stratification results

## Validation

The pipeline includes validation steps:
- Known ERBB2+ signature genes verification
- Data quality checks
- Sample matching verification
- Survival data completeness checks

## Usage Example

```R
# Load required libraries
source("required_libraries.R")

# Run analysis pipeline
source("main_analysis.R")

# Generate visualizations
source("visualization.R")

# Perform survival analysis
source("survival_analysis.R")
```

## Notes

- Ensure all input files are in the correct format
- Monitor memory usage with large datasets
- Consider using parallel processing for large-scale analyses
- Verify sample IDs match across all input files

## Troubleshooting

Common issues and solutions:
1. Sample ID mismatches: Use the standardize_ids_strict function
2. Memory issues: Filter low-count genes early
3. Missing survival data: Check completeness of clinical data
4. Zero-variance genes: Remove before LASSO regression

## Contributing

Feel free to submit issues and enhancement requests!
