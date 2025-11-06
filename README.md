# TCGA Breast Cancer Gene Expression Analysis

[![R](https://img.shields.io/badge/R-4.0%2B-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive R-based bioinformatics pipeline for analyzing differential gene expression in HER2/ERBB2-amplified breast cancer using TCGA Pan-Cancer Atlas RNA-seq data.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Data Requirements](#data-requirements)
- [Usage](#usage)
- [Pipeline Components](#pipeline-components)
- [Output Files](#output-files)
- [Results](#results)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)
- [Contributing](#contributing)

## Overview

This project provides an end-to-end analysis pipeline for investigating gene expression patterns in HER2+ breast cancer. The pipeline integrates multiple bioinformatics approaches including differential expression analysis, pathway enrichment, and survival prediction using machine learning.

**Key Objectives:**
- Identify differentially expressed genes in ERBB2-amplified vs. non-amplified breast tumors
- Characterize biological pathways associated with HER2+ breast cancer
- Build predictive survival models using gene expression signatures
- Generate publication-ready visualizations

## Features

- **Differential Expression Analysis**: DESeq2-based identification of significant genes
- **Multi-dimensional Visualization**: PCA plots, heatmaps, volcano plots
- **Pathway Analysis**: GO enrichment and GSEA using clusterProfiler
- **Survival Analysis**: LASSO-regularized Cox regression for risk stratification
- **Comprehensive QC**: Data validation and sample matching across datasets
- **Reproducible**: Fully documented code with standardized workflows

## Prerequisites

### System Requirements
- R version ≥ 4.0.0
- Minimum 8GB RAM (16GB recommended for large datasets)
- Operating System: Windows, macOS, or Linux

### Required R Packages

#### Bioconductor Packages
```R
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c(
  "DESeq2",            # Differential expression analysis
  "clusterProfiler",   # Pathway and GO enrichment analysis
  "pheatmap",          # Heatmap visualization
  "ComplexHeatmap",    # Advanced heatmap features
  "org.Hs.eg.db"       # Human gene annotation database
))
```

#### CRAN Packages
```R
install.packages(c(
  "glmnet",          # LASSO and elastic net regression
  "survival",        # Survival analysis
  "survminer",       # Survival plot visualization
  "ggplot2",         # Data visualization
  "reshape2",        # Data manipulation
  "dplyr",           # Data manipulation
  "gridExtra",       # Multiple plot arrangement
  "Cairo"            # High-quality graphics output
))
```

## Installation

1. **Clone the repository**
```bash
git clone https://github.com/yourusername/Breast-Cancer-Gene-Expression-Analysis.git
cd Breast-Cancer-Gene-Expression-Analysis
```

2. **Install dependencies**
```R
# Run the installation script (if provided)
source("install_packages.R")

# Or install packages manually using the commands above
```

3. **Verify installation**
```R
# Check that all packages load correctly
library(DESeq2)
library(clusterProfiler)
library(survival)
```

## Data Requirements

The pipeline requires three input files from the TCGA Breast Cancer (BRCA) Pan-Cancer Atlas dataset:

1. **RNA-seq data**: `data_mrna_seq_v2_rsem.txt`
   - Gene expression values (RSEM normalized)
   - Format: Genes × Samples matrix

2. **Clinical data**: `data_clinical_patient.txt`
   - Patient clinical information including survival data
   - Required fields: Patient ID, Overall Survival (Months), Overall Survival Status

3. **Copy Number Alteration data**: `data_cna.txt`
   - Gene-level CNA data
   - Used to determine ERBB2 amplification status

### Data Download

Download the TCGA BRCA dataset from:
- **cBioPortal**: [Breast Invasive Carcinoma (TCGA, PanCancer Atlas)](https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018)
- Download the complete dataset and extract the three required files

**Important**: Update the file paths in `Structured code.R` (lines 20-22) to point to your data location:
```R
rna_seq_path <- "path/to/your/data_mrna_seq_v2_rsem.txt"
clinical_path <- "path/to/your/data_clinical_patient.txt"
cna_path <- "path/to/your/data_cna.txt"
```

## Usage

### Quick Start

```R
# Load the main analysis script
source("Structured code.R")
```

The script will automatically:
1. Load and preprocess data
2. Perform differential expression analysis
3. Generate visualizations
4. Conduct pathway analysis
5. Build survival models
6. Save all outputs

### Step-by-Step Execution

For more control, you can run individual sections:

```R
# 1. Data Loading
source("Structured code.R", echo=TRUE, max.deparse.length=1000)
# Stop after section 1 for QC

# 2. Run DESeq2 Analysis
# Continue with sections 2-4

# 3. Generate Visualizations
# Sections 5

# 4. Pathway Analysis
# Section 6

# 5. Survival Analysis
# Section 7
```

## Pipeline Components

### 1. Data Loading and Preprocessing
- ID standardization across datasets
- Removal of metadata rows
- Sample matching and filtering
- Data quality checks

**Key Function**: `standardize_ids_strict()`
- Ensures consistent patient ID format across all datasets
- Handles various ID formats from TCGA

### 2. ERBB2 Status Classification
Samples are classified based on ERBB2 copy number alteration:
- **Amplified**: CNA > 0
- **Not Amplified**: CNA ≤ 0

### 3. Differential Expression Analysis
Uses DESeq2 methodology:
- Normalization for sequencing depth and composition
- Negative binomial modeling
- Variance stabilizing transformation (VST)
- Multiple testing correction (Benjamini-Hochberg)

**Significance Thresholds**:
- Adjusted p-value < 0.05
- |Log2 Fold Change| > 1

### 4. Visualization Suite

#### PCA Plot
- Dimensionality reduction to visualize sample clustering
- Color-coded by ERBB2 status
- Includes variance explained by each PC

#### Heatmap
- Top 50 differentially expressed genes
- Hierarchical clustering of samples and genes
- Z-score normalized expression values

#### Volcano Plot
- Visualizes significance vs. fold change
- Highlights significant genes (padj < 0.05, |log2FC| > 1)

#### Top DE Genes Bar Plot
- Top 10 most significant genes
- Directional change (up/down regulated)

### 5. Pathway Analysis
GO enrichment using clusterProfiler:
- Gene Set Enrichment Analysis (GSEA)
- Biological Process, Molecular Function, Cellular Component
- Network visualization of enriched pathways

### 6. Survival Analysis

#### LASSO Cox Regression
- Feature selection from thousands of genes
- 10-fold cross-validation for lambda optimization
- Risk score calculation
- Patient stratification into High/Low risk groups

**Outputs**:
- Survival curves with log-rank test
- Risk tables
- Gene coefficients
- Cross-validation plot

## Output Files

The pipeline generates the following output files:

### Differential Expression
- `differential_expression_analysis.pdf` - Combined DE plots
- `top_10_de_genes.pdf` - Top 10 DE genes bar plot
- `volcano_plot.pdf` - Volcano plot of all genes
- `combined_plots.pdf` - Side-by-side comparison

### Visualization
- `PCA_plot.pdf` - Principal component analysis
- `heatmap.pdf` - Clustered heatmap of top DE genes

### Pathway Analysis
- `go_enrichment_dotplot.pdf` - GO enrichment results
- `go_enrichment_gene_network.pdf` - Gene-pathway network

### Survival Analysis
- `survival_analysis.pdf` - Basic survival curves
- `lasso_survival_analysis.pdf` - LASSO-based survival curves
- `lasso_coefficients.pdf` - Selected gene coefficients
- `lasso_cv_plot.pdf` - Cross-validation plot
- `lasso_selected_genes.csv` - List of genes selected by LASSO

## Results

### Expected Outcomes

**Differential Expression**:
- Hundreds to thousands of significantly differentially expressed genes
- Strong enrichment of known HER2-pathway genes (ERBB2, GRB7, etc.)

**Pathway Analysis**:
- Enrichment in growth factor signaling pathways
- Cell proliferation and migration processes
- Immune response pathways

**Survival Analysis**:
- Significant stratification of patients into risk groups
- Gene signature predictive of overall survival
- Typical hazard ratios > 2 for high vs. low risk

### Quality Metrics
- PCA should show separation of ERBB2+ and ERBB2- samples
- Top DE genes should include ERBB2 and neighboring genes
- P-values for survival analysis typically < 0.05

## Troubleshooting

### Common Issues

**Issue 1: Sample ID Mismatches**
```
Error: Samples don't match across datasets
```
**Solution**: The `standardize_ids_strict()` function should handle this. Verify your data files are from the same TCGA cohort.

**Issue 2: Memory Issues**
```
Error: cannot allocate vector of size...
```
**Solution**:
- Filter low-count genes earlier (line 93-94)
- Reduce number of genes in heatmap
- Increase R memory limit: `memory.limit(size=16000)`

**Issue 3: Missing Survival Data**
```
Warning: Removed X rows with missing survival data
```
**Solution**: This is normal. The script filters out incomplete cases automatically. Check that you have >50 complete cases for survival analysis.

**Issue 4: File Path Errors**
```
Error: cannot open file 'data_mrna_seq_v2_rsem.txt'
```
**Solution**: Update file paths in lines 20-22 to match your data location. Use absolute paths for clarity.

**Issue 5: Package Installation Failures**
```
Installation of package 'X' had non-zero exit status
```
**Solution**:
- Ensure R version ≥ 4.0
- Try installing packages individually
- For Bioconductor packages, ensure BiocManager is up to date:
  ```R
  BiocManager::install(version = "3.18")
  ```

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{breast_cancer_analysis,
  title = {TCGA Breast Cancer Gene Expression Analysis Pipeline},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/yourusername/Breast-Cancer-Gene-Expression-Analysis}
}
```

### Data Citation
Please also cite the TCGA Pan-Cancer Atlas:

```bibtex
@article{tcga2018,
  title={The Cancer Genome Atlas Pan-Cancer analysis project},
  author={Cancer Genome Atlas Research Network and others},
  journal={Nature genetics},
  volume={45},
  number={10},
  pages={1113--1120},
  year={2013},
  publisher={Nature Publishing Group}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes:

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

### Areas for Contribution
- Additional visualization methods
- Integration of other cancer types
- Performance optimizations
- Additional pathway databases
- Docker containerization
- Automated testing

## Acknowledgments

- **TCGA Research Network** for providing the data
- **Bioconductor** community for excellent R packages
- **cBioPortal** for data distribution and visualization tools

## Contact

For questions, issues, or suggestions, please:
- Open an issue on GitHub
- Contact: [Your Email]

---

**Note**: This pipeline is for research purposes only. Results should be validated before clinical application.

Last Updated: November 2024
