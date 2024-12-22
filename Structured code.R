# TCGA Breast Cancer Analysis Pipeline
# This script analyzes breast cancer data from TCGA Pan-Cancer Atlas

# Load required libraries
library(DESeq2)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(glmnet)
library(org.Hs.eg.db)
library(survival)
library(ggplot2)
library(survminer)

###########################################
# 1. Data Loading and Initial Processing
###########################################

# Define file paths
rna_seq_path <- "C:/Users/mahen/Desktop/bioprinciples assignment/assignment 2/brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt"
clinical_path <- "C:/Users/mahen/Desktop/bioprinciples assignment/assignment 2/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt"
cna_path <- "C:/Users/mahen/Desktop/bioprinciples assignment/assignment 2/brca_tcga_pan_can_atlas_2018/data_cna.txt"

# Read data files with proper parameters
rna_seq <- read.delim(rna_seq_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
clinical <- read.delim(clinical_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cna <- read.delim(cna_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)


# Clean clinical data - remove metadata rows
clinical_clean <- clinical[!grepl("^#", clinical$X.Patient.Identifier), ]
clinical_clean <- clinical_clean[clinical_clean$X.Patient.Identifier != "PATIENT_ID", ]

###########################################
# 2. Patient ID Matching
###########################################

# Function to standardize IDs
standardize_ids_strict <- function(ids) {
  ids <- toupper(ids)
  ids <- gsub("-", ".", ids)  # Replace dashes with dots
  ids <- gsub("[^A-Z0-9.]", ".", ids)
  ids <- gsub("\\.\\d+$", "", ids)
  ids <- gsub("\\.+", ".", ids)
  ids <- sub("\\.$", "", ids)
  return(ids)
}

# Standardize IDs across all datasets
rna_ids <- standardize_ids_strict(colnames(rna_seq)[-1])  # Exclude first column
cna_ids <- standardize_ids_strict(colnames(cna)[-1])      # Exclude first column
clinical_ids <- standardize_ids_strict(clinical_clean$X.Patient.Identifier)

# Find common IDs across all datasets
common_ids <- Reduce(intersect, list(rna_ids, cna_ids, clinical_ids))
cat("Number of common patient IDs:", length(common_ids), "\n")

# Filter datasets to include only common IDs
rna_seq_filtered <- rna_seq[, c(1, match(common_ids, rna_ids) + 1)]
cna_filtered <- cna[, c(1, match(common_ids, cna_ids) + 1)]
clinical_filtered <- clinical[match(common_ids, clinical_ids), ]

###########################################
# 3. ERBB2 Status Determination
###########################################

# Extract ERBB2 CNA levels
erbb2_row <- which(cna_filtered$Hugo_Symbol == "ERBB2")
erbb2_cna <- as.numeric(as.character(cna_filtered[erbb2_row, -1]))

# Create metadata with ERBB2 status
metadata <- data.frame(
  Patient_ID = common_ids,
  ERBB2_Status = ifelse(erbb2_cna > 0, "Amplified", "Not_Amplified")
)

###########################################
# 4. DESeq2 Analysis
###########################################

# Prepare count matrix
rna_matrix <- as.matrix(rna_seq_filtered[,-1])
rownames(rna_matrix) <- rna_seq_filtered[,1]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(rna_matrix),
  colData = metadata,
  design = ~ ERBB2_Status
)

# Filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast=c("ERBB2_Status", "Amplified", "Not_Amplified"))
res_ordered <- res[order(res$padj),]

# Get variance stabilized transformation
vst <- varianceStabilizingTransformation(dds)

###########################################
# 5. Visualization
###########################################

# PCA Plot
pca_data <- plotPCA(vst, intgroup="ERBB2_Status", returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color=ERBB2_Status)) +
  geom_point(size=3) +
  geom_density2d() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  scale_color_manual(values=c("blue", "red")) +
  ggtitle("PCA of Gene Expression by ERBB2 Status")

# Heatmap of top DE genes
significant_genes <- which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)
mat <- assay(vst)[significant_genes[1:50], ]
metadata_factors <- data.frame(ERBB2_Status=colData(vst)$ERBB2_Status)

pheatmap(mat,
         annotation_col=metadata_factors,
         scale="row",
         show_rownames=TRUE,
         cluster_cols=TRUE,
         fontsize_row=8)

###########################################
# 6. Pathway Analysis
###########################################

# Prepare gene list
res_df <- as.data.frame(res_ordered)
symbols <- rownames(res_df)

# Map symbols to Entrez IDs
mapped_genes <- select(org.Hs.eg.db,
                      keys=symbols,
                      columns=c("ENTREZID", "SYMBOL"),
                      keytype="SYMBOL")

# Create ranked gene list
res_df$entrez <- mapped_genes$ENTREZID[match(symbols, mapped_genes$SYMBOL)]
gene_list <- res_df$log2FoldChange
names(gene_list) <- res_df$entrez
gene_list <- sort(gene_list[!is.na(names(gene_list))], decreasing=TRUE)

# Perform GO enrichment
go_enrich <- gseGO(geneList = gene_list,
                   ont = "ALL",
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   verbose = TRUE)

###########################################
# 7. Survival Analysis
###########################################
# Load required libraries
library(DESeq2)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(glmnet)
library(org.Hs.eg.db)
library(survival)
library(ggplot2)
library(reshape2)
library(survminer)

# Define file paths
rna_seq_path <- "C:/Users/mahen/Desktop/bioprinciples assignment/assignment 2/brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt"
clinical_path <- "C:/Users/mahen/Desktop/bioprinciples assignment/assignment 2/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt"
cna_path <- "C:/Users/mahen/Desktop/bioprinciples assignment/assignment 2/brca_tcga_pan_can_atlas_2018/data_cna.txt"

# Read data files with proper parameters
rna_seq <- read.delim(rna_seq_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
clinical <- read.delim(clinical_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cna <- read.delim(cna_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Clean clinical data - remove metadata rows
clinical_clean <- clinical[!grepl("^#", clinical$X.Patient.Identifier), ]
clinical_clean <- clinical_clean[clinical_clean$X.Patient.Identifier != "PATIENT_ID", ]

# Function to standardize IDs
standardize_ids_strict <- function(ids) {
  ids <- toupper(ids)
  ids <- gsub("-", ".", ids)  # Replace dashes with dots
  ids <- gsub("[^A-Z0-9.]", ".", ids)
  ids <- gsub("\\.\\d+$", "", ids)
  ids <- gsub("\\.+", ".", ids)
  ids <- sub("\\.$", "", ids)
  return(ids)
}

# Standardize IDs across all datasets
rna_ids <- standardize_ids_strict(colnames(rna_seq)[-1])  # Exclude first column
cna_ids <- standardize_ids_strict(colnames(cna)[-1])      # Exclude first column
clinical_ids <- standardize_ids_strict(clinical_clean$X.Patient.Identifier)

# Update the original data with standardized IDs
colnames(rna_seq)[-1] <- rna_ids
colnames(cna)[-1] <- cna_ids
clinical_clean$X.Patient.Identifier <- clinical_ids

# Find common IDs
common_ids <- Reduce(intersect, list(rna_ids, cna_ids, clinical_ids))
print(paste("Number of common samples:", length(common_ids)))

# Filter datasets to keep only common samples
rna_seq_filtered <- rna_seq[, c(1, match(common_ids, rna_ids) + 1)]
cna_filtered <- cna[, c(1, match(common_ids, cna_ids) + 1)]
clinical_filtered <- clinical_clean[match(common_ids, clinical_ids), ]

# Create ERBB2 metadata
erbb2_col <- which(cna_filtered$Hugo_Symbol == "ERBB2")
erbb2_cna <- cna_filtered[erbb2_col, -1]
erbb2_cna_numeric <- suppressWarnings(as.numeric(as.character(erbb2_cna)))

metadata <- data.frame(
  Patient_ID = colnames(erbb2_cna),
  ERBB2_Status = ifelse(!is.na(erbb2_cna_numeric) & erbb2_cna_numeric > 0, 
                        "Amplified", "Not_Amplified")
)

# Prepare RNA matrix for DESeq2
rna_matrix <- as.matrix(rna_seq_filtered[,-1])
rownames(rna_matrix) <- rna_seq_filtered[,1]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(rna_matrix),
  colData = metadata,
  design = ~ ERBB2_Status
)

# Filter low count genes and run DESeq2
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("ERBB2_Status", "Amplified", "Not_Amplified"))
res_ordered <- res[order(res$padj),]

# Get VST transformation
vst <- varianceStabilizingTransformation(dds)

# Prepare clinical data for survival analysis
clinical_filtered$Overall.Survival..Months. <- as.numeric(as.character(clinical_filtered$Overall.Survival..Months.))
clinical_filtered$Overall.Survival.Status <- as.character(clinical_filtered$Overall.Survival.Status)

# Filter out missing or invalid survival data
complete_cases <- !is.na(clinical_filtered$Overall.Survival..Months.) & 
  !is.na(clinical_filtered$Overall.Survival.Status) &
  clinical_filtered$Overall.Survival..Months. > 0 &
  clinical_filtered$Overall.Survival.Status %in% c("DECEASED", "LIVING")

clinical_survival <- clinical_filtered[complete_cases, ]

# Clean survival status by removing prefixes
clinical_filtered$Overall.Survival.Status <- gsub("^\\d+:", "", clinical_filtered$Overall.Survival.Status)

# Now let's verify the cleaning worked
print("Cleaned survival status:")
print(table(clinical_filtered$Overall.Survival.Status))

# Rerun the complete cases filter
complete_cases <- !is.na(clinical_filtered$Overall.Survival..Months.) & 
  !is.na(clinical_filtered$Overall.Survival.Status) &
  clinical_filtered$Overall.Survival..Months. > 0 &
  clinical_filtered$Overall.Survival.Status %in% c("DECEASED", "LIVING")

clinical_survival <- clinical_filtered[complete_cases, ]

# Verify we have data now
print("Clinical survival dimensions:")
print(dim(clinical_survival))

# Create survival model with the fixed data
vst_matrix <- assay(vst)
common_survival_samples <- intersect(colnames(vst_matrix), clinical_survival$X.Patient.Identifier)

# Align matrices for survival analysis
vst_matrix_aligned <- vst_matrix[, common_survival_samples]
clinical_aligned <- clinical_survival[match(common_survival_samples, clinical_survival$X.Patient.Identifier), ]

# Create design matrix
X <- t(vst_matrix_aligned)
var_predictors <- apply(X, 2, var) > 0
X <- X[, var_predictors]

# Create survival object
y <- Surv(
  time = clinical_aligned$Overall.Survival..Months.,
  event = ifelse(clinical_aligned$Overall.Survival.Status == "DECEASED", 1, 0)
)

# Verify dimensions before model fitting
print("Final dimensions:")
print(paste("X dimensions:", dim(X)[1], "samples x", dim(X)[2], "genes"))
print(paste("Y dimensions:", length(y)))
print(paste("Number of events (deaths):", sum(y[,2])))

# Fit survival model
if(nrow(X) > 0 && ncol(X) > 0 && !any(is.na(y))) {
  fit <- cv.glmnet(X, y, family = "cox", alpha = 1, standardize = TRUE)
  
  # Calculate risk scores
  coef_matrix <- coef(fit, s = "lambda.min")
  risk_scores <- as.numeric(X %*% coef_matrix)
  
  # Add risk scores to clinical data
  clinical_aligned$risk_score <- risk_scores
  clinical_aligned$risk_group <- ifelse(risk_scores > median(risk_scores), 
                                        "High Risk", "Low Risk")
  
  # Create survival curves
  surv_obj <- Surv(time = clinical_aligned$Overall.Survival..Months.,
                   event = clinical_aligned$Overall.Survival.Status == "DECEASED")
  
  fit_surv <- survfit(surv_obj ~ risk_group, data = clinical_aligned)
  
  # Plot survival curves
  pdf("survival_curves.pdf")
  print(ggsurvplot(fit_surv,
                   data = clinical_aligned,
                   pval = TRUE,
                   conf.int = TRUE,
                   risk.table = TRUE,
                   title = "Survival Analysis by Risk Group"))
  dev.off()
  
  # Print some model summary statistics
  print("Number of genes with non-zero coefficients:")
  print(sum(coef_matrix != 0))
  
  # Distribution of risk groups
  print("Risk group distribution:")
  print(table(clinical_aligned$risk_group))
}



# Create enhanced survival plot
# First, create the survival object and fit
surv_obj <- Surv(time = clinical_aligned$Overall.Survival..Months.,
                 event = clinical_aligned$Overall.Survival.Status == "DECEASED")

fit_surv <- survfit(surv_obj ~ risk_group, data = clinical_aligned)

# Create enhanced survival plot with multiple components
surv_plot <- ggsurvplot(
  fit = fit_surv,
  data = clinical_aligned,
  
  # Main plot customization
  pval = TRUE,              # Add p-value
  pval.method = TRUE,       # Add method name
  conf.int = TRUE,          # Add confidence intervals
  conf.int.alpha = 0.2,     # Confidence interval transparency
  
  # Risk table customization
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Color table by groups
  risk.table.height = 0.25, # Table height relative to plot
  
  # Add number at risk table
  tables.height = 0.3,      # Height of tables
  tables.theme = theme_cleantable(),
  
  # Customize lines
  linetype = "strata",      # Different line types
  surv.median.line = "hv",  # Add median survival lines
  
  # Labels and title
  title = "Overall Survival Analysis by Risk Group",
  subtitle = "Based on Gene Expression Signature",
  xlab = "Time (Months)",
  ylab = "Overall Survival Probability",
  legend.title = "Risk Groups",
  legend.labs = c("High Risk", "Low Risk"),
  
  # Color scheme
  palette = c("#E7B800", "#2E9FDF"),
  
  # Theme customization
  ggtheme = theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(face = "italic", size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 10)
    ),
  
  # Additional features
  ncensor.plot = TRUE,      # Add censoring plot
  ncensor.plot.height = 0.25
)

# Add statistical summary to the plot
surv_plot$plot <- surv_plot$plot +
  annotate(
    "text",
    x = max(clinical_aligned$Overall.Survival..Months.) * 0.2,
    y = 0.2,
    label = sprintf(
      "n = %d\nEvents = %d\nMedian Survival (High Risk) = %.1f months\nMedian Survival (Low Risk) = %.1f months",
      nrow(clinical_aligned),
      sum(clinical_aligned$Overall.Survival.Status == "DECEASED"),
      median_surv[1, "median"],
      median_surv[2, "median"]
    ),
    hjust = 0,
    size = 3.5
  )

# Save the plot
pdf("survival_analysis.pdf", width = 10, height = 12)
print(surv_plot)
dev.off()

# Also display the plot in the current device
print(surv_plot)

