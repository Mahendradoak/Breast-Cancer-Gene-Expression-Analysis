# ERBB2 (HER2) amplification in TCGA breast cancer

Compares the transcriptome of **ERBB2-amplified versus non-amplified** primary breast tumours in TCGA (cBioPortal `brca_tcga_pan_can_atlas_2018`). It then builds a **cross-validated LASSO-Cox expression signature** for overall survival.

## What it does

1. **Loads the cBioPortal study files.** It skips the clinical metadata lines, keeps primary tumours only (barcode type `01`), and matches patients across expression, copy number and clinical data.
2. **Calls ERBB2 status from GISTIC2.** *Amplified* means a high-level amplification (GISTIC = 2), which corresponds to clinical HER2 positivity. Low-level gains (GISTIC = 1) count as not amplified.
3. **Runs differential expression with limma:**
   - limma-trend on log2(RSEM + 1) by default, because cBioPortal RSEM values are already normalised rather than counts;
   - limma-voom if a raw count matrix is supplied with `--counts`.
4. **Reports a positive control.** Up-regulation of the 17q12 amplicon genes (ERBB2, GRB7, PGAP3, STARD3, MIEN1, …) confirms the comparison works. It is co-amplification, not a finding.
5. **Runs GO gene-set enrichment (GSEA)** on genes ranked by the moderated t statistic.
6. **Builds a survival signature without in-sample evaluation.**
   - Outer 5-fold cross-validation, with gene filtering and λ tuning done inside each training fold.
   - Every patient's risk score comes from a model that never saw that patient.
   - Reports the out-of-fold C-index (± SE), a Kaplan–Meier split on out-of-fold risk, and how often each gene is selected across folds.

## Run

```bash
Rscript install_packages.R
# Download brca_tcga_pan_can_atlas_2018 from cBioPortal and unpack it
Rscript run_analysis.R --study brca_tcga_pan_can_atlas_2018 --out results
Rscript tests/testthat.R      # 10 expectations on a simulated study, ~20 s
```

Outputs in `results/`:

| File | Contents |
|---|---|
| `DE_ERBB2_amplified_vs_not.csv` | limma results for every expressed gene |
| `positive_control_17q12_amplicon.csv` | amplicon genes, which should be strongly up |
| `GSEA_GO_BP.csv` | GO biological process enrichment |
| `volcano.pdf` | volcano plot with amplicon genes highlighted |
| `cv_risk_scores.csv` | out-of-fold risk score and fold for each patient |
| `lasso_selection_frequency.csv` | how often each gene was selected across folds |
| `km_cross_validated_risk.pdf` | Kaplan–Meier split on out-of-fold risk |
| `summary.json` | patients, amplified count, DE count, positive-control check, CV C-index, log-rank p |

## What was wrong before (`archive/`)

| Problem | Effect | Fix |
|---|---|---|
| ERBB2 "Amplified" = GISTIC > 0 | Single-copy gains counted as HER2-amplified, diluting the comparison | GISTIC ≥ 2 |
| DESeq2 on rounded RSEM-normalised values | Violates DESeq2's raw-count negative-binomial model; normalisation applied twice | limma-trend (or limma-voom on real counts) |
| LASSO-Cox fitted on all patients, then the same patients split at the median fitted risk | Kaplan–Meier separation and the 8-gene signature (`archive/original_results/lasso_selected_genes.csv`) are in-sample and optimistic | Nested cross-validation; out-of-fold C-index |
| GSEA ranked by log fold change | Ignores variance; noisy low-expression genes rank highly | Ranked by moderated t |
| Sample IDs collapsed to patients without checking sample type | Metastatic or normal samples could replace primaries | Primary tumours only |
| In sections 1–3, clinical rows indexed into the *uncleaned* table with indices computed on the cleaned table | Clinical data shifted by the metadata rows. The later survival code re-read the data, so the final plots were not affected | Single, tested loader |
| ERBB2 and 17q12 genes presented as top differential findings | Circular: they are co-amplified with ERBB2 | Reported as a positive control |
| `install.packages()` inside the analysis; hard-coded `C:/Users/...` paths; duplicated sections | Not reproducible | `install_packages.R`, command-line paths, one script |

## Status

Tested on simulated data. **The real cBioPortal run still needs to be done**; it takes a few minutes on a laptop.
