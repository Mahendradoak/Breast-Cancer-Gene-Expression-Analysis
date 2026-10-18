# ERBB2 (HER2) amplification in TCGA breast cancer: expression consequences and
# a cross-validated survival signature. Inputs are the cBioPortal
# "brca_tcga_pan_can_atlas_2018" study files.

# ---- loading -----------------------------------------------------------------
read_cbio_matrix <- function(path) {
  d <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  d <- d[!is.na(d$Hugo_Symbol) & d$Hugo_Symbol != "" & !duplicated(d$Hugo_Symbol), ]
  m <- as.matrix(d[, setdiff(names(d), c("Hugo_Symbol", "Entrez_Gene_Id")), drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- d$Hugo_Symbol
  m
}

# cBioPortal clinical files start with four '#' metadata lines; read past them.
read_cbio_clinical <- function(path) {
  d <- utils::read.delim(path, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  d
}

# Keep primary tumours (barcode type 01) and key every matrix by 12-char patient ID.
# The original collapsed sample IDs to patient IDs without checking sample type.
to_patient_columns <- function(m) {
  st <- substr(colnames(m), 14, 15)
  m <- m[, st == "01", drop = FALSE]
  colnames(m) <- substr(colnames(m), 1, 12)
  m[, !duplicated(colnames(m)), drop = FALSE]
}

# ---- ERBB2 status ---------------------------------------------------------------
# GISTIC2 thresholded calls: -2 deep deletion, -1 shallow loss, 0 neutral,
# 1 low-level gain, 2 high-level amplification. Clinical HER2 positivity
# corresponds to high-level amplification (2), not any gain (> 0) as before.
erbb2_status <- function(cna, threshold = 2) {
  v <- cna["ERBB2", ]
  factor(ifelse(v >= threshold, "Amplified", "Not_Amplified"), levels = c("Not_Amplified", "Amplified"))
}

# ---- differential expression -------------------------------------------------
# cBioPortal RSEM values are already normalised (not integer counts), so DESeq2's
# negative-binomial model does not apply. Use limma-trend on log2(RSEM + 1).
# If raw counts are supplied (e.g. GDC STAR counts), use limma-voom instead.
de_limma <- function(expr, group, counts = FALSE, min_expr = 1) {
  design <- stats::model.matrix(~ group)
  if (counts) {
    dge <- edgeR::DGEList(expr)
    keep <- edgeR::filterByExpr(dge, design)
    dge <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
    v <- limma::voom(dge, design)
    fit <- limma::eBayes(limma::lmFit(v, design))
  } else {
    y <- log2(expr + 1)
    y <- y[rowMeans(y) >= min_expr, , drop = FALSE]
    fit <- limma::eBayes(limma::lmFit(y, design), trend = TRUE)
  }
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")
  data.frame(gene = rownames(tt), log2FC = tt$logFC, t = tt$t, p_value = tt$P.Value, adj_p_value = tt$adj.P.Val)
}

# Genes of the 17q12 ERBB2 amplicon. Their up-regulation is a POSITIVE CONTROL
# (co-amplification), not a discovery.
AMPLICON_17Q12 <- c("ERBB2", "GRB7", "PGAP3", "STARD3", "MIEN1", "PNMT", "TCAP", "PPP1R1B", "NEUROD2", "IKZF3")

# ---- GSEA ----------------------------------------------------------------------
# Rank by moderated t (accounts for variance), not by raw log fold change.
gsea_go <- function(de) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message("clusterProfiler/org.Hs.eg.db not installed; skipping GSEA")
    return(NULL)
  }
  map <- suppressMessages(clusterProfiler::bitr(de$gene, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db))
  r <- de$t[match(map$SYMBOL, de$gene)]
  names(r) <- map$ENTREZID
  r <- sort(r[!duplicated(names(r))], decreasing = TRUE)
  as.data.frame(clusterProfiler::gseGO(r, ont = "BP", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID",
                                       minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05, seed = TRUE, verbose = FALSE))
}

# ---- survival signature --------------------------------------------------------
parse_os <- function(clin) {
  status <- sub("^\\d+:", "", clin$OS_STATUS)
  data.frame(patient_id = clin$PATIENT_ID, time = suppressWarnings(as.numeric(clin$OS_MONTHS)),
             event = as.integer(status == "DECEASED"))
}

# The original fitted LASSO-Cox on ALL patients and then split the SAME patients
# at the median of the fitted risk score, so its Kaplan-Meier separation was
# in-sample. Here every patient's risk score comes from a model that never saw
# them (outer K-fold), with lambda tuned inside each training fold.
cv_lasso_cox <- function(x, time, event, outer_folds = 5, n_genes = 2000, seed = 42) {
  set.seed(seed)
  fold <- sample(rep(seq_len(outer_folds), length.out = nrow(x)))
  lp <- rep(NA_real_, nrow(x))
  selected <- list()
  for (k in seq_len(outer_folds)) {
    tr <- fold != k
    v <- apply(x[tr, , drop = FALSE], 2, stats::var)
    keep <- order(v, decreasing = TRUE)[seq_len(min(n_genes, sum(v > 0)))]   # filter on training fold only
    y <- survival::Surv(time[tr], event[tr])
    fit <- glmnet::cv.glmnet(x[tr, keep, drop = FALSE], y, family = "cox", alpha = 1, nfolds = 10)
    lp[!tr] <- as.numeric(stats::predict(fit, x[!tr, keep, drop = FALSE], s = "lambda.min", type = "link"))
    b <- as.matrix(stats::coef(fit, s = "lambda.min"))
    selected[[k]] <- rownames(b)[b[, 1] != 0]
  }
  cidx <- survival::concordance(survival::Surv(time, event) ~ lp, reverse = TRUE)
  freq <- sort(table(unlist(selected)), decreasing = TRUE)
  list(linear_predictor = lp, fold = fold, c_index = unname(cidx$concordance), c_index_se = sqrt(cidx$var),
       selection_frequency = data.frame(gene = names(freq), folds_selected = as.integer(freq)))
}
