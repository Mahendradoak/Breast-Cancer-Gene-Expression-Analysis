#!/usr/bin/env Rscript
# Usage:
#   Rscript run_analysis.R --study brca_tcga_pan_can_atlas_2018 --out results
#     (expects data_mrna_seq_v2_rsem.txt, data_cna.txt, data_clinical_patient.txt in --study)
#   Optional: --counts FILE  raw count matrix (Hugo_Symbol + sample columns) -> limma-voom
#             --cna-threshold 2

if (!exists("de_limma", mode = "function")) {
  here <- local({
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(f)) dirname(normalizePath(f)) else getwd()
  })
  source(file.path(here, "R", "analysis.R"))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  a <- as.list(setNames(args[seq(2, length(args), 2)], sub("^--", "", args[seq(1, length(args), 2)])))
  if (is.null(a$study) || is.null(a$out)) stop("--study and --out are required")
  thr <- if (is.null(a$`cna-threshold`)) 2 else as.numeric(a$`cna-threshold`)
  dir.create(a$out, showWarnings = FALSE, recursive = TRUE)

  use_counts <- !is.null(a$counts)
  expr <- to_patient_columns(read_cbio_matrix(if (use_counts) a$counts else file.path(a$study, "data_mrna_seq_v2_rsem.txt")))
  cna <- to_patient_columns(read_cbio_matrix(file.path(a$study, "data_cna.txt")))
  os <- parse_os(read_cbio_clinical(file.path(a$study, "data_clinical_patient.txt")))
  ids <- sort(Reduce(intersect, list(colnames(expr), colnames(cna), os$patient_id)))
  expr <- expr[, ids]; cna <- cna[, ids]; os <- os[match(ids, os$patient_id), ]
  status <- erbb2_status(cna, thr)
  message(sprintf("%d patients: %d ERBB2-amplified (GISTIC >= %g)", length(ids), sum(status == "Amplified"), thr))

  de <- de_limma(expr, status, counts = use_counts)
  utils::write.csv(de, file.path(a$out, "DE_ERBB2_amplified_vs_not.csv"), row.names = FALSE)
  amp <- de[de$gene %in% AMPLICON_17Q12, ]
  utils::write.csv(amp, file.path(a$out, "positive_control_17q12_amplicon.csv"), row.names = FALSE)
  gs <- gsea_go(de)
  if (!is.null(gs)) utils::write.csv(gs, file.path(a$out, "GSEA_GO_BP.csv"), row.names = FALSE)

  grDevices::pdf(file.path(a$out, "volcano.pdf"), 7, 6)
  with(de, plot(log2FC, -log10(adj_p_value), pch = 20, cex = 0.4, col = ifelse(gene %in% AMPLICON_17Q12, "firebrick", "grey55"),
                main = "ERBB2-amplified vs not (limma); red = 17q12 amplicon"))
  grDevices::dev.off()

  ok <- !is.na(os$time) & os$time > 0 & !is.na(os$event)
  x <- t(if (use_counts) log2(edgeR::cpm(expr) + 1) else log2(expr + 1))[ok, ]
  sig <- cv_lasso_cox(x, os$time[ok], os$event[ok])
  grp <- ifelse(sig$linear_predictor > stats::median(sig$linear_predictor), "High", "Low")
  lr <- survival::survdiff(survival::Surv(os$time[ok], os$event[ok]) ~ grp)
  lr_p <- stats::pchisq(lr$chisq, 1, lower.tail = FALSE)
  utils::write.csv(data.frame(patient_id = ids[ok], cv_risk = sig$linear_predictor, fold = sig$fold, group = grp),
                   file.path(a$out, "cv_risk_scores.csv"), row.names = FALSE)
  utils::write.csv(sig$selection_frequency, file.path(a$out, "lasso_selection_frequency.csv"), row.names = FALSE)
  grDevices::pdf(file.path(a$out, "km_cross_validated_risk.pdf"), 7, 6)
  plot(survival::survfit(survival::Surv(os$time[ok], os$event[ok]) ~ grp), col = c("firebrick", "steelblue"),
       xlab = "Months", ylab = "Overall survival",
       main = sprintf("Out-of-fold LASSO-Cox risk: C = %.3f, log-rank p = %.3g", sig$c_index, lr_p))
  graphics::legend("topright", c("High", "Low"), col = c("firebrick", "steelblue"), lty = 1)
  grDevices::dev.off()

  s <- list(n_patients = length(ids), n_amplified = sum(status == "Amplified"), cna_threshold = thr,
            de_method = if (use_counts) "limma-voom (counts)" else "limma-trend on log2(RSEM+1)",
            n_de_fdr05_lfc1 = sum(de$adj_p_value < 0.05 & abs(de$log2FC) >= 1),
            amplicon_genes_up = sum(amp$adj_p_value < 0.05 & amp$log2FC > 0), amplicon_genes_tested = nrow(amp),
            survival_n = sum(ok), survival_events = sum(os$event[ok]), cv_c_index = sig$c_index,
            cv_c_index_se = sig$c_index_se, cv_logrank_p = lr_p)
  jsonlite::write_json(s, file.path(a$out, "summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 4)
  invisible(s)
}

if (sys.nframe() == 0) main()
