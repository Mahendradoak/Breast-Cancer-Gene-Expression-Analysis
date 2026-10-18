root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
source(file.path(root, "R", "analysis.R"))
source(file.path(root, "run_analysis.R"))

make_study <- function(dir, n = 160, g = 400, seed = 1) {
  set.seed(seed)
  pats <- sprintf("TCGA-A%d-%04d", 1:n %% 9, 1:n)
  samples <- paste0(pats, "-01")
  genes <- c(AMPLICON_17Q12, sprintf("GENE%03d", seq_len(g - length(AMPLICON_17Q12))))
  amp <- rbinom(n, 1, 0.2) == 1
  cna <- matrix(sample(c(-1, 0, 1), g * n, TRUE, c(.2, .6, .2)), g, n, dimnames = list(genes, samples))
  cna["ERBB2", ] <- ifelse(amp, 2, sample(c(0, 1), n, TRUE))         # gains (1) are NOT amplification
  rsem <- matrix(2^rnorm(g * n, 7, 1), g, n, dimnames = dimnames(cna))
  rsem[AMPLICON_17Q12[1:5], amp] <- rsem[AMPLICON_17Q12[1:5], amp] * 8
  risk <- log2(rsem["GENE001", ] + 1) - 7
  time <- rexp(n, 0.02 * exp(1.2 * risk)); cens <- runif(n, 10, 150)
  # an extra metastatic sample that must be ignored
  rsem <- cbind(rsem, `TCGA-A1-0001-06` = rsem[, 1]); cna <- cbind(cna, `TCGA-A1-0001-06` = cna[, 1])
  wm <- function(m, f) data.table::fwrite(data.frame(Hugo_Symbol = rownames(m), Entrez_Gene_Id = seq_len(nrow(m)), m,
                                                     check.names = FALSE), f, sep = "\t")
  wm(rsem, file.path(dir, "data_mrna_seq_v2_rsem.txt")); wm(cna, file.path(dir, "data_cna.txt"))
  clin <- data.frame(PATIENT_ID = pats, OS_STATUS = ifelse(time < cens, "1:DECEASED", "0:LIVING"),
                     OS_MONTHS = round(pmin(time, cens), 2))
  writeLines(c("#Patient Identifier\tOverall Survival Status\tOverall Survival (Months)", "#x\tx\tx", "#STRING\tSTRING\tNUMBER", "#1\t1\t1"),
             file.path(dir, "data_clinical_patient.txt"))
  suppressWarnings(utils::write.table(clin, file.path(dir, "data_clinical_patient.txt"), sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE))
  list(amp = setNames(amp, pats))
}

test_that("GISTIC gain (1) is not amplification", {
  cna <- matrix(c(-1, 0, 1, 2), 1, dimnames = list("ERBB2", paste0("s", 1:4)))
  expect_equal(as.character(erbb2_status(cna)), c("Not_Amplified", "Not_Amplified", "Not_Amplified", "Amplified"))
})

test_that("only primary tumour columns are kept", {
  m <- matrix(1:3, 1, dimnames = list("g", c("TCGA-AA-0001-01", "TCGA-AA-0001-06", "TCGA-AA-0002-11")))
  expect_equal(colnames(to_patient_columns(m)), "TCGA-AA-0001")
})

test_that("cBioPortal clinical metadata lines are skipped", {
  d <- withr::local_tempdir(); make_study(d)
  cl <- read_cbio_clinical(file.path(d, "data_clinical_patient.txt"))
  expect_equal(names(cl), c("PATIENT_ID", "OS_STATUS", "OS_MONTHS"))
  expect_true(all(grepl("^TCGA-", cl$PATIENT_ID)))
})

test_that("full run: amplicon recovered, CV signature uses held-out predictions", {
  d <- withr::local_tempdir(); sim <- make_study(d)
  out <- file.path(d, "out")
  s <- suppressMessages(main(c("--study", d, "--out", out)))
  expect_equal(s$n_patients, 160)
  expect_equal(s$n_amplified, sum(sim$amp))
  expect_gte(s$amplicon_genes_up, 5)
  expect_gt(s$cv_c_index, 0.6)
  r <- read.csv(file.path(out, "cv_risk_scores.csv"))
  expect_equal(sort(unique(r$fold)), 1:5)
  freq <- read.csv(file.path(out, "lasso_selection_frequency.csv"))
  expect_equal(freq$gene[1], "GENE001")
})
