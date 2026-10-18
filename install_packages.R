# Tested with R 4.3.3, limma 3.58.1, edgeR 4.0.16, glmnet 4.1-8, survival 3.5-8.
install.packages(c("data.table", "jsonlite", "survival", "glmnet", "testthat", "withr", "BiocManager"),
                 repos = "https://cloud.r-project.org")
BiocManager::install(c("limma", "edgeR", "clusterProfiler", "org.Hs.eg.db"), update = FALSE, ask = FALSE)
