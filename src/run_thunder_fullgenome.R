#!/usr/bin/env Rscript
# THUNDER full-genome pipeline: initial NMF -> feature selection -> second NMF
# Requires: NMF, thundeR (github: brycerowland/thundeR), dplyr
# Input: output/simulation_100000/thunder_input.tsv
# Output: output/simulation_100000/thunder_H_fullgenome.tsv
# Runtime: ~28 min on full genome (1.17M features, nrun=10)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_thunder_fullgenome.R <input.tsv> <output_H.tsv>")
}
input_path <- args[1]
output_path <- args[2]

user_lib <- "/home/koki/R/x86_64-conda-linux-gnu-library/4.5"
.libPaths(c(user_lib, .libPaths()))
library(NMF)
library(thundeR)
library(dplyr)

cat("=== THUNDER Full Genome Pipeline ===\n\n")

# Load full data
cat("[load] reading", input_path, "...\n")
raw_mix <- read.delim(input_path, check.names=FALSE)
cat("  ", nrow(raw_mix), "features x", ncol(raw_mix)-2, "samples\n")

mix_mat <- as.matrix(raw_mix[, 3:ncol(raw_mix)])
rs <- rowSums(mix_mat)
keep <- rs > 0
cat("  nonzero rows:", sum(keep), "\n")

# Step 1: Initial NMF
cat("\n[Step 1] Initial NMF (K=4, nrun=10, sequential)...\n")
t0 <- proc.time()
fit_init <- NMF::nmf(mix_mat[keep,], rank=4, method="brunet",
                      seed="random", nrun=10, .pbackend="seq")
elapsed <- (proc.time() - t0)[3]
cat("  Done in", round(elapsed/60, 1), "min\n")

# Step 2: Feature selection
cat("\n[Step 2] THUNDER feature selection...\n")
subset_mix <- thundeR::subset_init_nmf(
  mix_mat[keep,], fit_init,
  feature_name = raw_mix$feature_name[keep],
  contact_type = raw_mix$contact_type[keep]
)
cat("  Selected features:", nrow(subset_mix), "/", sum(keep), "\n")

if (nrow(subset_mix) == 0) {
  cat("  WARNING: No features selected! Using initial NMF result.\n")
  H_final <- coef(fit_init)
} else {
  # Step 3: Second NMF
  cat("\n[Step 3] Second NMF on selected features (K=4, nrun=10)...\n")
  subset_mat <- as.matrix(subset_mix %>% select(-feature_name, -contact_type))

  t0 <- proc.time()
  fit_final <- NMF::nmf(subset_mat, rank=4, method="brunet",
                         seed="random", nrun=10, .pbackend="seq")
  elapsed <- (proc.time() - t0)[3]
  cat("  Done in", round(elapsed/60, 1), "min\n")

  H_final <- coef(fit_final)
}

# Proportions
H_norm <- sweep(H_final, 2, colSums(H_final), "/")
cat("\n=== THUNDER Estimated Proportions (Full Genome) ===\n")
print(round(t(H_norm), 3))

write.table(t(H_norm), output_path, sep="\t", quote=FALSE)
cat("\n[saved]", output_path, "\n")
