#!/usr/bin/env Rscript
# Extract THUNDER W and H matrices via single-run NMF (for comparison plotting).
# The full THUNDER pipeline (run_thunder_fullgenome.R) uses nrun=10 for H estimation
# but does not save W. This script runs nrun=1 to get W for downstream analysis.
# Input: output/simulation_100000/thunder_input.tsv
# Output: output/simulation_100000/thunder_W_k4.tsv, thunder_H_k4_raw.tsv

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_thunder_extract_W.R <input.tsv> <output_W.tsv> <output_H.tsv>")
}
input_path <- args[1]
output_W_path <- args[2]
output_H_path <- args[3]

user_lib <- "/home/koki/R/x86_64-conda-linux-gnu-library/4.5"
.libPaths(c(user_lib, .libPaths()))
library(NMF)

cat("Loading THUNDER input...\n")

raw_mix <- read.delim(input_path, check.names=FALSE)
mix_mat <- as.matrix(raw_mix[, 3:ncol(raw_mix)])
rs <- rowSums(mix_mat)
keep <- rs > 0

cat("Running NMF for W extraction (K=4, nrun=1)...\n")
fit <- NMF::nmf(mix_mat[keep,], rank=4, method="brunet", seed=42, nrun=1)

W_thunder <- basis(fit)
H_thunder <- coef(fit)

cat("W:", nrow(W_thunder), "x", ncol(W_thunder), "\n")
cat("H:", nrow(H_thunder), "x", ncol(H_thunder), "\n")

# Save W with feature info
feature_names <- raw_mix$feature_name[keep]
out <- data.frame(feature_name=feature_names, W_thunder)
write.table(out, output_W_path, sep="\t", quote=FALSE, row.names=FALSE)

# Also save H
write.table(H_thunder, output_H_path, sep="\t", quote=FALSE)

cat("[saved]", output_W_path, output_H_path, "\n")
