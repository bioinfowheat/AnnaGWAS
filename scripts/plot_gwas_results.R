#!/usr/bin/env Rscript
# =============================================================================
# GWAS Results Visualization
# Generates Manhattan plot, QQ plot, and PCA plot
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(qqman)
library(RColorBrewer)

# -----------------------------------------------------------------------------
# Read GWAS results
# -----------------------------------------------------------------------------
gwas <- read.table(snakemake@input[["gwas"]], header = TRUE, sep = "\t",
                   comment.char = "", stringsAsFactors = FALSE)

# Rename columns for qqman compatibility
gwas_clean <- gwas %>%
  filter(!is.na(P)) %>%
  select(CHR = `#CHROM`, BP = POS, SNP = ID, P) %>%
  mutate(
    # Convert chromosome to numeric (handle non-numeric chromosomes)
    CHR = as.numeric(factor(CHR, levels = unique(CHR)))
  ) %>%
  filter(!is.na(P) & P > 0 & P <= 1)

cat("GWAS results loaded:", nrow(gwas_clean), "variants\n")

# -----------------------------------------------------------------------------
# Manhattan Plot
# -----------------------------------------------------------------------------
cat("Generating Manhattan plot...\n")

pdf(snakemake@output[["manhattan"]], width = 12, height = 6)

manhattan(gwas_clean,
          chr = "CHR",
          bp = "BP",
          snp = "SNP",
          p = "P",
          main = "Manhattan Plot - GWAS Results",
          col = c("#3498db", "#2c3e50"),
          suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-8),
          cex = 0.6,
          cex.axis = 0.8)

dev.off()

# -----------------------------------------------------------------------------
# QQ Plot
# -----------------------------------------------------------------------------
cat("Generating QQ plot...\n")

pdf(snakemake@output[["qq"]], width = 6, height = 6)

qq(gwas_clean$P,
   main = "QQ Plot - GWAS Results",
   cex = 0.6,
   col = "#3498db")

# Calculate genomic inflation factor
chisq <- qchisq(1 - gwas_clean$P, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
legend("topleft",
       legend = bquote(lambda == .(round(lambda, 3))),
       bty = "n",
       cex = 1.2)

dev.off()

cat("Genomic inflation factor (lambda):", round(lambda, 3), "\n")

# -----------------------------------------------------------------------------
# PCA Plot
# -----------------------------------------------------------------------------
cat("Generating PCA plot...\n")

# Read PCA results
eigenvec <- read.table(snakemake@input[["pca_vec"]], header = TRUE)
eigenval <- read.table(snakemake@input[["pca_val"]], header = FALSE)$V1

# Calculate variance explained
var_explained <- eigenval / sum(eigenval) * 100

# Create PCA plot
pca_plot <- ggplot(eigenvec, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.7, color = "#3498db") +
  labs(
    title = "PCA - Population Structure",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(snakemake@output[["pca"]], pca_plot, width = 8, height = 6)

# -----------------------------------------------------------------------------
# Summary Statistics
# -----------------------------------------------------------------------------
cat("\n=== GWAS Summary ===\n")
cat("Total variants tested:", nrow(gwas_clean), "\n")
cat("Suggestive hits (p < 1e-5):", sum(gwas_clean$P < 1e-5), "\n")
cat("Genome-wide significant (p < 5e-8):", sum(gwas_clean$P < 5e-8), "\n")

# Top hits
cat("\nTop 10 associations:\n")
top_hits <- gwas_clean %>%
  arrange(P) %>%
  head(10)
print(top_hits)

cat("\nPlots saved successfully!\n")
