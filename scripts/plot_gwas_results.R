#!/usr/bin/env Rscript
# =============================================================================
# GWAS Results Visualization
# Generates Manhattan plot, QQ plot, and PCA plot
#
# Usage: Rscript plot_gwas_results.R <gwas_file> <pca_eigenvec> <pca_eigenval> \
#            <manhattan_out> <qq_out> <pca_out>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript plot_gwas_results.R <gwas_file> <pca_eigenvec> <pca_eigenval> <manhattan_out> <qq_out> <pca_out>")
}

gwas_file     <- args[1]
pca_vec_file  <- args[2]
pca_val_file  <- args[3]
manhattan_out <- args[4]
qq_out        <- args[5]
pca_out       <- args[6]

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Ensure output directory exists
dir.create(dirname(manhattan_out), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Read GWAS results
# -----------------------------------------------------------------------------
gwas <- read.table(gwas_file, header = TRUE, sep = "\t",
                   comment.char = "", stringsAsFactors = FALSE)

gwas_clean <- gwas %>%
  filter(!is.na(P) & P > 0 & P <= 1) %>%
  select(CHR = X.CHROM, BP = POS, SNP = ID, P)

cat("GWAS results loaded:", nrow(gwas_clean), "variants\n")

# -----------------------------------------------------------------------------
# Manhattan Plot (ggplot2 — clear chromosome separation)
# -----------------------------------------------------------------------------
cat("Generating Manhattan plot...\n")

# Define chromosome order: autosomes 1–27, then Z, W, then scaffolds grouped
autosomes <- as.character(1:27)
sex_chr <- c("Z", "W")
scaffolds <- sort(unique(gwas_clean$CHR[grepl("^NW_", gwas_clean$CHR)]))

chr_order <- c(autosomes, sex_chr, scaffolds)

# Keep only chromosomes present in the data
chr_order <- chr_order[chr_order %in% unique(gwas_clean$CHR)]

gwas_clean <- gwas_clean %>%
  filter(CHR %in% chr_order) %>%
  mutate(CHR = factor(CHR, levels = chr_order))

# Build cumulative genomic positions for x axis
chr_lengths <- gwas_clean %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP), .groups = "drop") %>%
  mutate(cum_start = cumsum(as.numeric(lag(chr_len, default = 0))),
         gap = (row_number() - 1) * 2e6)   # small gap between chromosomes

# Add cumulative position to every SNP
gwas_plot <- gwas_clean %>%
  inner_join(chr_lengths %>% select(CHR, cum_start, gap), by = "CHR") %>%
  mutate(BP_cum = BP + cum_start + gap,
         logp = -log10(P))

# Midpoint of each chromosome for axis labels
axis_df <- gwas_plot %>%
  group_by(CHR) %>%
  summarise(center = (min(BP_cum) + max(BP_cum)) / 2, .groups = "drop")

# Label: show "1"–"27", "Z", "W"; collapse scaffolds into one label
axis_df <- axis_df %>%
  mutate(label = ifelse(grepl("^NW_", CHR), "", as.character(CHR)))

# Merge scaffolds into one label at their midpoint
scaff_rows <- which(axis_df$label == "")
if (length(scaff_rows) > 0) {
  mid_scaff <- scaff_rows[ceiling(length(scaff_rows) / 2)]
  axis_df$label[mid_scaff] <- "Scaff"
}

# Alternating colors per chromosome
gwas_plot$chr_num <- as.numeric(gwas_plot$CHR)

# Two-tone palette
col_even <- "#2166ac"
col_odd  <- "#67a9cf"

p_manh <- ggplot(gwas_plot, aes(x = BP_cum, y = logp,
                                colour = factor(chr_num %% 2))) +
  geom_point(alpha = 0.6, size = 0.8, shape = 16) +
  scale_colour_manual(values = c("0" = col_even, "1" = col_odd), guide = "none") +
  scale_x_continuous(breaks = axis_df$center,
                     labels = axis_df$label,
                     expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  geom_hline(yintercept = -log10(5e-8), colour = "red",
             linetype = "solid", linewidth = 0.4) +
  geom_hline(yintercept = -log10(1e-5), colour = "blue",
             linetype = "dashed", linewidth = 0.3) +
  labs(x = "Chromosome", y = expression(-log[10](p)),
       title = "Manhattan Plot — GWAS Results") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.text.x  = element_text(size = 8, angle = 0, vjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 8)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border = element_blank(),
    plot.margin  = margin(10, 15, 10, 10)
  )

ggsave(manhattan_out, p_manh, width = 14, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# QQ Plot
# -----------------------------------------------------------------------------
cat("Generating QQ plot...\n")

n <- nrow(gwas_clean)
qq_df <- data.frame(
  observed = sort(-log10(gwas_clean$P)),
  expected = sort(-log10(ppoints(n)))
)

# Genomic inflation factor
chisq  <- qchisq(1 - gwas_clean$P, 1)
lambda <- median(chisq) / qchisq(0.5, 1)

p_qq <- ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.4) +
  geom_point(alpha = 0.5, size = 0.8, colour = "#3498db") +
  annotate("text", x = 0.5, y = max(qq_df$observed) * 0.95,
           label = bquote(lambda == .(round(lambda, 3))),
           hjust = 0, size = 4.5) +
  labs(x = expression(Expected ~ -log[10](p)),
       y = expression(Observed ~ -log[10](p)),
       title = "QQ Plot — GWAS Results") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  coord_fixed()

ggsave(qq_out, p_qq, width = 6, height = 6, dpi = 300)

cat("Genomic inflation factor (lambda):", round(lambda, 3), "\n")

# -----------------------------------------------------------------------------
# PCA Plot
# -----------------------------------------------------------------------------
cat("Generating PCA plot...\n")

eigenvec <- read.table(pca_vec_file, header = TRUE, comment.char = "")
eigenval <- read.table(pca_val_file, header = FALSE)$V1
var_explained <- eigenval / sum(eigenval) * 100

pca_plot <- ggplot(eigenvec, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.7, color = "#3498db") +
  labs(
    title = "PCA — Population Structure",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(pca_out, pca_plot, width = 8, height = 6)

# -----------------------------------------------------------------------------
# Summary Statistics
# -----------------------------------------------------------------------------
cat("\n=== GWAS Summary ===\n")
cat("Total variants tested:", nrow(gwas_clean), "\n")
cat("Suggestive hits (p < 1e-5):", sum(gwas_clean$P < 1e-5), "\n")
cat("Genome-wide significant (p < 5e-8):", sum(gwas_clean$P < 5e-8), "\n")

cat("\nTop 10 associations:\n")
top_hits <- gwas_clean %>%
  arrange(P) %>%
  head(10)
print(top_hits)

cat("\nPlots saved successfully!\n")
