#want to arcsin transform the allelic percent data for HK and HU KO in order to run a wilcoxon test on each gene

library(tidyverse)
library(readxl)

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

#  File with B6, CAST and % allelic numbers
df <- read_excel("/Users/murvinmm/Downloads/HK_HU_allelic_table_genes_XI_A_K.xlsx")

# second table with chr, start, stop columns
coords <- read.delim("/Users/murvinmm/Downloads/2025_04_17_TSC_HU_HK_KD_no_dox_RNAseq_SNP_reads_all_data.txt")  

# Join on gene name 
df <- df %>%
  left_join(coords %>% select(gene, chr, start, end), by = "gene") %>%
  relocate(chr, start, end, .after = gene)

#label Allelic % columns for tests (This is b6/b6+cast as a decimal (0-1)) 
HK_nodox_cols <- c("HK_nodox_F8_allelic_pct", "HK_nodox_F14_allelic_pct",
                   "HK_nodox_F18_r1_allelic_pct", "HK_nodox_F18_r2_allelic_pct")

HK_dox_cols   <- c("HK_dox_F8_allelic_pct", "HK_dox_F14_allelic_pct",
                   "HK_dox_F18_r1_allelic_pct", "HK_dox_F18_r2_allelic_pct")

HU_nodox_cols <- c("HU_nodox_F8_allelic_pct", "HU_nodox_F14_allelic_pct",
                   "HU_nodox_F18_r1_allelic_pct", "HU_nodox_F18_r2_allelic_pct")

HU_dox_cols   <- c("HU_dox_F8_allelic_pct", "HU_dox_F14_allelic_pct",
                   "HU_dox_F18_r1_allelic_pct", "HU_dox_F18_r2_allelic_pct")

# Arcsine transform 
arcsine_transform <- function(p) {
  p <- pmax(0, pmin(1, p))  # clamp to [0,1]
  asin(sqrt(p))
}

# Wilcoxon test
run_wilcoxon <- function(nodox_vals, dox_vals) {
  nodox_t <- arcsine_transform(as.numeric(nodox_vals))
  dox_t   <- arcsine_transform(as.numeric(dox_vals))
  result  <- tryCatch(
    wilcox.test(nodox_t, dox_t, alternative = "less", exact = FALSE),
    error = function(e) list(p.value = NA)
  )
  result$p.value
}

# Apply per gene
df <- df %>%
  rowwise() %>%
  mutate(
    pval_HK = run_wilcoxon(c_across(all_of(HK_nodox_cols)),
                           c_across(all_of(HK_dox_cols))),
    pval_HU = run_wilcoxon(c_across(all_of(HU_nodox_cols)),
                           c_across(all_of(HU_dox_cols)))
  ) %>%
  ungroup()

# Export with p vals
write.csv(df, "/Users/murvinmm/Downloads/HK_HU_allelic_table_genes_XI_A_K_p_vals_one_side_coords.csv", row.names = FALSE)

# Filter: p < 0.06 in HK OR HU, plus always include Xist and Airn
anchor_genes <- c("Xist", "Airn")

sig <- df %>%
  filter(pval_HK < 0.06 | pval_HU < 0.06 | gene %in% anchor_genes)

# Heatmap function
make_heatmap <- function(dat, title, row_split_col = NULL) {  # FIX: added row_split_col parameter
  
  # Order by start position
  dat <- dat %>% arrange(start)
  
  gene_names <- dat$gene
  
  # Label coloring: pink = HK only, orange = HU only, blend = both
  hk_sig  <- !is.na(dat$pval_HK) & dat$pval_HK < 0.06
  hu_sig  <- !is.na(dat$pval_HU) & dat$pval_HU < 0.06
  anchor  <- gene_names %in% anchor_genes
  

  
  gene_colors <- case_when(
    anchor            ~ "black",
    hk_sig & hu_sig   ~ "black",
    hk_sig            ~ "turquoise3",
    hu_sig            ~ "#FFA500",
    TRUE              ~ "grey50"
  )
  
  # Build matrix: HK nodox | HK dox | HU nodox | HU dox
  mat_HK_nodox <- as.matrix(dat[, HK_nodox_cols])
  mat_HK_dox   <- as.matrix(dat[, HK_dox_cols])
  mat_HU_nodox <- as.matrix(dat[, HU_nodox_cols])
  mat_HU_dox   <- as.matrix(dat[, HU_dox_cols])
  
  mat <- cbind(mat_HK_nodox, mat_HK_dox, mat_HU_nodox, mat_HU_dox)
  rownames(mat) <- gene_names
  
  # Row-wise min/max scaling
  mat_scaled <- t(apply(mat, 1, function(row) {
    rng <- range(row, na.rm = TRUE)
    if (diff(rng) == 0) return(row - rng[1])
    (row - rng[1]) / diff(rng)
  }))
  
  # Color scale: blue (low) → white (mid) → red (high)
  col_fun <- colorRamp2(c(0, 0.5, 1), c("navy", "white", "salmon"))
  
  # Column split to enforce HK | gap | HU grouping
  col_split <- factor(
    c(rep("HK nodox", 4), rep("HK dox", 4), rep("HU nodox", 4), rep("HU dox", 4)),
    levels = c("HK nodox", "HK dox", "HU nodox", "HU dox")
  )
  
  # row_split_col  as a function parameter
  row_split_vec <- if (!is.null(row_split_col)) {
    factor(dat[[row_split_col]], levels = unique(dat[[row_split_col]]))
  } else NULL
  
  cell_size <- unit(6, "mm")  # adjust to resize everything
  
  Heatmap(
    mat_scaled,
    name            = "Row-scaled\nallelic %",
    col             = col_fun,
    cluster_rows    = FALSE,
    cluster_columns = FALSE,
    show_row_names  = TRUE,
    row_names_gp    = gpar(col = gene_colors, fontsize = convertUnit(cell_size, "pt"), fontfamily = "Arial"),
    column_split    = col_split,
    column_title    = title,
    column_gap      = unit(c(2, 6, 2), "mm"),  # 6mm gap between HK dox and HU nodox
    row_split       = row_split_vec,
    row_gap         = unit(4, "mm"),
    row_title_gp    = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
    na_col          = "grey80",
    border          = TRUE,
    row_names_side  = "left",
    width           = ncol(mat_scaled) * cell_size,
    height          = nrow(mat_scaled) * cell_size
  )
}

# Build combined dataset ordered by chr section then start position
sig_combined <- sig %>%
  mutate(section = case_when(
    chr == "chrX"  ~ "X chromosome",
    chr == "chr17" ~ "Airn domain",
    TRUE           ~ "Other"
  )) %>%
  filter(section != "Other") %>%
  arrange(factor(section, levels = c("X chromosome", "Airn domain")), start)

# Draw and save single combined heatmap
pdf("/Users/murvinmm/Downloads/HK_HU_heatmap_combined.pdf", width = 8, height = 18)
draw(make_heatmap(sig_combined, "HK and HU Allelic %", row_split_col = "section"))
dev.off()

# Export table with p-values
write.csv(df, "/Users/murvinmm/Downloads/HK_HU_allelic_table_genes_XI_A_K_p_vals.csv", row.names = FALSE)