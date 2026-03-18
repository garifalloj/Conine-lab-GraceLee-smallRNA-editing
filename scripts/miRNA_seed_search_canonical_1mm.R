#!/usr/bin/env Rscript

# ===============================================================
# ONE-SHOT miRNA seed site analysis (canonical + 7mer-1mm) in 3'UTRs
#
# What this script does (single run):
#  1) Loads DESeq2 results (gene_id, log2FoldChange, pvalue; optional baseMean)
#  2) Fetches Ensembl 3'UTRs (version-pinned) and cleans sequences
#  3) Counts canonical sites per transcript (8mer, 7mer-m8, 7mer-A1, 6mer)
#  4) Counts 7mer-1mm (any) = EXACTLY 1 mismatch within seed(2–8) complement
#  5) Collapses to gene-level counts + flags; writes gene-level CSV
#  6) Fisher enrichment (UP vs BG_rest; DOWN vs BG_rest) for:
#        6mer, 7mer-m8, 7mer-A1, 8mer, 7mer-1mm (any)
#     -> writes fisher_results_combined.csv
#  7) Makes ONE grouped barplot including those 5 classes,
#     with per-bar gene counts and p-values annotated.
#  8) ECDF plots ONLY using "No site" baselines:
#        A) Canonical exclusive class ECDF vs baseline = No canonical site
#        B) Mismatch-category ECDF vs baseline = no canonical & no 1mm
#
# Outputs:
#   - gene_level_site_counts_and_flags_canon_plus_1mm.csv
#   - fisher_results_combined.csv
#   - barplot_site_classes_canon_plus_1mm_annotated.png
#   - ecdf_canonical_exclusive_vs_NoCanonical.png
#   - ecdf_mismatch_categories_vs_NoSite.png
#   - ecdf_stats_canonical_exclusive_vs_NoCanonical.csv
#   - ecdf_stats_mismatch_categories_vs_NoSite.csv
#   - session_info.txt
# ===============================================================

setwd("/2025 Dec Analysis/miR465_parthenote_morula")

# -------------------- PACKAGES --------------------
suppressPackageStartupMessages({
  req <- c("biomaRt","Biostrings","dplyr","ggplot2","readr","tidyr","purrr","scales","stringr","tibble")
  to_install <- req[!sapply(req, requireNamespace, quietly=TRUE)]
  if (length(to_install)) install.packages(to_install, repos="https://cloud.r-project.org")
  if (!suppressWarnings(requireNamespace("BiocManager", quietly=TRUE))) install.packages("BiocManager")
  if (!suppressWarnings(requireNamespace("Biostrings", quietly=TRUE))) BiocManager::install("Biostrings", ask=FALSE, update=FALSE)
  
  library(biomaRt); library(Biostrings); library(dplyr); library(ggplot2)
  library(readr); library(tidyr); library(purrr); library(scales); library(stringr); library(tibble)
})

# -------------------- USER INPUTS --------------------
deseq_file <- NULL  # set path OR leave NULL to choose interactively
outdir     <- "miRNA_site_outputs_CANON_plus_1mm_low"

mirna_seq  <- "GAUCAGGGCCUUUCUAAGUAGA"  # miR-465c-3p (RNA 5'->3')
mirna_name <- "miR-465c-3p"
#mirna_seq  <- "UAAUACUGCCGGGUAAUGAUGGA"
#mirna_name <- "miR-200c-3p"

ensembl_ver <- 110
species_ds  <- "mmusculus_gene_ensembl"

alpha <- 0.05  # DE p-value threshold for UP/DOWN definitions

# -------------------- I/O SETUP --------------------
if (is.null(deseq_file)) {
  message("Choose your DESeq2 CSV (must include columns: gene_id, log2FoldChange, pvalue; optional baseMean).")
  deseq_file <- file.choose()
}
if (!file.exists(deseq_file)) stop("DESeq2 results file not found: ", deseq_file)

OUTDIR <- normalizePath(outdir, mustWork = FALSE)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

message("==> DESeq2 file: ", deseq_file)
message("==> Output dir:  ", OUTDIR)
message("==> miRNA:       ", mirna_name, " (", mirna_seq, ")")

# -------------------- LOAD DESeq2 RESULTS --------------------
deseq_data <- suppressMessages(readr::read_csv(deseq_file, show_col_types = FALSE))

needed <- c("gene_id","log2FoldChange","pvalue")
if (!all(needed %in% names(deseq_data))) {
  stop("Input DESeq2 CSV must contain: gene_id, log2FoldChange, pvalue. Found: ",
       paste(names(deseq_data), collapse=", "))
}

deseq_data <- deseq_data %>%
  mutate(
    gene_id        = as.character(gene_id),
    log2FoldChange = as.numeric(log2FoldChange),
    pvalue         = as.numeric(pvalue)
  ) %>%
  filter(!is.na(gene_id), is.finite(log2FoldChange))

has_baseMean <- "baseMean" %in% names(deseq_data)
if (has_baseMean) deseq_data <- deseq_data %>% mutate(baseMean = as.numeric(baseMean))

sig_data  <- deseq_data %>% filter(!is.na(pvalue), pvalue < alpha)
up_genes   <- unique(sig_data$gene_id[sig_data$log2FoldChange > 0])
down_genes <- unique(sig_data$gene_id[sig_data$log2FoldChange < 0])
all_genes  <- unique(deseq_data$gene_id)

message(sprintf("Genes in DE table: all=%d, up(sig p<%.2f)=%d, down(sig)=%d",
                length(all_genes), alpha, length(up_genes), length(down_genes)))

# -------------------- CONNECT TO ENSEMBL (version pinned) --------------------
message("Connecting to Ensembl...")
ensembl <- tryCatch({
  if (is.na(ensembl_ver)) useEnsembl(biomart="ensembl") else useEnsembl(biomart="ensembl", version=ensembl_ver)
}, error=function(e) {
  message("Falling back to default Ensembl mirror: ", e$message)
  useMart("ensembl")
})
mart_ds <- useDataset(species_ds, mart = ensembl)

safe_getBM <- function(...) {
  for (i in 1:3) {
    out <- tryCatch(getBM(...), error=function(e) NULL)
    if (!is.null(out)) return(out)
    Sys.sleep(1 + i)
  }
  stop("biomaRt getBM failed repeatedly.")
}

is_ensembl_id <- function(x) {
  if (!length(x)) return(FALSE)
  mean(grepl("^ENSMUSG\\d+", x)) > 0.6
}

input_are_ensembl <- is_ensembl_id(all_genes)
id_filter <- if (input_are_ensembl) "ensembl_gene_id" else "external_gene_name"
message("ID auto-detect: treating gene_id as ", if (input_are_ensembl) "ENSEMBL IDs" else "GENE SYMBOLS")

# Map input IDs <-> Ensembl IDs
message("Mapping input IDs <-> Ensembl...")
gene_map <- safe_getBM(
  attributes = c("ensembl_gene_id","external_gene_name"),
  filters    = id_filter,
  values     = all_genes,
  mart       = mart_ds
) %>% distinct()

if (!nrow(gene_map)) stop("No mappings returned from Ensembl for your gene_id values.")

gene_map <- gene_map %>%
  mutate(gene_id_input = if (input_are_ensembl) ensembl_gene_id else external_gene_name)

# Fetch UTRs
message("Fetching 3'UTRs...")
utr <- safe_getBM(
  attributes = c("ensembl_gene_id","ensembl_transcript_id","external_gene_name","3utr"),
  filters    = "ensembl_gene_id",
  values     = unique(gene_map$ensembl_gene_id),
  mart       = mart_ds
)

utr <- utr %>% filter(!is.na(`3utr`), nchar(`3utr`) > 0)
if (!nrow(utr)) stop("No 3'UTR sequences returned from biomaRt.")

clean_utr <- function(x) {
  x <- gsub("\\s+", "", x)
  x <- toupper(x)
  x <- chartr("U","T", x)
  gsub("[^ACGTRYMKSWBDHVN]", "N", x)
}
utr <- utr %>% mutate(`3utr` = clean_utr(`3utr`), utr_len = nchar(`3utr`)) %>% filter(utr_len > 0)

message("UTRs retrieved for ", nrow(utr), " transcripts across ",
        length(unique(utr$ensembl_gene_id)), " genes.")

# -------------------- BUILD SEED PATTERNS --------------------
mirna_to_patterns <- function(mirna_rna_5to3) {
  mir <- toupper(gsub("U","T", mirna_rna_5to3))
  if (nchar(mir) < 8) stop("miRNA sequence too short; need >= 8 nt.")
  seed_2_8 <- substr(mir, 2, 8)  # miRNA positions 2–8
  seed_2_7 <- substr(mir, 2, 7)
  
  seed_2_8_rc <- as.character(reverseComplement(DNAString(seed_2_8)))
  seed_2_7_rc <- as.character(reverseComplement(DNAString(seed_2_7)))
  
  list(
    eightmer     = paste0(seed_2_8_rc, "A"),  # 8mer-A1
    seven_m8     = seed_2_8_rc,               # 7mer-m8
    seven_A1     = paste0(seed_2_7_rc, "A"),  # 7mer-A1
    sixmer_2to7  = seed_2_7_rc,               # 6mer
    seed7        = seed_2_8_rc                # 7mer for mismatch scan (exactly 1 mismatch)
  )
}
patterns <- mirna_to_patterns(mirna_seq)
message("Seed patterns (DNA, UTR-sense):")
print(patterns)

# -------------------- COUNT SITES PER TRANSCRIPT --------------------
count_sites_per_transcript <- function(utr_df, patterns) {
  s <- DNAStringSet(utr_df$`3utr`)
  
  # canonical exact matches
  vcount0 <- function(pat) vcountPattern(pat, s, fixed=TRUE)
  
  cnt_8mer    <- vcount0(patterns$eightmer)
  cnt_7mer_m8 <- vcount0(patterns$seven_m8)
  cnt_7mer_A1 <- vcount0(patterns$seven_A1)
  cnt_6mer    <- vcount0(patterns$sixmer_2to7)
  cnt_any_canon <- cnt_8mer + cnt_7mer_m8 + cnt_7mer_A1 + cnt_6mer
  
  # EXACTLY 1 mismatch vs the 7mer seed
  cnt_7mer_max1 <- vcountPattern(patterns$seed7, s,
                                 max.mismatch=1, with.indels=FALSE, fixed=TRUE)
  cnt_7mer_0mm  <- vcountPattern(patterns$seed7, s,
                                 max.mismatch=0, with.indels=FALSE, fixed=TRUE)
  cnt_7mer_1mm  <- pmax(cnt_7mer_max1 - cnt_7mer_0mm, 0L)
  
  utr_df %>%
    mutate(
      cnt_8mer      = cnt_8mer,
      cnt_7mer_m8   = cnt_7mer_m8,
      cnt_7mer_A1   = cnt_7mer_A1,
      cnt_6mer      = cnt_6mer,
      cnt_any_canon = cnt_any_canon,
      cnt_7mer_1mm  = cnt_7mer_1mm
    )
}

message("Counting canonical + 7mer-1mm per transcript...")
utr_cnt <- count_sites_per_transcript(utr, patterns)

# -------------------- COLLAPSE TO GENE LEVEL --------------------
gene_level <- utr_cnt %>%
  group_by(ensembl_gene_id, external_gene_name) %>%
  summarize(
    cnt_8mer      = sum(cnt_8mer,      na.rm=TRUE),
    cnt_7mer_m8   = sum(cnt_7mer_m8,   na.rm=TRUE),
    cnt_7mer_A1   = sum(cnt_7mer_A1,   na.rm=TRUE),
    cnt_6mer      = sum(cnt_6mer,      na.rm=TRUE),
    cnt_any_canon = sum(cnt_any_canon, na.rm=TRUE),
    cnt_7mer_1mm  = sum(cnt_7mer_1mm,  na.rm=TRUE),
    utr_len_med   = median(utr_len,    na.rm=TRUE),
    .groups="drop"
  ) %>%
  left_join(gene_map, by=c("ensembl_gene_id","external_gene_name")) %>%
  mutate(
    site_8mer      = cnt_8mer      > 0,
    site_7mer_m8   = cnt_7mer_m8   > 0,
    site_7mer_A1   = cnt_7mer_A1   > 0,
    site_6mer      = cnt_6mer      > 0,
    site_any_canon = cnt_any_canon > 0,
    site_7mer_1mm  = cnt_7mer_1mm  > 0,
    site_any_canon_or_1mm = site_any_canon | site_7mer_1mm,
    site_no_canon  = !site_any_canon,
    site_no_any    = !site_any_canon_or_1mm
  )

readr::write_csv(
  gene_level %>% arrange(external_gene_name),
  file.path(OUTDIR, "gene_level_site_counts_and_flags_canon_plus_1mm.csv")
)

message("Gene-level rows with UTRs: ", nrow(gene_level))

# -------------------- DEFINE UP/DOWN/BG USING INPUT KEY --------------------
up_gene <- gene_level %>% semi_join(tibble(gene_id = up_genes), by=c("gene_id_input"="gene_id"))
down_gene <- gene_level %>% semi_join(tibble(gene_id = down_genes), by=c("gene_id_input"="gene_id"))
bg_gene <- gene_level

message(sprintf("Up genes with UTRs:   %d", nrow(up_gene)))
message(sprintf("Down genes with UTRs: %d", nrow(down_gene)))
message(sprintf("BG genes with UTRs:   %d", nrow(bg_gene)))

# -------------------- FISHER + BARPLOT (ONE PLOT) --------------------
fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 1e-4 ~ formatC(p, format="e", digits=2),
    p < 1e-3 ~ sprintf("%.3g", p),
    TRUE     ~ sprintf("%.3f", p)
  )
}

# Fisher compares group vs BG_rest (BG excluding group genes)
fisher_for_flag <- function(df_group, df_bg, flag, contrast_label) {
  if (!nrow(df_group) || !nrow(df_bg)) {
    return(tibble(flag=flag, contrast=contrast_label, p=NA_real_, OR=NA_real_, lo=NA_real_, hi=NA_real_))
  }
  bg_rest <- df_bg %>% anti_join(df_group %>% distinct(external_gene_name), by="external_gene_name")
  a <- sum(df_group[[flag]], na.rm=TRUE)
  b <- nrow(df_group) - a
  c <- sum(bg_rest[[flag]], na.rm=TRUE)
  d <- nrow(bg_rest) - c
  ft <- fisher.test(matrix(c(a,b,c,d), nrow=2, byrow=TRUE))
  tibble(flag=flag, contrast=contrast_label, p=unname(ft$p.value),
         OR=unname(ft$estimate), lo=unname(ft$conf.int[1]), hi=unname(ft$conf.int[2]))
}

# Site classes requested (exactly these 5)
flags  <- c("site_6mer","site_7mer_m8","site_7mer_A1","site_8mer","site_7mer_1mm")
labels <- c("6mer","7mer-m8","7mer-A1","8mer","7mer-1mm (any)")

# Fractions + counts for bars
frac_tbl <- function(df, group_label) {
  tibble(
    group = group_label,
    flag  = flags,
    class = labels[match(flags, flags)],
    hits  = vapply(flags, function(fl) sum(df[[fl]], na.rm=TRUE), numeric(1)),
    total = nrow(df),
    frac  = if (nrow(df) == 0) NA_real_ else hits / total
  )
}

bars <- bind_rows(
  frac_tbl(up_gene,   "Upregulated"),
  frac_tbl(down_gene, "Downregulated"),
  frac_tbl(bg_gene,   "Background")
) %>%
  mutate(
    class = factor(class, levels=labels),
    group = factor(group, levels=c("Upregulated","Downregulated","Background"))
  )

fisher_res <- bind_rows(
  map_dfr(flags, ~ fisher_for_flag(up_gene,   bg_gene, .x, "UP_vs_BG")),
  map_dfr(flags, ~ fisher_for_flag(down_gene, bg_gene, .x, "DOWN_vs_BG"))
) %>%
  mutate(q = p.adjust(p, method="BH"))

readr::write_csv(fisher_res, file.path(OUTDIR, "fisher_results_combined.csv"))

# Join p-values onto bars (only for Up/Down)
bars_annot <- bars %>%
  mutate(contrast = recode(as.character(group),
                           "Upregulated"="UP_vs_BG",
                           "Downregulated"="DOWN_vs_BG",
                           "Background"=NA_character_)) %>%
  left_join(fisher_res, by=c("contrast"="contrast","flag"="flag")) %>%
  mutate(
    p_str = if_else(is.na(p), NA_character_, fmt_p(p)),
    # annotate each bar with counts ALWAYS; add p for Up/Down
    label = case_when(
      group == "Background" ~ sprintf("n=%d/%d", hits, total),
      TRUE ~ sprintf("n=%d/%d\np=%s", hits, total, p_str)
    )
  )

# Place labels a bit above each bar
y_pad <- 0.06
p_bar <- ggplot(bars_annot, aes(x=group, y=frac, fill=class)) +
  geom_col(position=position_dodge(width=0.85)) +
  geom_text(aes(label=label, y=pmin(1, frac + y_pad)),
            position=position_dodge(width=0.85),
            vjust=0, size=3.4, lineheight=0.95, na.rm=TRUE) +
  scale_y_continuous(labels=scales::percent_format(),
                     limits=c(0, 1.20),
                     expand=expansion(mult=c(0.02, 0.12))) +
  coord_cartesian(clip="off") +
  theme_bw(base_size=12) +
  theme(plot.margin = margin(8, 28, 8, 8),
        legend.position="right") +
  labs(
    title = paste0(mirna_name, " site-class presence (canonical + 7mer-1mm)"),
    subtitle = "Bar labels show gene counts (hits/total). P-values are Fisher tests for Up vs BG_rest and Down vs BG_rest.",
    x = NULL,
    y = "Fraction of genes with ≥1 site",
    fill = NULL
  )

ggsave(file.path(OUTDIR, "barplot_site_classes_canon_plus_1mm_annotated.png"),
       p_bar, width=9.0, height=5.3, dpi=300)

# -------------------- ECDF A: CANONICAL EXCLUSIVE vs No-canonical baseline --------------------
# exclusive canonical class (ignore mismatch here): 8mer > 7m8 > 7A1 > 6mer > No canonical
canon_excl <- gene_level %>%
  transmute(
    gene_id_input,
    canon_class = case_when(
      site_8mer      ~ "8mer",
      site_7mer_m8   ~ "7mer-m8",
      site_7mer_A1   ~ "7mer-A1",
      site_6mer      ~ "6mer",
      TRUE           ~ "No canonical"
    )
  )

canon_df <- deseq_data %>%
  inner_join(canon_excl, by=c("gene_id"="gene_id_input")) %>%
  filter(is.finite(log2FoldChange)) %>%
  mutate(canon_class = factor(canon_class, levels=c("No canonical","6mer","7mer-A1","7mer-m8","8mer")))

baseline_canon <- canon_df %>% filter(canon_class=="No canonical") %>% pull(log2FoldChange)

wcx_less <- function(x, y) {
  if (length(x) < 5 || length(y) < 5) return(NA_real_)
  suppressWarnings(wilcox.test(x, y, alternative="less")$p.value)
}
# ===============================================================
# ECDF C: NON-EXCLUSIVE classes vs "No canonical" baseline
# - Non-exclusive = a gene can appear in multiple curves
#   (e.g., 8mer genes are also counted in 7mer-m8, 7mer-A1, 6mer if those flags are true)
# - Baseline (black) = No canonical site at all (no 6/7/8)
# ===============================================================

nonex_df <- gene_level %>%
  transmute(
    gene_id_input,
    site_6mer,
    site_7mer_A1,
    site_7mer_m8,
    site_8mer,
    site_any_canon
  ) %>%
  inner_join(deseq_data, by = c("gene_id_input" = "gene_id")) %>%
  filter(is.finite(log2FoldChange))

baseline_nonex <- nonex_df %>%
  filter(!site_any_canon) %>%
  pull(log2FoldChange)

# Long format: each gene contributes to every class it belongs to (non-exclusive)
nonex_long <- bind_rows(
  nonex_df %>% filter(site_6mer)      %>% mutate(site_class = "6mer (non-exclusive)"),
  nonex_df %>% filter(site_7mer_A1)   %>% mutate(site_class = "7mer-A1 (non-exclusive)"),
  nonex_df %>% filter(site_7mer_m8)   %>% mutate(site_class = "7mer-m8 (non-exclusive)"),
  nonex_df %>% filter(site_8mer)      %>% mutate(site_class = "8mer (non-exclusive)")
) %>%
  mutate(site_class = factor(site_class,
                             levels = c("6mer (non-exclusive)",
                                        "7mer-A1 (non-exclusive)",
                                        "7mer-m8 (non-exclusive)",
                                        "8mer (non-exclusive)")))

# Stats vs baseline (left shift => more downregulated)
ecdf_stats_nonex <- nonex_long %>%
  group_by(site_class) %>%
  summarize(
    n    = n(),
    dmed = median(log2FoldChange, na.rm = TRUE) - median(baseline_nonex, na.rm = TRUE),
    p    = wcx_less(log2FoldChange, baseline_nonex),
    .groups = "drop"
  ) %>%
  mutate(q = p.adjust(p, method = "BH"))

readr::write_csv(ecdf_stats_nonex,
                 file.path(OUTDIR, "ecdf_stats_canonical_nonexclusive_vs_NoCanonical.csv"))

lab_nonex <- setNames(
  sprintf("%s  Δmed=%.3f  p=%s  q=%s  (n=%d)",
          ecdf_stats_nonex$site_class,
          ecdf_stats_nonex$dmed,
          formatC(ecdf_stats_nonex$p, format="e", digits=2),
          formatC(ecdf_stats_nonex$q, format="e", digits=2),
          ecdf_stats_nonex$n),
  as.character(ecdf_stats_nonex$site_class)
)

# X-limits based on central 98% of baseline + all plotted points
all_vals_nonex <- c(baseline_nonex, nonex_long$log2FoldChange)
x_min3 <- quantile(all_vals_nonex, 0.01, na.rm = TRUE)
x_max3 <- quantile(all_vals_nonex, 0.99, na.rm = TRUE)

p_ecdf_nonex <- ggplot() +
  stat_ecdf(data = tibble(log2FoldChange = baseline_nonex),
            aes(x = log2FoldChange), linewidth = 1, color = "black") +
  stat_ecdf(data = nonex_long,
            aes(x = log2FoldChange, color = site_class), linewidth = 1) +
  coord_cartesian(xlim = c(x_min3, x_max3)) +
  scale_color_manual(
    values = c("6mer (non-exclusive)"    = "#2CA02C",
               "7mer-A1 (non-exclusive)" = "#1F77B4",
               "7mer-m8 (non-exclusive)" = "#D62728",
               "8mer (non-exclusive)"    = "#7B61FF"),
    labels = lab_nonex
  ) +
  theme_bw(base_size = 12) +
  theme(legend.text = element_text(size = 9)) +
  labs(
    title    = paste0("ECDF: canonical NON-exclusive classes vs No-canonical baseline (", mirna_name, ")"),
    subtitle = sprintf("Baseline (black): No canonical site (n=%d)", sum(!nonex_df$site_any_canon)),
    x        = "log2 fold change",
    y        = "Cumulative fraction",
    color    = NULL
  )

ggsave(file.path(OUTDIR, "ecdf_canonical_nonexclusive_vs_NoCanonical.png"),
       p_ecdf_nonex, width = 7.8, height = 5.0, dpi = 300)

ecdf_stats_canon <- canon_df %>%
  group_by(canon_class) %>%
  summarize(
    n    = sum(is.finite(log2FoldChange)),
    dmed = median(log2FoldChange, na.rm=TRUE) - median(baseline_canon, na.rm=TRUE),
    p    = wcx_less(log2FoldChange, baseline_canon),
    .groups="drop"
  ) %>%
  mutate(q = p.adjust(p, method="BH"))

readr::write_csv(ecdf_stats_canon, file.path(OUTDIR, "ecdf_stats_canonical_exclusive_vs_NoCanonical.csv"))

lab_canon <- setNames(
  sprintf("%s  Δmed=%.3f  q=%s  (n=%d)",
          ecdf_stats_canon$canon_class,
          ecdf_stats_canon$dmed,
          formatC(ecdf_stats_canon$q, format="e", digits=2),
          ecdf_stats_canon$n),
  as.character(ecdf_stats_canon$canon_class)
)

all_vals_canon <- canon_df$log2FoldChange
x_min <- quantile(all_vals_canon, 0.01, na.rm=TRUE)
x_max <- quantile(all_vals_canon, 0.99, na.rm=TRUE)

p_ecdf_canon <- ggplot() +
  stat_ecdf(data=canon_df %>% filter(canon_class=="No canonical"),
            aes(x=log2FoldChange), linewidth=1, color="black") +
  stat_ecdf(data=canon_df %>% filter(canon_class!="No canonical"),
            aes(x=log2FoldChange, color=canon_class), linewidth=1) +
  coord_cartesian(xlim=c(x_min, x_max)) +
  scale_color_manual(
    values=c("6mer"="#2CA02C","7mer-A1"="#1F77B4","7mer-m8"="#D62728","8mer"="#7B61FF"),
    labels=lab_canon[names(lab_canon) != "No canonical"]
  ) +
  theme_bw(base_size=12) +
  theme(legend.text = element_text(size=9)) +
  labs(
    title = paste0("ECDF: canonical exclusive classes vs No-canonical (", mirna_name, ")"),
    subtitle = sprintf("Baseline (black): No canonical site (n=%d)", sum(canon_df$canon_class=="No canonical")),
    x = "log2 fold change",
    y = "Cumulative fraction",
    color = NULL
  )

ggsave(file.path(OUTDIR, "ecdf_canonical_exclusive_vs_NoCanonical.png"),
       p_ecdf_canon, width=7.8, height=5.0, dpi=300)

# -------------------- ECDF B: MISMATCH CATEGORIES vs No-site baseline --------------------
# baseline = no canonical AND no 1mm
mm_merge <- gene_level %>%
  transmute(
    gene_id_input,
    site_any_canon,
    site_7mer_1mm,
    site_any_canon_or_1mm
  ) %>%
  inner_join(deseq_data, by=c("gene_id_input"="gene_id")) %>%
  filter(is.finite(log2FoldChange)) %>%
  mutate(
    mm_class = case_when(
      !site_any_canon_or_1mm                 ~ "No site (no canon & no 1mm)",
      site_any_canon                         ~ "Canonical (any)",
      !site_any_canon & site_7mer_1mm         ~ "7mer-1mm only",
      site_any_canon & site_7mer_1mm          ~ "Canonical + 1mm"
    )
  )

baseline_mm <- mm_merge %>% filter(mm_class=="No site (no canon & no 1mm)") %>% pull(log2FoldChange)

mm_stats <- mm_merge %>%
  filter(mm_class != "No site (no canon & no 1mm)") %>%
  group_by(mm_class) %>%
  summarize(
    n    = n(),
    dmed = median(log2FoldChange, na.rm=TRUE) - median(baseline_mm, na.rm=TRUE),
    p    = wcx_less(log2FoldChange, baseline_mm),
    .groups="drop"
  ) %>%
  mutate(q = p.adjust(p, method="BH"))

readr::write_csv(mm_stats, file.path(OUTDIR, "ecdf_stats_mismatch_categories_vs_NoSite.csv"))

lab_mm <- setNames(
  sprintf("%s  Δmed=%.3f  p=%s  q=%s  (n=%d)",
          mm_stats$mm_class,
          mm_stats$dmed,
          formatC(mm_stats$p, format="e", digits=2),
          formatC(mm_stats$q, format="e", digits=2),
          mm_stats$n),
  mm_stats$mm_class
)

# build long df for plotting
mm_long <- mm_merge %>%
  mutate(mm_class = factor(mm_class,
                           levels=c("No site (no canon & no 1mm)",
                                    "7mer-1mm only",
                                    "Canonical (any)",
                                    "Canonical + 1mm")))

all_vals_mm <- mm_long$log2FoldChange
x_min2 <- quantile(all_vals_mm, 0.01, na.rm=TRUE)
x_max2 <- quantile(all_vals_mm, 0.99, na.rm=TRUE)

p_ecdf_mm <- ggplot() +
  stat_ecdf(data=mm_long %>% filter(mm_class=="No site (no canon & no 1mm)"),
            aes(x=log2FoldChange), linewidth=1, color="black") +
  stat_ecdf(data=mm_long %>% filter(mm_class!="No site (no canon & no 1mm)"),
            aes(x=log2FoldChange, color=mm_class), linewidth=1) +
  coord_cartesian(xlim=c(x_min2, x_max2)) +
  scale_color_manual(
    values=c("7mer-1mm only"="#1F77B4",
             "Canonical (any)"="#D62728",
             "Canonical + 1mm"="#7B61FF"),
    labels=lab_mm
  ) +
  theme_bw(base_size=12) +
  theme(legend.text = element_text(size=9)) +
  labs(
    title = paste0("ECDF: mismatch categories vs No-site baseline (", mirna_name, ")"),
    subtitle = sprintf("Baseline (black): no canonical & no 1mm (n=%d)",
                       sum(mm_long$mm_class=="No site (no canon & no 1mm)")),
    x = "log2 fold change",
    y = "Cumulative fraction",
    color = NULL
  )

ggsave(file.path(OUTDIR, "ecdf_mismatch_categories_vs_NoSite.png"),
       p_ecdf_mm, width=7.8, height=5.0, dpi=300)

# -------------------- SESSION INFO --------------------
sink(file.path(OUTDIR, "session_info.txt"))
cat("miRNA: ", mirna_name, "\n")
cat("miRNA seq: ", mirna_seq, "\n\n")
cat("DESeq file: ", deseq_file, "\n")
cat("Ensembl version: ", ensembl_ver, "\n")
cat("Dataset: ", species_ds, "\n")
cat("Input IDs treated as: ", if (input_are_ensembl) "ENSEMBL" else "SYMBOLS", "\n\n")
sessionInfo()
sink()

message("All done. Outputs written to: ", OUTDIR)

