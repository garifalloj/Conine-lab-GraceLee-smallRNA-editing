#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

# ------------------ USER EDITS ------------------
INDIR  <- "/abe_ago_200_465_analysis/zero_MARINE_outputs"
OUTDIR_AG <- file.path("/abe_ago_200_465_analysis/one_filter_AtoG")
dir.create(OUTDIR_AG, showWarnings = FALSE, recursive = TRUE)

FILES <- list.files(INDIR, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE)
stopifnot(length(FILES) > 0)

# Optional: keep consistent with earlier logic (set FALSE if you don't want this filter)
REQUIRE_STRAND_CONVERSION_NONEMPTY <- TRUE
# ------------------------------------------------

read_one <- function(file) {
  suppressWarnings(readr::read_tsv(
    file,
    col_types = cols(
      site_id = col_character(),
      barcode = col_character(),
      contig = col_character(),
      position = col_double(),
      ref = col_character(),
      alt = col_character(),
      strand = col_character(),
      count = col_double(),
      coverage = col_double(),
      conversion = col_double(),
      strand_conversion = col_character(),
      feature_name = col_character(),
      feature_type = col_character(),
      feature_strand = col_character()
    ),
    progress = FALSE
  ))
}

filter_criteria <- function(df) {
  out <- df %>%
    dplyr::filter(
      !is.na(feature_name),
      stringr::str_trim(feature_name) != "."
    )
  
  if (REQUIRE_STRAND_CONVERSION_NONEMPTY) {
    out <- out %>%
      dplyr::filter(!is.na(strand_conversion), stringr::str_trim(strand_conversion) != ".")
  }
  
  # Keep only A>G edits
  out %>%
    dplyr::filter(
      !is.na(ref), !is.na(alt),
      toupper(ref) == "A",
      toupper(alt) == "G"
    )
}

write_filtered_file <- function(file) {
  df <- read_one(file)
  df_filt <- filter_criteria(df)
  
  base <- base::basename(file)
  base_noext <- stringr::str_replace(base, "\\.tsv(\\.gz)?$", "")
  out_path <- file.path(OUTDIR_AG, paste0(base_noext, "_AtoG.tsv"))
  
  readr::write_tsv(df_filt, out_path)
  tibble::tibble(file = base, n_rows_written = nrow(df_filt), out_path = out_path)
}

message("Filtering & writing ", length(FILES), " files...")
log_df <- purrr::map_dfr(FILES, write_filtered_file)

readr::write_csv(log_df, file.path(OUTDIR_AG, "write_log.csv"))
message("Done. Filtered files written to:\n  ", OUTDIR_AG)
