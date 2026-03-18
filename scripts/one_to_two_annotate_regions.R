#!/usr/bin/env Rscript
# ============================================================
# Annotate per-site TSVs with region label (UTR3/CDS/UTR5/INTRON/OTHER)
# using a mm10 GTF (e.g., gencode.vM10.annotation.gtf).
#
# Filters applied:
#   - feature_name is present (not NA/""/"."), and
#   - site is A>G based on ref/alt columns (ref=="A" & alt=="G")
#
# INTRON definition (fast):
#   - per gene_name, per chromosome:
#       exon blocks are reduced
#       gene body = [min exon start, max exon end] on that chromosome
#       introns = gaps between reduced exon blocks within that span
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(IRanges)
  library(S4Vectors)
})

# ---------------------------
# USER SETTINGS
# ---------------------------
sites_dir   <- "/abe_ago_200_465_analysis/one_filter_AtoG"
gtf_path    <- "/abe_ago_200_465_analysis/gencode.vM10.annotation.gtf"
out_dir     <- "/abe_ago_200_465_analysis/two_annotate_by_region"
summary_csv <- file.path(out_dir, "annotation_summary.csv")

# Your example position looks 1-based already
coord_is_0based <- FALSE

# Region priority if multiple overlaps occur (leftmost wins)
region_priority <- c("UTR3", "CDS", "UTR5", "INTRON", "OTHER")

# ---------------------------
# Helpers
# ---------------------------
reduce_gr <- GenomicRanges::reduce

harmonize_chr_gr <- function(gr) {
  # add "chr" prefix if absent; chrMT -> chrM
  if (!any(grepl("^chr", seqlevels(gr)))) {
    old_levels <- seqlevels(gr)
    new_levels <- paste0("chr", old_levels)
    gr <- renameSeqlevels(gr, setNames(new_levels, old_levels))
  }
  if ("chrMT" %in% seqlevels(gr)) {
    gr <- renameSeqlevels(gr, c(chrMT = "chrM"))
  }
  gr
}

harmonize_contig_vec <- function(x) {
  x <- as.character(x)
  x <- ifelse(grepl("^chr", x), x, paste0("chr", x))
  x <- sub("^chrMT$", "chrM", x)
  x
}

hits_any <- function(query_gr, subject_gr) {
  if (length(subject_gr) == 0L || length(query_gr) == 0L) {
    return(rep(FALSE, length(query_gr)))
  }
  h <- GenomicRanges::findOverlaps(query_gr, subject_gr, ignore.strand = TRUE)
  if (length(h) == 0L) return(rep(FALSE, length(query_gr)))
  out <- rep(FALSE, length(query_gr))
  out[unique(S4Vectors::queryHits(h))] <- TRUE
  out
}

assign_by_priority <- function(hit_named, priority) {
  for (lab in priority) {
    if (!is.na(hit_named[lab]) && isTRUE(hit_named[lab])) return(lab)
  }
  "OTHER"
}

make_utr_from_generic <- function(gtf) {
  gtf_type <- as.character(gtf$type)
  
  utr_explicit_3 <- gtf[gtf_type %in% c("three_prime_utr", "3UTR", "3'UTR", "UTR3")]
  utr_explicit_5 <- gtf[gtf_type %in% c("five_prime_utr",  "5UTR", "5'UTR", "UTR5")]
  
  if (length(utr_explicit_3) > 0 || length(utr_explicit_5) > 0) {
    return(list(utr3 = utr_explicit_3, utr5 = utr_explicit_5))
  }
  
  utr_generic <- gtf[gtf_type == "UTR"]
  if (length(utr_generic) == 0) {
    return(list(utr3 = GRanges(), utr5 = GRanges()))
  }
  
  tx <- mcols(utr_generic)$transcript_id
  cds <- gtf[gtf_type == "CDS"]
  cds_tx <- mcols(cds)$transcript_id
  
  if (is.null(tx) || is.null(cds_tx)) {
    warning("GTF UTRs are generic and transcript_id missing; cannot split UTR5 vs UTR3. Returning empty UTR5/UTR3.")
    return(list(utr3 = GRanges(), utr5 = GRanges()))
  }
  
  tx <- as.character(tx)
  cds_tx <- as.character(cds_tx)
  
  cds <- cds[!is.na(cds_tx) & cds_tx != ""]
  cds_tx <- as.character(mcols(cds)$transcript_id)
  
  cds_dt <- data.table(
    transcript_id = cds_tx,
    cds_start = start(cds),
    cds_end   = end(cds),
    strand    = as.character(strand(cds))
  )
  cds_bounds <- cds_dt[, .(
    cds_min = min(cds_start),
    cds_max = max(cds_end),
    strand  = strand[1]
  ), by = transcript_id]
  
  utr_dt <- data.table(
    idx = seq_along(utr_generic),
    transcript_id = as.character(mcols(utr_generic)$transcript_id),
    u_start = start(utr_generic),
    u_end   = end(utr_generic)
  )
  utr_dt <- utr_dt[!is.na(transcript_id) & transcript_id != ""]
  utr_dt <- merge(utr_dt, cds_bounds, by = "transcript_id", all.x = TRUE)
  
  utr_dt[, utr_class := NA_character_]
  utr_dt[strand == "+" & !is.na(cds_min) & u_end   < cds_min, utr_class := "UTR5"]
  utr_dt[strand == "+" & !is.na(cds_max) & u_start > cds_max, utr_class := "UTR3"]
  utr_dt[strand == "-" & !is.na(cds_max) & u_start > cds_max, utr_class := "UTR5"]
  utr_dt[strand == "-" & !is.na(cds_min) & u_end   < cds_min, utr_class := "UTR3"]
  
  utr5_idx <- utr_dt[utr_class == "UTR5", idx]
  utr3_idx <- utr_dt[utr_class == "UTR3", idx]
  
  list(
    utr3 = utr_generic[utr3_idx],
    utr5 = utr_generic[utr5_idx]
  )
}

build_introns_gene_chr <- function(exon_gr) {
  if (length(exon_gr) == 0) return(GRanges())
  
  if (!("gene_name" %in% names(mcols(exon_gr)))) {
    stop("Exon GRanges lacks gene_name; cannot derive introns.")
  }
  
  gene <- as.character(mcols(exon_gr)$gene_name)
  chr  <- as.character(seqnames(exon_gr))
  keep <- !is.na(gene) & gene != "" & !is.na(chr) & chr != ""
  exon_gr <- exon_gr[keep]
  gene <- gene[keep]
  chr  <- chr[keep]
  
  exon2 <- exon_gr
  mcols(exon2) <- NULL
  strand(exon2) <- "*"
  
  key <- paste(gene, chr, sep = "||")
  exon2$key <- key
  exon_by_key <- split(exon2, exon2$key)
  
  exon_by_key_red <- S4Vectors::endoapply(exon_by_key, function(x) reduce_gr(x, ignore.strand = TRUE))
  
  intr_list <- S4Vectors::endoapply(exon_by_key_red, function(ex) {
    if (length(ex) == 0) return(GRanges())
    body <- range(ex)
    intr <- gaps(ex, start = start(body), end = end(body))
    intr[width(intr) > 0]
  })
  
  intr_gr <- unlist(intr_list, use.names = FALSE)
  if (length(intr_gr) == 0) return(GRanges())
  reduce_gr(intr_gr, ignore.strand = TRUE)
}

# ---------------------------
# 1) Load GTF and make feature ranges
# ---------------------------
message("[1/4] Importing GTF: ", gtf_path)
gtf <- rtracklayer::import(gtf_path)
gtf <- harmonize_chr_gr(gtf)

gtf_type <- as.character(gtf$type)

cds_gr <- gtf[gtf_type == "CDS"]

utr_parts <- make_utr_from_generic(gtf)
utr3_gr <- utr_parts$utr3
utr5_gr <- utr_parts$utr5

message("  CDS intervals:   ", length(cds_gr))
message("  3'UTR intervals: ", length(utr3_gr))
message("  5'UTR intervals: ", length(utr5_gr))

if (length(cds_gr) == 0) stop("No CDS features found in GTF.")

exon_gr <- gtf[gtf_type == "exon"]
if (length(exon_gr) == 0) stop("No exon features found in GTF; cannot derive introns.")
if (!("gene_name" %in% names(mcols(exon_gr)))) stop("GTF exons missing gene_name; cannot derive introns.")

message("[1/4] Deriving introns from exons (gene_name + chr; fast)...")
intron_red <- build_introns_gene_chr(exon_gr)

cds_red  <- reduce_gr(cds_gr,  ignore.strand = TRUE)
utr3_red <- if (length(utr3_gr) > 0) reduce_gr(utr3_gr, ignore.strand = TRUE) else GRanges()
utr5_red <- if (length(utr5_gr) > 0) reduce_gr(utr5_gr, ignore.strand = TRUE) else GRanges()

message("  Reduced intervals: CDS=", length(cds_red),
        " | UTR3=", length(utr3_red),
        " | UTR5=", length(utr5_red),
        " | INTRON=", length(intron_red))

# ---------------------------
# 2) Iterate over TSV files and annotate (FILTERED)
# ---------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(sites_dir, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE)
if (length(files) == 0) stop("No .tsv or .tsv.gz files found in: ", sites_dir)

message("[2/4] Found ", length(files), " files to process.")

summary_dt <- data.table(
  file = character(),
  n_rows_in = integer(),
  n_rows_after_filter = integer(),
  n_valid_coords = integer(),
  n_CDS = integer(),
  n_UTR3 = integer(),
  n_UTR5 = integer(),
  n_INTRON = integer(),
  n_OTHER = integer(),
  out_file = character()
)

for (f in files) {
  message("\nProcessing: ", basename(f))
  dt <- fread(f)
  
  # ---- REQUIRED COLUMNS (YOUR ACTUAL HEADERS) ----
  required <- c("contig", "position", "feature_name", "ref", "alt", "count", "coverage")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("Missing required columns in ", basename(f), ": ", paste(missing, collapse = ", "))
  }
  
  n_in <- nrow(dt)
  
  # ---- FILTER: feature_name present AND A>G ----
  dt <- dt[
    !is.na(feature_name) & feature_name != "" & feature_name != "." &
      !is.na(ref) & !is.na(alt) &
      toupper(ref) == "A" & toupper(alt) == "G"
  ]
  n_after <- nrow(dt)
  
  if (n_after == 0) {
    warning("  All rows removed by filter (feature_name present & ref=A & alt=G): ", basename(f))
    out_file <- file.path(out_dir, basename(f))
    fwrite(dt, out_file, sep = "\t")
    summary_dt <- rbind(summary_dt, data.table(
      file = basename(f),
      n_rows_in = n_in,
      n_rows_after_filter = n_after,
      n_valid_coords = 0L,
      n_CDS=0L, n_UTR3=0L, n_UTR5=0L, n_INTRON=0L, n_OTHER=0L,
      out_file = out_file
    ))
    next
  }
  
  # Harmonize contigs to chr*
  dt[, contig := harmonize_contig_vec(contig)]
  
  dt[, row_id := .I]
  dt[, region := NA_character_]
  dt[, region_source := NA_character_]
  
  chr <- as.character(dt$contig)
  pos <- suppressWarnings(as.integer(dt$position))
  if (isTRUE(coord_is_0based)) pos <- pos + 1L
  
  ok <- !is.na(chr) & !is.na(pos) & pos >= 1L
  n_ok <- sum(ok)
  
  if (n_ok == 0) {
    warning("  No valid coordinates after filter; writing with region=NA: ", basename(f))
    out_file <- file.path(out_dir, basename(f))
    fwrite(dt[, !"row_id"], out_file, sep = "\t")
    summary_dt <- rbind(summary_dt, data.table(
      file = basename(f),
      n_rows_in = n_in,
      n_rows_after_filter = n_after,
      n_valid_coords = n_ok,
      n_CDS=0L, n_UTR3=0L, n_UTR5=0L, n_INTRON=0L, n_OTHER=0L,
      out_file = out_file
    ))
    next
  }
  
  dt_ok <- dt[ok, .(row_id, contig, position)]
  gr_sites <- GRanges(
    seqnames = dt_ok$contig,
    ranges   = IRanges(start = dt_ok$position, end = dt_ok$position)
  )
  
  gtf_levels <- unique(c(
    seqlevels(cds_red),
    seqlevels(utr3_red),
    seqlevels(utr5_red),
    seqlevels(intron_red)
  ))
  common <- intersect(unique(as.character(seqnames(gr_sites))), gtf_levels)
  
  if (length(common) == 0) {
    stop(
      "No matching chromosome names between TSV and GTF for file ", basename(f), ".\n",
      "Example TSV contig: ", paste(head(unique(as.character(seqnames(gr_sites))), 3), collapse = ", "), "\n",
      "Example GTF seqlevels: ", paste(head(gtf_levels, 3), collapse = ", "), "\n",
      "Fix by making both use chr* or no-chr consistently."
    )
  }
  
  gr_sites2 <- keepSeqlevels(gr_sites, common, pruning.mode = "coarse")
  keep_ok2 <- dt_ok$contig %in% common
  dt_ok2 <- dt_ok[keep_ok2]
  
  cds_use    <- keepSeqlevels(cds_red,    common, pruning.mode = "coarse")
  utr3_use   <- if (length(utr3_red)   > 0) keepSeqlevels(utr3_red,   common, pruning.mode = "coarse") else GRanges()
  utr5_use   <- if (length(utr5_red)   > 0) keepSeqlevels(utr5_red,   common, pruning.mode = "coarse") else GRanges()
  intron_use <- if (length(intron_red) > 0) keepSeqlevels(intron_red, common, pruning.mode = "coarse") else GRanges()
  
  hit_cds    <- hits_any(gr_sites2, cds_use)
  hit_utr3   <- hits_any(gr_sites2, utr3_use)
  hit_utr5   <- hits_any(gr_sites2, utr5_use)
  hit_intron <- hits_any(gr_sites2, intron_use)
  
  region_vec <- vapply(seq_along(gr_sites2), function(i) {
    assign_by_priority(
      c(CDS = hit_cds[i], UTR3 = hit_utr3[i], UTR5 = hit_utr5[i], INTRON = hit_intron[i], OTHER = TRUE),
      region_priority
    )
  }, character(1))
  
  annot <- data.table(row_id = dt_ok2$row_id, region = region_vec)
  dt[annot, on = "row_id", region := i.region]
  dt[ok, region_source := "GTF_overlap"]
  
  out_file <- file.path(out_dir, basename(f))
  fwrite(dt[, !"row_id"], out_file, sep = "\t")
  
  n_cds    <- sum(dt$region == "CDS",    na.rm = TRUE)
  n_utr3   <- sum(dt$region == "UTR3",   na.rm = TRUE)
  n_utr5   <- sum(dt$region == "UTR5",   na.rm = TRUE)
  n_intron <- sum(dt$region == "INTRON", na.rm = TRUE)
  n_other  <- sum(dt$region == "OTHER",  na.rm = TRUE)
  
  summary_dt <- rbind(
    summary_dt,
    data.table(
      file = basename(f),
      n_rows_in = n_in,
      n_rows_after_filter = n_after,
      n_valid_coords = n_ok,
      n_CDS = n_cds,
      n_UTR3 = n_utr3,
      n_UTR5 = n_utr5,
      n_INTRON = n_intron,
      n_OTHER = n_other,
      out_file = out_file
    )
  )
  
  message("  total rows (in): ", n_in,
          " | after filter: ", n_after,
          " | valid coords: ", n_ok,
          " | CDS=", n_cds,
          " | UTR3=", n_utr3,
          " | UTR5=", n_utr5,
          " | INTRON=", n_intron,
          " | OTHER=", n_other,
          " -> ", out_file)
}

fwrite(summary_dt, summary_csv)
message("\n[4/4] Done.")
message("Output folder: ", out_dir)
message("Summary CSV: ", summary_csv)
