#!/usr/bin/env Rscript
# ============================================================
# Compute RNA-editing EPKM per gene for:
#   UTR3, CDS, INTRON   (+ OTHER tracked, no EPKM)
#
# Hardened for "quirky" Bioconductor installs:
#  - Filters out NA/empty gene_name + transcript_id early
#  - Avoids overlapsAny() (uses findOverlaps-based helper)
#  - Prevents NA gene IDs from entering GRanges lists / SAF
#  - Guards empty outputs more explicitly
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(Rsubread)
  library(S4Vectors)
})

# ---------------------------
# USER SETTINGS (EDIT THESE)
# ---------------------------
sites_dir <- "/abe_ago_200_465_analysis/two_annotate_by_region"
bam_dir   <- "/20251028_GL_single_embryo_abe_ago_miR465_miR200_2cell/rsem"
gtf_path  <- "/abe_ago_200_465_analysis/gencode.vM10.annotation.gtf"

out_dir <- file.path("/abe_ago_200_465_analysis/three_EPKM_outputs_UTR3_CDS_INTRON")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# featureCounts settings
is_paired <- TRUE
strandSpecific <- 0   # 0=unstranded, 1=forward, 2=reverse
nthreads <- 8

# per-site coordinate system
coord_is_0based <- FALSE

# optional stability filter
min_mapped_reads_region <- 0

# Prefer TSV region column if present
use_tsv_region_if_present <- TRUE

# ---------------------------
# Helpers
# ---------------------------
normalize_sample_id <- function(x) {
  x <- basename(x)
  x <- sub("\\.tsv(\\.gz)?$", "", x)
  x <- sub("\\.bam$", "", x)
  x <- sub("^rsem\\.out\\.", "", x)
  x <- sub("\\.genome\\.sorted$", "", x)
  x <- sub("(_S\\d+_L\\d+).*$", "", x)
  x <- sub("_filtered_site.*$", "", x)
  x
}

harmonize_contig_vec <- function(x, want_chr_prefix=TRUE) {
  x <- as.character(x)
  if (want_chr_prefix) {
    need <- !grepl("^chr", x)
    x[need] <- paste0("chr", x[need])
  }
  x <- sub("^chrMT$", "chrM", x)
  x <- sub("^MT$", "chrM", x)
  x
}

# findOverlaps-based "overlapsAny" for compatibility
overlaps_any_safe <- function(query_gr, subject_gr, ignore.strand=TRUE) {
  if (length(query_gr) == 0) return(logical(0))
  if (length(subject_gr) == 0) return(rep(FALSE, length(query_gr)))
  hits <- GenomicRanges::findOverlaps(query_gr, subject_gr, ignore.strand=ignore.strand)
  out <- rep(FALSE, length(query_gr))
  out[unique(S4Vectors::queryHits(hits))] <- TRUE
  out
}

# Reduce per gene_name using index splitting (safe)
reduce_by_gene_safe <- function(gr) {
  if (length(gr) == 0) return(list())
  mcn <- names(S4Vectors::mcols(gr))
  if (!("gene_name" %in% mcn)) stop("GRanges missing gene_name mcol.")
  gene <- as.character(S4Vectors::mcols(gr)[, "gene_name"])
  
  keep <- !is.na(gene) & nzchar(gene) & gene != "."
  gr <- gr[keep]
  gene <- gene[keep]
  if (length(gr) == 0) return(list())
  
  S4Vectors::mcols(gr) <- NULL
  GenomicRanges::strand(gr) <- "*"
  idx_by_gene <- split(seq_along(gr), gene)
  lapply(idx_by_gene, function(idx) GenomicRanges::reduce(gr[idx], ignore.strand = TRUE))
}

# Reduce per transcript_id using index splitting (safe)
reduce_by_tx_safe <- function(gr) {
  if (length(gr) == 0) return(list())
  mcn <- names(S4Vectors::mcols(gr))
  if (!("transcript_id" %in% mcn)) stop("GRanges missing transcript_id mcol.")
  tx <- as.character(S4Vectors::mcols(gr)[, "transcript_id"])
  
  keep <- !is.na(tx) & nzchar(tx) & tx != "."
  gr <- gr[keep]
  tx <- tx[keep]
  if (length(gr) == 0) return(list())
  
  idx_by_tx <- split(seq_along(gr), tx)
  lapply(idx_by_tx, function(idx) GenomicRanges::reduce(gr[idx], ignore.strand = TRUE))
}

# Global reduced ranges from list-of-GRanges (robust; never returns a list)
global_reduce_from_list <- function(gr_list) {
  if (length(gr_list) == 0) return(GenomicRanges::GRanges())
  
  acc <- GenomicRanges::GRanges()
  
  for (i in seq_along(gr_list)) {
    g <- gr_list[[i]]
    if (is.null(g) || length(g) == 0) next
    acc <- base::c(acc, g)  # <- key change
  }
  
  if (length(acc) == 0) return(GenomicRanges::GRanges())
  GenomicRanges::reduce(acc, ignore.strand = TRUE)
}

# lengths per gene per region (bp), list-of-GRanges
len_dt <- function(gr_list, region) {
  if (length(gr_list) == 0) return(data.table(gene=character(), region=character(), region_bp=numeric()))
  lens <- vapply(gr_list, function(x) {
    if (is.null(x) || length(x) == 0) return(0)
    sum(GenomicRanges::width(x))
  }, numeric(1))
  # remove NA/empty gene IDs
  genes <- names(lens)
  keep <- !is.na(genes) & nzchar(genes) & genes != "."
  lens <- lens[keep]
  data.table(gene = names(lens), region = region, region_bp = as.numeric(lens))
}

# SAF writer for list-of-GRanges
write_saf <- function(gr_list, saf_path) {
  if (length(gr_list) == 0) stop("SAF would be empty: ", saf_path)
  # filter bad names
  nm <- names(gr_list)
  keep <- !is.na(nm) & nzchar(nm) & nm != "."
  gr_list <- gr_list[keep]
  if (length(gr_list) == 0) stop("SAF would be empty after dropping NA/empty gene IDs: ", saf_path)
  
  dt_list <- lapply(names(gr_list), function(gene) {
    g <- gr_list[[gene]]
    if (is.null(g) || length(g) == 0) return(NULL)
    data.table(
      GeneID = gene,
      Chr    = as.character(GenomicRanges::seqnames(g)),
      Start  = GenomicRanges::start(g),
      End    = GenomicRanges::end(g),
      Strand = as.character(GenomicRanges::strand(g))
    )
  })
  
  dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  if (nrow(dt) == 0) stop("SAF would be empty after filtering: ", saf_path)
  fwrite(dt, saf_path, sep="\t")
  saf_path
}

normalize_region_label <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("3'UTR","UTR3","THREE_PRIME_UTR","THREEPRIMEUTR","THREEPRIME_UTR")] <- "UTR3"
  x[x %in% c("CDS","CODING","CODING_SEQUENCE","PROTEIN_CODING_CDS")]            <- "CDS"
  x[x %in% c("INTRON","INTRONIC")]                                              <- "INTRON"
  x
}

# fallback overlap priority: UTR3 > CDS > INTRON > OTHER
pick_region_priority <- function(hit_utr3, hit_cds, hit_intron) {
  region <- rep("OTHER", length(hit_cds))
  region[hit_intron] <- "INTRON"
  region[hit_cds]    <- "CDS"
  region[hit_utr3]   <- "UTR3"
  region
}

# ---------------------------
# [1/6] Import GTF; derive CDS, UTR3, EXON, GENE; build INTRON = gene - exon
# ---------------------------
message("[1/6] Importing GTF: ", gtf_path)
gtf <- rtracklayer::import(gtf_path)

mcols_names <- names(S4Vectors::mcols(gtf))
if (!("gene_name" %in% mcols_names)) stop("GTF missing gene_name column.")
if (!("transcript_id" %in% mcols_names)) stop("GTF missing transcript_id column.")

# drop rows without gene_name upfront (prevents NA names later)
gn <- as.character(S4Vectors::mcols(gtf)[, "gene_name"])
keep_gn <- !is.na(gn) & nzchar(gn) & gn != "."
gtf <- gtf[keep_gn]

# detect whether GTF uses chr prefix
gtf_has_chr <- any(grepl("^chr", as.character(GenomicRanges::seqnames(gtf))))
message("GTF chr-prefix: ", gtf_has_chr)

# Extract features
cds_gr  <- gtf[gtf$type == "CDS"]
utr_gr  <- gtf[gtf$type == "UTR"]
exon_gr <- gtf[gtf$type == "exon"]
gene_gr <- gtf[gtf$type == "gene"]

message("Raw features: CDS=", length(cds_gr),
        " UTR=", length(utr_gr),
        " exon=", length(exon_gr),
        " gene=", length(gene_gr))

# Reduce per transcript for CDS and UTR
cds_by_tx <- reduce_by_tx_safe(cds_gr)
utr_by_tx <- reduce_by_tx_safe(utr_gr)

utr_tx_id <- as.character(S4Vectors::mcols(utr_gr)[, "transcript_id"])
keep_u <- !is.na(utr_tx_id) & nzchar(utr_tx_id) & utr_tx_id != "."
utr_tx_id <- utr_tx_id[keep_u]
utr_gr_idx <- which(keep_u)
utr_tx_idx <- split(utr_gr_idx, utr_tx_id)

cds_tx_id <- as.character(S4Vectors::mcols(cds_gr)[, "transcript_id"])
keep_c <- !is.na(cds_tx_id) & nzchar(cds_tx_id) & cds_tx_id != "."
cds_tx_id <- cds_tx_id[keep_c]
cds_gr_idx <- which(keep_c)
cds_tx_idx <- split(cds_gr_idx, cds_tx_id)

# tx -> gene_name map (drop NA/empty)
tx_gene_map <- unique(data.table(
  transcript_id = c(
    as.character(S4Vectors::mcols(cds_gr)[, "transcript_id"]),
    as.character(S4Vectors::mcols(utr_gr)[, "transcript_id"])
  ),
  gene_name = c(
    as.character(S4Vectors::mcols(cds_gr)[, "gene_name"]),
    as.character(S4Vectors::mcols(utr_gr)[, "gene_name"])
  )
))
tx_gene_map <- tx_gene_map[!is.na(transcript_id) & nzchar(transcript_id) & transcript_id != "."]
tx_gene_map <- tx_gene_map[!is.na(gene_name) & nzchar(gene_name) & gene_name != "."]
tx_gene_map <- tx_gene_map[!duplicated(transcript_id)]

# Build UTR3 per transcript: UTR segments strictly 3' of CDS bounds (strand-aware)
utr3_list <- list()
common_txs <- intersect(names(utr_by_tx), names(cds_by_tx))

if (length(common_txs) == 0) {
  message("WARNING: No transcripts have both UTR and CDS. UTR3 will be empty.")
}

for (tx in common_txs) {
  u <- utr_by_tx[[tx]]
  cds_tx <- cds_by_tx[[tx]]
  if (length(u) == 0 || length(cds_tx) == 0) next
  
  cds_start <- min(GenomicRanges::start(cds_tx))
  cds_end   <- max(GenomicRanges::end(cds_tx))
  
  # strand inference (prefer UTR rows, else CDS)
  s <- character()
  idxu <- utr_tx_idx[[tx]]
  if (!is.null(idxu) && length(idxu) > 0) {
    s <- unique(as.character(GenomicRanges::strand(utr_gr[idxu])))
  } else {
    idxc <- cds_tx_idx[[tx]]
    if (!is.null(idxc) && length(idxc) > 0) {
      s <- unique(as.character(GenomicRanges::strand(cds_gr[idxc])))
    }
  }
  s <- s[s %in% c("+","-")]
  if (length(s) == 0) next
  s <- s[1]
  
  if (s == "+") {
    u3 <- u[GenomicRanges::start(u) > cds_end]
  } else {
    u3 <- u[GenomicRanges::end(u) < cds_start]
  }
  
  if (length(u3) > 0) {
    utr3_list[[tx]] <- GenomicRanges::reduce(u3, ignore.strand = TRUE)
  }
}

# Attach gene_name to UTR3 transcript pieces and combine
utr3_pieces <- lapply(names(utr3_list), function(tx) {
  gr_piece <- utr3_list[[tx]]
  g <- tx_gene_map[match(tx, transcript_id), gene_name]
  if (is.na(g) || !nzchar(g) || g == ".") return(NULL)
  S4Vectors::mcols(gr_piece)$gene_name <- g
  gr_piece
})
utr3_pieces <- utr3_pieces[!vapply(utr3_pieces, is.null, logical(1))]
utr3_gr <- GenomicRanges::GRanges()
if (length(utr3_pieces) > 0) {
  for (k in seq_along(utr3_pieces)) {
    utr3_gr <- c(utr3_gr, utr3_pieces[[k]])
  }
}
message("Derived UTR3 pieces: ", length(utr3_gr))

# Reduce per gene (safe)
cds_by_gene   <- reduce_by_gene_safe(cds_gr)
utr3_by_gene  <- reduce_by_gene_safe(utr3_gr)
exon_by_gene  <- reduce_by_gene_safe(exon_gr)
gene_by_gene  <- reduce_by_gene_safe(gene_gr)

# Introns = gene body - exons
message("Computing introns: gene - exon ...")
genes <- names(gene_by_gene)
intron_by_gene <- vector("list", length(genes))
names(intron_by_gene) <- genes

for (g in genes) {
  gb <- gene_by_gene[[g]]
  ex <- if (g %in% names(exon_by_gene)) exon_by_gene[[g]] else GenomicRanges::GRanges()
  ex <- GenomicRanges::reduce(ex, ignore.strand = TRUE)
  intr <- GenomicRanges::setdiff(gb, ex, ignore.strand = TRUE)
  intron_by_gene[[g]] <- intr
}

message("Genes with regions:",
        " CDS=", length(cds_by_gene),
        " UTR3=", length(utr3_by_gene),
        " INTRON=", length(intron_by_gene))

# Length denominators
gene_len <- rbind(
  len_dt(utr3_by_gene,  "UTR3"),
  len_dt(cds_by_gene,   "CDS"),
  len_dt(intron_by_gene,"INTRON")
)
fwrite(gene_len, file.path(out_dir, "gene_region_lengths_bp.csv"))

# Global reduced ranges for optional fallback overlap labeling
utr3_red   <- if (length(utr3_by_gene)   > 0) global_reduce_from_list(utr3_by_gene)   else GenomicRanges::GRanges()
cds_red    <- if (length(cds_by_gene)    > 0) global_reduce_from_list(cds_by_gene)    else GenomicRanges::GRanges()
intron_red <- if (length(intron_by_gene) > 0) global_reduce_from_list(intron_by_gene) else GenomicRanges::GRanges()

# ---------------------------
# [2/6] featureCounts denominators per region
# ---------------------------
message("[2/6] Running featureCounts for region-mapped reads (UTR3/CDS/INTRON)")

bam_files <- list.files(bam_dir, pattern="genome\\.sorted\\.bam$", recursive=TRUE, full.names=TRUE)
if (length(bam_files) == 0) stop("No genome.sorted.bam files found under: ", bam_dir)
message("Found BAMs: ", length(bam_files))

saf_utr3   <- file.path(out_dir, "UTR3.saf")
saf_cds    <- file.path(out_dir, "CDS.saf")
saf_intron <- file.path(out_dir, "INTRON.saf")

if (length(utr3_by_gene)   > 0) write_saf(utr3_by_gene,   saf_utr3)
if (length(cds_by_gene)    > 0) write_saf(cds_by_gene,    saf_cds)
if (length(intron_by_gene) > 0) write_saf(intron_by_gene, saf_intron)

run_fc_saf <- function(saf_path, tag) {
  message("  featureCounts SAF: ", tag)
  fc <- featureCounts(
    files = bam_files,
    annot.ext = saf_path,
    isGTFAnnotationFile = FALSE,
    isPairedEnd = is_paired,
    strandSpecific = strandSpecific,
    nthreads = nthreads
  )
  colnames(fc$counts) <- sub("\\.genome\\.sorted\\.bam$", "", basename(colnames(fc$counts)))
  list(counts = fc$counts, stat = fc$stat)
}

fc_list <- list()
if (file.exists(saf_utr3))   fc_list$UTR3   <- run_fc_saf(saf_utr3,   "UTR3")
if (file.exists(saf_cds))    fc_list$CDS    <- run_fc_saf(saf_cds,    "CDS")
if (file.exists(saf_intron)) fc_list$INTRON <- run_fc_saf(saf_intron, "INTRON")
if (length(fc_list) == 0) stop("No SAFs were created; cannot run featureCounts.")

samples_raw <- colnames(fc_list[[1]]$counts)
mapped_dt <- data.table(sample_raw = samples_raw)
for (reg in names(fc_list)) mapped_dt[, (reg) := as.numeric(colSums(fc_list[[reg]]$counts))]
mapped_dt[, sample := normalize_sample_id(sample_raw)]
mapped_dt[, sample_raw := NULL]
fwrite(mapped_dt, file.path(out_dir, "mapped_reads_per_region_per_sample.csv"))

# ---------------------------
# [3/6] Numerators: edit counts from TSVs by gene+region
# ---------------------------
message("[3/6] Reading per-site TSVs and aggregating edit counts")
site_files <- list.files(sites_dir, pattern="\\.tsv(\\.gz)?$", full.names=TRUE)
if (length(site_files) == 0) stop("No TSVs found in: ", sites_dir)

edit_agg_list <- vector("list", length(site_files))

for (i in seq_along(site_files)) {
  f <- site_files[i]
  samp <- normalize_sample_id(f)
  message("  [", i, "/", length(site_files), "] ", samp)
  
  dt <- fread(f)
  req <- c("contig","position","count","feature_name")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) stop("Missing columns in ", basename(f), ": ", paste(miss, collapse=", "))
  
  dt[, contig2 := harmonize_contig_vec(contig, want_chr_prefix = gtf_has_chr)]
  pos2 <- suppressWarnings(as.integer(dt$position))
  if (isTRUE(coord_is_0based)) pos2 <- pos2 + 1L
  dt[, position2 := pos2]
  
  # drop bad gene names too
  dt[, gene := as.character(feature_name)]
  dt <- dt[
    !is.na(contig2) & !is.na(position2) & position2 >= 1L &
      !is.na(count) &
      !is.na(gene) & nzchar(gene) & gene != "."
  ]
  if (nrow(dt) == 0) next
  
  have_region_col <- "region" %in% names(dt)
  if (isTRUE(use_tsv_region_if_present) && have_region_col) {
    dt[, region2 := normalize_region_label(region)]
    dt[!region2 %in% c("UTR3","CDS","INTRON"), region2 := "OTHER"]
  } else {
    # fallback overlap labeling
    gr <- GRanges(seqnames = dt$contig2, ranges = IRanges(dt$position2, dt$position2))
    hit_utr3   <- overlaps_any_safe(gr, utr3_red,   ignore.strand=TRUE)
    hit_cds    <- overlaps_any_safe(gr, cds_red,    ignore.strand=TRUE)
    hit_intron <- overlaps_any_safe(gr, intron_red, ignore.strand=TRUE)
    dt[, region2 := pick_region_priority(hit_utr3, hit_cds, hit_intron)]
  }
  
  dt[, sample := samp]
  
  edit_agg_list[[i]] <- dt[, .(edit_counts = sum(as.numeric(count), na.rm=TRUE)),
                           by=.(sample, gene, region=region2)]
}

edit_counts <- rbindlist(edit_agg_list, use.names=TRUE, fill=TRUE)
fwrite(edit_counts, file.path(out_dir, "edit_counts_by_gene_region_sample.csv"))

# ---------------------------
# [4/6] Compute EPKM
# ---------------------------
message("[4/6] Computing EPKM")

mapped_long <- melt(mapped_dt, id.vars="sample", variable.name="region", value.name="mapped_reads")
mapped_long <- mapped_long[!is.na(mapped_reads)]
mapped_long[, region := toupper(region)]

epkm <- merge(edit_counts, gene_len, by=c("gene","region"), all.x=TRUE)
epkm <- merge(epkm, mapped_long, by=c("sample","region"), all.x=TRUE)

epkm <- epkm[region %in% c("UTR3","CDS","INTRON")]
epkm <- epkm[!is.na(region_bp) & region_bp > 0 & !is.na(mapped_reads) & mapped_reads > 0]
epkm <- epkm[mapped_reads >= min_mapped_reads_region]

epkm[, EPKM := (edit_counts * 1e9) / (mapped_reads * region_bp)]
fwrite(epkm, file.path(out_dir, "EPKM_long_gene_region_sample.csv"))

# ---------------------------
# [5/6] Wide matrices per region
# ---------------------------
message("[5/6] Writing wide EPKM matrices per region")
for (reg in sort(unique(epkm$region))) {
  w <- dcast(epkm[region == reg], gene ~ sample, value.var="EPKM", fill=0)
  fwrite(w, file.path(out_dir, paste0("EPKM_", reg, "_wide.csv")))
}

# ---------------------------
# [6/6] Done + sanity checks
# ---------------------------
message("[6/6] Done.")
message("Output folder: ", out_dir)

edit_counts2 <- fread(file.path(out_dir, "edit_counts_by_gene_region_sample.csv"))
mapped_dt2   <- fread(file.path(out_dir, "mapped_reads_per_region_per_sample.csv"))
gene_len2    <- fread(file.path(out_dir, "gene_region_lengths_bp.csv"))

message("Sanity:")
message("  n_edit_rows: ", nrow(edit_counts2))
message("  n_gene_len_rows: ", nrow(gene_len2))
message("  n_mapped_samples: ", nrow(mapped_dt2))
message("  frac genes in gene_len: ", mean(unique(edit_counts2$gene) %in% gene_len2$gene))
message("  samples in both (intersection): ", length(intersect(unique(edit_counts2$sample), mapped_dt2$sample)))

miss1 <- setdiff(unique(edit_counts2$sample), mapped_dt2$sample)
miss2 <- setdiff(mapped_dt2$sample, unique(edit_counts2$sample))
if (length(miss1) > 0) { message("  WARNING: samples in TSVs but not BAM denominators:"); print(miss1) }
if (length(miss2) > 0) { message("  NOTE: samples in BAM denominators but no TSV edits observed:"); print(miss2) }

