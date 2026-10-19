#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 04_individual_heterozygosity.R ALL_SAMPLES.vcf.gz SAMPLE_MAP.tsv OUTPUT.tsv")
}

vcf_file <- args[1]
map_file <- args[2]
output_file <- args[3]

bcftools <- Sys.which("bcftools")
if (!nzchar(bcftools)) stop("bcftools was not found on PATH")
if (!file.exists(vcf_file)) stop("VCF not found: ", vcf_file)
if (!file.exists(map_file)) stop("Sample map not found: ", map_file)

# bcftools stats reports per-sample counts without loading the full genotype
# matrix into R. For the biallelic SNP VCF used here, observed individual
# heterozygosity is nHets / (nRefHom + nNonRefHom + nHets).
stderr_file <- tempfile(fileext = ".bcftools.stderr")
on.exit(unlink(stderr_file), add = TRUE)
stats_lines <- system2(
  bcftools,
  c("stats", "-s", "-", vcf_file),
  stdout = TRUE, stderr = stderr_file
)
status <- attr(stats_lines, "status")
if (!is.null(status) && status != 0) {
  details <- if (file.exists(stderr_file)) {
    paste(readLines(stderr_file, warn = FALSE), collapse = "\n")
  } else {
    ""
  }
  stop("bcftools stats failed", if (nzchar(details)) paste0(":\n", details))
}

psc_lines <- grep("^PSC\\t", stats_lines, value = TRUE)
if (!length(psc_lines)) stop("bcftools stats returned no PSC sample records")
fields <- strsplit(psc_lines, "\t", fixed = TRUE)
if (any(lengths(fields) < 14)) stop("Unexpected bcftools PSC record format")

value <- function(position) {
  suppressWarnings(as.numeric(vapply(fields, `[`, character(1), position)))
}
sample_stats <- data.frame(
  # PSC columns are: code, set ID, sample, nRefHom, nNonRefHom, nHets,
  # nTransitions, nTransversions, nIndels, average depth, nSingletons,
  # nHapRef, nHapAlt, and nMissing.
  raw_id = vapply(fields, `[`, character(1), 3),
  homozygous_reference_snps = value(4),
  homozygous_alternate_snps = value(5),
  heterozygous_snps = value(6),
  missing_snps = value(14),
  stringsAsFactors = FALSE
)
if (anyDuplicated(sample_stats$raw_id)) {
  stop("bcftools returned duplicated sample identifiers")
}
count_columns <- setdiff(names(sample_stats), "raw_id")
if (any(!is.finite(as.matrix(sample_stats[count_columns])))) {
  stop("Non-numeric genotype counts in bcftools PSC records")
}

map <- read.delim(map_file, check.names = FALSE, stringsAsFactors = FALSE)
needed_map <- c(
  "raw_id", "sampleID", "phyloID", "pop_gen", "pop_gen_population",
  "Plate", "selected"
)
absent <- setdiff(needed_map, names(map))
if (length(absent)) {
  stop("Sample map lacks columns: ", paste(absent, collapse = ", "))
}
if (anyDuplicated(map$raw_id)) stop("Sample map contains duplicated raw_id values")

map <- map[match(sample_stats$raw_id, map$raw_id), , drop = FALSE]
if (any(is.na(map$raw_id))) stop("VCF samples are missing from the sample map")

called <- sample_stats$homozygous_reference_snps +
  sample_stats$homozygous_alternate_snps + sample_stats$heterozygous_snps
total <- called + sample_stats$missing_snps

out <- data.frame(
  raw_id = sample_stats$raw_id,
  sampleID = map$sampleID,
  phyloID = map$phyloID,
  pop_gen = map$pop_gen,
  population = map$pop_gen_population,
  Plate = map$Plate,
  selected = map$selected,
  total_snps = total,
  called_snps = called,
  call_rate = ifelse(total > 0, called / total, NA_real_),
  heterozygous_snps = sample_stats$heterozygous_snps,
  homozygous_reference_snps = sample_stats$homozygous_reference_snps,
  homozygous_alternate_snps = sample_stats$homozygous_alternate_snps,
  snp_heterozygosity = ifelse(
    called > 0, sample_stats$heterozygous_snps / called, NA_real_
  ),
  stringsAsFactors = FALSE
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.table(
  out, output_file, sep = "\t", quote = FALSE,
  row.names = FALSE, na = ""
)
message(
  "Individual SNP heterozygosity: ", output_file,
  " (", nrow(out), " samples)"
)

