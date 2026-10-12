#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 01_validate_metadata.R METADATA.csv INPUT.vcf.gz OUTPUT_DIR")
}

metadata_file <- args[1]
vcf_file <- args[2]
output_dir <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(
  metadata_file, check.names = FALSE, stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)

required <- c(
  "sampleID", "phyloID", "s384", "oldfilename", "pop_gen",
  "pop_gen_population", "spID", "geoID", "Latitude", "Longitude", "Plate"
)
missing_columns <- setdiff(required, names(metadata))
if (length(missing_columns)) {
  stop("Metadata columns not found: ", paste(missing_columns, collapse = ", "))
}
if (anyNA(metadata$sampleID) || any(!nzchar(trimws(metadata$sampleID)))) {
  stop("Every metadata row must have a non-empty sampleID")
}
if (anyDuplicated(metadata$sampleID)) {
  stop("Duplicated canonical sampleID values: ",
       paste(unique(metadata$sampleID[duplicated(metadata$sampleID)]), collapse = ", "))
}

read_vcf_ids <- function(path) {
  bcftools <- Sys.which("bcftools")
  if (nzchar(bcftools)) {
    ids <- system2(bcftools, c("query", "-l", path), stdout = TRUE, stderr = TRUE)
    status <- attr(ids, "status")
    if (is.null(status) || status == 0) return(ids)
  }
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) stop("VCF has no #CHROM header: ", path)
    if (startsWith(line, "#CHROM")) return(strsplit(line, "\t", fixed = TRUE)[[1]][-(1:9)])
  }
}

vcf_ids <- read_vcf_ids(vcf_file)
if (!length(vcf_ids)) stop("No sample identifiers found in VCF")
if (anyDuplicated(vcf_ids)) stop("VCF contains duplicated sample identifiers")

# The raw files may contain current sampleID values, s384 names, or historical
# oldfilename values. Build one explicit alias-to-canonical mapping. "REMOVED"
# is a status label in oldfilename, not a biological sample identifier.
alias_columns <- c("sampleID", "s384", "oldfilename")
alias_rows <- do.call(rbind, lapply(alias_columns, function(column) {
  value <- trimws(as.character(metadata[[column]]))
  keep <- !is.na(value) & nzchar(value) & toupper(value) != "REMOVED"
  data.frame(
    alias_type = column,
    raw_id = value[keep],
    sampleID = metadata$sampleID[keep],
    stringsAsFactors = FALSE
  )
}))
alias_rows <- unique(alias_rows)

alias_targets <- split(alias_rows$sampleID, alias_rows$raw_id)
ambiguous <- names(alias_targets)[vapply(alias_targets, function(x) length(unique(x)) > 1, logical(1))]
if (length(ambiguous)) {
  stop("Aliases map to multiple sampleID values: ", paste(ambiguous, collapse = ", "))
}
alias_lookup <- vapply(alias_targets, function(x) unique(x)[1], character(1))

unmatched_vcf <- setdiff(vcf_ids, names(alias_lookup))
if (length(unmatched_vcf)) {
  stop("VCF identifiers absent from all metadata aliases: ",
       paste(unmatched_vcf, collapse = ", "))
}
canonical_ids <- unname(alias_lookup[vcf_ids])
if (anyDuplicated(canonical_ids)) {
  stop("Multiple VCF identifiers resolve to the same canonical sampleID")
}

metadata_vcf <- metadata[match(canonical_ids, metadata$sampleID), , drop = FALSE]
if (anyNA(metadata_vcf$sampleID)) stop("Internal metadata matching failure")

valid_class <- metadata_vcf$pop_gen %in% c("uniform", "admixed")
has_population <- !is.na(metadata_vcf$pop_gen_population) &
  nzchar(trimws(metadata_vcf$pop_gen_population))
if (any(xor(valid_class, has_population))) {
  bad <- metadata_vcf$sampleID[xor(valid_class, has_population)]
  stop("pop_gen and pop_gen_population disagree for: ", paste(bad, collapse = ", "))
}
assigned <- valid_class & has_population

crosswalk <- data.frame(
  vcf_order = seq_along(vcf_ids),
  raw_id = vcf_ids,
  sampleID = metadata_vcf$sampleID,
  phyloID = metadata_vcf$phyloID,
  pop_gen = metadata_vcf$pop_gen,
  pop_gen_population = metadata_vcf$pop_gen_population,
  spID = metadata_vcf$spID,
  geoID = metadata_vcf$geoID,
  Latitude = metadata_vcf$Latitude,
  Longitude = metadata_vcf$Longitude,
  Plate = metadata_vcf$Plate,
  selected = assigned,
  stringsAsFactors = FALSE
)

selected_map <- crosswalk[assigned, , drop = FALSE]
population_summary <- aggregate(
  selected_map$sampleID,
  by = list(
    pop_gen = selected_map$pop_gen,
    pop_gen_population = selected_map$pop_gen_population
  ),
  FUN = length
)
names(population_summary)[3] <- "n_samples"
population_summary <- population_summary[
  order(population_summary$pop_gen, population_summary$pop_gen_population),,
  drop = FALSE
]

if (any(population_summary$n_samples < 3 | population_summary$n_samples > 5)) {
  bad <- population_summary[
    population_summary$n_samples < 3 | population_summary$n_samples > 5,, drop = FALSE
  ]
  stop("Selected populations must contain 3-5 samples: ",
       paste(paste0(bad$pop_gen_population, "=", bad$n_samples), collapse = ", "))
}

write.table(
  crosswalk, file.path(output_dir, "sample_to_population.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  alias_rows[order(alias_rows$sampleID, alias_rows$alias_type), ],
  file.path(output_dir, "identifier_aliases.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE, na = ""
)
writeLines(selected_map$raw_id, file.path(output_dir, "assigned_raw_ids.txt"))
writeLines(selected_map$sampleID, file.path(output_dir, "assigned_sampleIDs.txt"))
write.table(
  population_summary, file.path(output_dir, "population_summary.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE
)

dartr_metadata <- data.frame(
  id = selected_map$raw_id,
  pop = selected_map$pop_gen_population,
  lat = selected_map$Latitude,
  lon = selected_map$Longitude,
  sampleID = selected_map$sampleID,
  phyloID = selected_map$phyloID,
  pop_gen = selected_map$pop_gen,
  stringsAsFactors = FALSE
)
write.csv(dartr_metadata, file.path(output_dir, "dartR_metadata.csv"), row.names = FALSE, na = "")

expected_samples <- as.integer(Sys.getenv("EXPECTED_ASSIGNED_SAMPLES", "215"))
expected_populations <- as.integer(Sys.getenv("EXPECTED_POPULATIONS", "56"))
if (nrow(selected_map) != expected_samples) {
  warning("Selected sample count is ", nrow(selected_map), "; expected ", expected_samples)
}
if (nrow(population_summary) != expected_populations) {
  warning("Population count is ", nrow(population_summary), "; expected ", expected_populations)
}

report <- c(
  paste("Metadata rows:", nrow(metadata)),
  paste("VCF samples:", length(vcf_ids)),
  paste("Selected samples:", nrow(selected_map)),
  paste("Uniform samples:", sum(selected_map$pop_gen == "uniform")),
  paste("Admixed samples:", sum(selected_map$pop_gen == "admixed")),
  paste("Selected populations:", nrow(population_summary)),
  paste("Uniform populations:", sum(population_summary$pop_gen == "uniform")),
  paste("Admixed populations:", sum(population_summary$pop_gen == "admixed")),
  paste("VCF IDs resolved through sampleID:", sum(vcf_ids %in% metadata$sampleID)),
  paste("VCF IDs resolved through aliases:", sum(!(vcf_ids %in% metadata$sampleID)))
)
writeLines(report, file.path(output_dir, "validation_report.txt"))
message(paste(report, collapse = "\n"))

