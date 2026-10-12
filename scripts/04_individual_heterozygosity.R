#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 04_individual_heterozygosity.R ALL_SNPS.vcf.gz SAMPLE_MAP.tsv OUTPUT.tsv")
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
source(file.path(script_dir, "popgen_helpers.R"))

loaded <- load_dosage(args[1], args[2])
dos <- loaded$dos
map <- loaded$map
called <- rowSums(!is.na(dos))
heterozygous <- rowSums(dos == 1, na.rm = TRUE)
homozygous_reference <- rowSums(dos == 0, na.rm = TRUE)
homozygous_alternate <- rowSums(dos == 2, na.rm = TRUE)

out <- data.frame(
  raw_id = map$raw_id,
  sampleID = map$sampleID,
  phyloID = map$phyloID,
  pop_gen = map$pop_gen,
  population = map$pop_gen_population,
  Plate = map$Plate,
  total_snps = ncol(dos),
  called_snps = called,
  call_rate = called / ncol(dos),
  heterozygous_snps = heterozygous,
  homozygous_reference_snps = homozygous_reference,
  homozygous_alternate_snps = homozygous_alternate,
  snp_heterozygosity = ifelse(called > 0, heterozygous / called, NA_real_),
  stringsAsFactors = FALSE
)
dir.create(dirname(args[3]), recursive = TRUE, showWarnings = FALSE)
write.table(out, args[3], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Individual SNP heterozygosity: ", args[3])

