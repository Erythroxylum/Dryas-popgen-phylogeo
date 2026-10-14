#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste(
    "Usage: 11_summarize_refugia_metrics.R",
    "POPULATION_DIVERSITY_DIR WITHIN_LINEAGE_DIR OUTPUT_DIR"
  ))
}

diversity_dir <- args[1]
lineage_dir <- args[2]
output_dir <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_required <- function(path, columns) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  absent <- setdiff(columns, names(x))
  if (length(absent)) {
    stop("Columns absent from ", path, ": ", paste(absent, collapse = ", "))
  }
  x
}

snp <- read_required(
  file.path(diversity_dir, "population_snp_diversity.tsv"),
  c("population", "pop_gen", "n_samples", "Ho", "Hs")
)
snp_ci <- read_required(
  file.path(diversity_dir, "population_snp_diversity_bootstrap_ci.tsv"),
  c("population", "metric", "estimate", "lower_95", "upper_95", "nboot")
)
pi <- read_required(
  file.path(diversity_dir, "filtered_sequence_population_pi.tsv"),
  c("population", "pop_gen", "n_loci", "callable_sites", "pi")
)
pi_ci <- read_required(
  file.path(diversity_dir, "filtered_sequence_pi_bootstrap_ci.tsv"),
  c("population", "pi", "lower_95", "upper_95", "nboot")
)
private <- read_required(
  file.path(lineage_dir, "within_lineage_private_alleles_by_population.tsv"),
  c(
    "population", "species_lineage", "n_samples",
    "n_uniform_populations_in_lineage",
    "within_lineage_private_allelic_richness", "status"
  )
)
private_ci <- read_required(
  file.path(lineage_dir, "within_lineage_private_alleles_bootstrap_ci.tsv"),
  c("population", "estimate", "lower_95", "upper_95", "nboot")
)

snp <- snp[snp$pop_gen == "uniform", , drop = FALSE]
pi <- pi[pi$pop_gen == "uniform", , drop = FALSE]

metric_ci <- function(metric, prefix) {
  z <- snp_ci[snp_ci$metric == metric, c(
    "population", "lower_95", "upper_95", "nboot"
  )]
  names(z)[-1] <- paste0(prefix, c("_lower_95", "_upper_95", "_nboot"))
  z
}

pi_ci_small <- pi_ci[c("population", "lower_95", "upper_95", "nboot")]
names(pi_ci_small)[-1] <- c("pi_lower_95", "pi_upper_95", "pi_nboot")
private_ci_small <- private_ci[c("population", "lower_95", "upper_95", "nboot")]
names(private_ci_small)[-1] <- c(
  "private_lower_95", "private_upper_95", "private_nboot"
)

population <- merge(
  private,
  snp[c("population", "Ho", "Hs")],
  by = "population", all.x = TRUE, sort = FALSE
)
population <- merge(
  population,
  pi[c("population", "n_loci", "callable_sites", "pi")],
  by = "population", all.x = TRUE, sort = FALSE
)
population <- merge(population, metric_ci("Ho", "Ho"), by = "population", all.x = TRUE, sort = FALSE)
population <- merge(population, metric_ci("Hs", "He"), by = "population", all.x = TRUE, sort = FALSE)
population <- merge(population, pi_ci_small, by = "population", all.x = TRUE, sort = FALSE)
population <- merge(population, private_ci_small, by = "population", all.x = TRUE, sort = FALSE)
names(population)[names(population) == "Hs"] <- "He"

population <- population[order(
  population$species_lineage, population$population
), , drop = FALSE]

if (anyNA(population$Ho) || anyNA(population$He) || anyNA(population$pi)) {
  stop("At least one uniform population is missing Ho, He, or pi")
}

summarize_metric <- function(x, name) {
  good <- is.finite(x)
  values <- x[good]
  if (!length(values)) {
    return(setNames(as.list(c(0, NA, NA, NA, NA)), paste0(
      name, c("_n", "_mean", "_median", "_minimum", "_maximum")
    )))
  }
  setNames(
    as.list(c(length(values), mean(values), median(values), min(values), max(values))),
    paste0(name, c("_n", "_mean", "_median", "_minimum", "_maximum"))
  )
}

lineage_rows <- lapply(split(population, population$species_lineage), function(x) {
  fields <- c(
    list(
      species_lineage = x$species_lineage[1],
      n_uniform_populations = nrow(x),
      n_samples = sum(x$n_samples)
    ),
    summarize_metric(x$pi, "pi"),
    summarize_metric(x$Ho, "Ho"),
    summarize_metric(x$He, "He"),
    summarize_metric(
      x$within_lineage_private_allelic_richness,
      "within_lineage_private_allelic_richness"
    )
  )
  as.data.frame(fields, stringsAsFactors = FALSE)
})
lineage <- do.call(rbind, lineage_rows)
rownames(lineage) <- NULL
lineage <- lineage[order(lineage$species_lineage), , drop = FALSE]

write.table(
  population,
  file.path(output_dir, "uniform_population_refugia_metrics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  lineage,
  file.path(output_dir, "uniform_lineage_diversity_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

message(
  "Combined metrics: ",
  file.path(output_dir, "uniform_population_refugia_metrics.tsv")
)
message(
  "Lineage summary: ",
  file.path(output_dir, "uniform_lineage_diversity_summary.tsv")
)

