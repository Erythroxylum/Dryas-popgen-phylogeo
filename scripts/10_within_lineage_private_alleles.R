#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste(
    "Usage: 10_within_lineage_private_alleles.R",
    "ONE_SNP_PER_LOCUS.vcf.gz SAMPLE_MAP.tsv OUTPUT_DIR"
  ))
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
source(file.path(script_dir, "popgen_helpers.R"))

required_packages(c("dartR", "adegenet", "vcfR"))
vcf_file <- args[1]
map_file <- args[2]
output_dir <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

nboot <- as.integer(Sys.getenv("NBOOT", "1000"))
seed <- as.integer(Sys.getenv("BOOTSTRAP_SEED", "20261012"))
min_gene_copies <- as.integer(Sys.getenv("RAREFACTION_GENE_COPIES", "6"))
min_outside_call_rate <- as.numeric(
  Sys.getenv("PRIVATE_MIN_OUTSIDE_CALL_RATE", "0.80")
)
expected_lineages <- as.integer(Sys.getenv("EXPECTED_UNIFORM_LINEAGES", "11"))

loaded <- load_dosage(vcf_file, map_file)
dos <- loaded$dos
map <- loaded$map

# The revised spID column is the evolutionary-lineage assignment. Restrict the
# comparison to previously delimited uniform geographic populations so that
# introgressed populations are not forced into a single lineage.
uniform <- map$pop_gen == "uniform"
if (!any(uniform)) stop("No samples have pop_gen == 'uniform'")
if (any(is.na(map$spID[uniform]) | !nzchar(trimws(map$spID[uniform])))) {
  stop("Every uniform sample must have a non-empty lineage assignment in spID")
}

uniform_map <- map[uniform, , drop = FALSE]
uniform_dos <- dos[uniform, , drop = FALSE]
uniform_map$species_lineage <- trimws(as.character(uniform_map$spID))
uniform_map$population <- trimws(as.character(uniform_map$pop_gen_population))

population_lineages <- split(
  uniform_map$species_lineage, uniform_map$population
)
mixed <- names(population_lineages)[vapply(
  population_lineages, function(x) length(unique(x)) != 1, logical(1)
)]
if (length(mixed)) {
  stop(
    "Uniform populations assigned to more than one spID lineage: ",
    paste(mixed, collapse = ", ")
  )
}

population_map <- unique(uniform_map[c("population", "species_lineage")])
population_map <- population_map[
  order(population_map$species_lineage, population_map$population),,
  drop = FALSE
]
sample_counts <- table(uniform_map$population)
population_map$n_samples <- as.integer(sample_counts[population_map$population])
lineage_counts <- table(population_map$species_lineage)
population_map$n_uniform_populations_in_lineage <- as.integer(
  lineage_counts[population_map$species_lineage]
)

if (length(unique(population_map$species_lineage)) != expected_lineages) {
  warning(
    "Found ", length(unique(population_map$species_lineage)),
    " lineages among uniform populations; expected ", expected_lineages,
    ". Check the revised spID assignments."
  )
}

rarefaction_probability <- function(total_copies, allele_copies, draw) {
  if (total_copies < draw || allele_copies < 1) return(NA_real_)
  nonallele <- total_copies - allele_copies
  miss <- if (nonallele < draw) 0 else exp(
    lchoose(nonallele, draw) - lchoose(total_copies, draw)
  )
  1 - miss
}

populations <- population_map$population
private_contribution <- matrix(
  NA_real_, nrow = ncol(uniform_dos), ncol = length(populations),
  dimnames = list(colnames(uniform_dos), populations)
)
private_rows <- list()
private_row_index <- 0L

message(
  "Calculating private alleles within ",
  length(unique(population_map$species_lineage)),
  " lineages and ", length(populations), " uniform populations..."
)

for (p in seq_along(populations)) {
  focal_population <- populations[p]
  lineage <- population_map$species_lineage[p]
  lineage_populations <- population_map$population[
    population_map$species_lineage == lineage
  ]

  # A single sampled population provides no within-lineage comparison. Its
  # private-allele estimate is undefined rather than zero.
  if (length(lineage_populations) < 2) next

  focal <- uniform_map$population == focal_population
  outside <- uniform_map$species_lineage == lineage & !focal
  outside_n <- sum(outside)
  if (!outside_n) next

  for (locus in seq_len(ncol(uniform_dos))) {
    focal_dos <- uniform_dos[focal, locus]
    outside_dos <- uniform_dos[outside, locus]
    focal_called <- sum(!is.na(focal_dos))
    outside_called <- sum(!is.na(outside_dos))
    focal_copies <- 2L * focal_called
    outside_call_rate <- outside_called / outside_n

    if (focal_copies < min_gene_copies ||
        outside_call_rate < min_outside_call_rate) next

    private_contribution[locus, p] <- 0
    focal_alt <- sum(focal_dos, na.rm = TRUE)
    outside_alt <- sum(outside_dos, na.rm = TRUE)
    focal_ref <- focal_copies - focal_alt
    outside_ref <- 2L * outside_called - outside_alt
    focal_counts <- c(REF = focal_ref, ALT = focal_alt)
    outside_counts <- c(REF = outside_ref, ALT = outside_alt)

    for (allele in names(focal_counts)) {
      if (focal_counts[allele] > 0 && outside_counts[allele] == 0) {
        expected <- rarefaction_probability(
          focal_copies, focal_counts[allele], min_gene_copies
        )
        private_contribution[locus, p] <-
          private_contribution[locus, p] + expected
        private_row_index <- private_row_index + 1L
        private_rows[[private_row_index]] <- data.frame(
          species_lineage = lineage,
          population = focal_population,
          locus = colnames(uniform_dos)[locus],
          allele = allele,
          focal_allele_copies = focal_counts[allele],
          focal_called_gene_copies = focal_copies,
          outside_called_samples = outside_called,
          outside_call_rate = outside_call_rate,
          rarefied_detection_probability = expected,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

private_detail <- if (length(private_rows)) do.call(rbind, private_rows) else {
  data.frame(
    species_lineage = character(), population = character(),
    locus = character(), allele = character(),
    focal_allele_copies = integer(), focal_called_gene_copies = integer(),
    outside_called_samples = integer(), outside_call_rate = numeric(),
    rarefied_detection_probability = numeric(), stringsAsFactors = FALSE
  )
}

eligible_loci <- colSums(!is.na(private_contribution))
private_richness <- vapply(seq_along(populations), function(i) {
  x <- private_contribution[, i]
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}, numeric(1))

private_summary <- population_map
private_summary$eligible_loci <- as.integer(eligible_loci[private_summary$population])
private_summary$observed_private_alleles <- vapply(
  private_summary$population,
  function(z) sum(private_detail$population == z), integer(1)
)
private_summary$rarefied_private_alleles <- colSums(
  private_contribution[, private_summary$population, drop = FALSE], na.rm = TRUE
)
private_summary$within_lineage_private_allelic_richness <- private_richness[
  match(private_summary$population, populations)
]
private_summary$rarefied_gene_copies <- min_gene_copies
private_summary$minimum_outside_call_rate <- min_outside_call_rate
private_summary$status <- ifelse(
  private_summary$n_uniform_populations_in_lineage < 2,
  "not_estimable_single_uniform_population", "estimated"
)
not_estimable <- private_summary$status != "estimated"
private_summary$eligible_loci[not_estimable] <- NA_integer_
private_summary$observed_private_alleles[not_estimable] <- NA_integer_
private_summary$rarefied_private_alleles[not_estimable] <- NA_real_
private_summary$within_lineage_private_allelic_richness[not_estimable] <- NA_real_

# Resample independent ipyrad-locus representatives. All populations use the
# same bootstrap locus draws, retaining covariance for later comparisons.
set.seed(seed)
n_loci <- nrow(private_contribution)
bootstrap_counts <- rmultinom(
  nboot, size = n_loci, prob = rep(1 / n_loci, n_loci)
)
observed <- !is.na(private_contribution)
values <- private_contribution
values[!observed] <- 0
numerator <- t(values) %*% bootstrap_counts
denominator <- t(observed * 1) %*% bootstrap_counts
boot <- numerator / denominator
boot[denominator == 0] <- NA_real_

private_ci <- do.call(rbind, lapply(seq_along(populations), function(i) {
  values <- boot[i, ]
  values <- values[is.finite(values)]
  interval <- if (length(values)) {
    quantile(values, c(0.025, 0.975), names = FALSE)
  } else c(NA_real_, NA_real_)
  data.frame(
    population = populations[i],
    species_lineage = population_map$species_lineage[
      match(populations[i], population_map$population)
    ],
    estimate = private_richness[i],
    lower_95 = interval[1], upper_95 = interval[2], nboot = nboot,
    stringsAsFactors = FALSE
  )
}))

write.table(
  population_map,
  file.path(output_dir, "uniform_population_lineage_map.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  private_summary,
  file.path(output_dir, "within_lineage_private_alleles_by_population.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  private_ci,
  file.path(output_dir, "within_lineage_private_alleles_bootstrap_ci.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  private_detail,
  file.path(output_dir, "within_lineage_private_alleles_by_locus.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

message("Within-lineage private-allele tables: ", output_dir)

