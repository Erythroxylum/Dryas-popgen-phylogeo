#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 03_population_diversity.R ONE_SNP_PER_LOCUS.vcf.gz SAMPLE_MAP.tsv OUTPUT_DIR")
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
source(file.path(script_dir, "popgen_helpers.R"))

required_packages(c("hierfstat", "dartR", "adegenet", "vcfR", "ggplot2"))
vcf_file <- args[1]
map_file <- args[2]
output_dir <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
nboot <- as.integer(Sys.getenv("NBOOT", "1000"))
seed <- as.integer(Sys.getenv("BOOTSTRAP_SEED", "20261012"))
min_gene_copies <- as.integer(Sys.getenv("RAREFACTION_GENE_COPIES", "6"))
min_outside_call_rate <- as.numeric(Sys.getenv("PRIVATE_MIN_OUTSIDE_CALL_RATE", "0.80"))

loaded <- load_dosage(vcf_file, map_file)
dos <- loaded$dos
map <- loaded$map
hf_object <- make_hierfstat(dos, map$pop_gen_population)
hf <- hf_object$data
population_levels <- hf_object$levels

message("Calculating hierfstat basic statistics across ", ncol(dos),
        " independent ipyrad-locus representatives...")
basic <- hierfstat::basic.stats(hf, diploid = TRUE)
Ho <- rename_population_matrix(basic$Ho, population_levels)
Hs <- rename_population_matrix(basic$Hs, population_levels)
Fis_locus <- rename_population_matrix(basic$Fis, population_levels)

allelic <- hierfstat::allelic.richness(
  hf, min.n = min_gene_copies, diploid = TRUE
)
Ar <- rename_population_matrix(allelic$Ar, population_levels)

population_metadata <- unique(map[c("pop_gen_population", "pop_gen")])
population_metadata <- population_metadata[
  match(population_levels, population_metadata$pop_gen_population),,
  drop = FALSE
]
sample_counts <- table(factor(map$pop_gen_population, levels = population_levels))
called_per_locus <- basic$n.ind.samp
called_per_locus <- rename_population_matrix(called_per_locus, population_levels)

point <- data.frame(
  population = population_levels,
  pop_gen = population_metadata$pop_gen,
  n_samples = as.integer(sample_counts),
  loci_with_data = colSums(called_per_locus > 0, na.rm = TRUE),
  mean_called_individuals = mean_na(called_per_locus),
  Ho = mean_na(Ho),
  Hs = mean_na(Hs),
  Fis_ratio = 1 - mean_na(Ho) / mean_na(Hs),
  mean_locus_Fis = mean_na(Fis_locus),
  allelic_richness = mean_na(Ar),
  rarefied_gene_copies = allelic$min.all,
  stringsAsFactors = FALSE
)

# Resample the one-SNP-per-ipyrad-locus rows. rmultinom produces bootstrap
# multiplicities without constructing NBOOT copies of each statistic matrix.
set.seed(seed)
n_loci <- nrow(Ho)
bootstrap_counts <- rmultinom(nboot, size = n_loci, prob = rep(1 / n_loci, n_loci))

bootstrap_means <- function(x) {
  observed <- !is.na(x)
  values <- x
  values[!observed] <- 0
  numerator <- t(values) %*% bootstrap_counts
  denominator <- t(observed * 1) %*% bootstrap_counts
  result <- numerator / denominator
  result[denominator == 0] <- NA_real_
  result
}

boot_Ho <- bootstrap_means(Ho)
boot_Hs <- bootstrap_means(Hs)
boot_Ar <- bootstrap_means(Ar)
boot_Fis <- 1 - boot_Ho / boot_Hs

ci_rows <- lapply(seq_along(population_levels), function(i) {
  metrics <- list(Ho = boot_Ho[i, ], Hs = boot_Hs[i, ],
                  Fis_ratio = boot_Fis[i, ], allelic_richness = boot_Ar[i, ])
  do.call(rbind, lapply(names(metrics), function(metric) {
    values <- metrics[[metric]]
    values <- values[is.finite(values)]
    interval <- if (length(values)) quantile(values, c(0.025, 0.975), names = FALSE) else c(NA, NA)
    data.frame(
      population = population_levels[i], metric = metric,
      estimate = point[[metric]][i], lower_95 = interval[1], upper_95 = interval[2],
      nboot = nboot, stringsAsFactors = FALSE
    )
  }))
})
confidence_intervals <- do.call(rbind, ci_rows)

# Calculate private alleles explicitly. A candidate allele must be observed in
# the focal population, absent from all called chromosomes outside it, have at
# least the requested outside-sample call rate, and have enough focal gene
# copies for rarefaction. The hypergeometric probability gives its expected
# representation after rarefying the focal population to min_gene_copies.
rarefaction_probability <- function(total_copies, allele_copies, draw) {
  if (total_copies < draw || allele_copies < 1) return(NA_real_)
  nonallele <- total_copies - allele_copies
  miss <- if (nonallele < draw) 0 else exp(
    lchoose(nonallele, draw) - lchoose(total_copies, draw)
  )
  1 - miss
}

private_contribution <- matrix(
  NA_real_, nrow = ncol(dos), ncol = length(population_levels),
  dimnames = list(colnames(dos), population_levels)
)
private_rows <- list()
private_row_index <- 0L

for (p in seq_along(population_levels)) {
  focal <- map$pop_gen_population == population_levels[p]
  outside <- !focal
  outside_n <- sum(outside)
  for (locus in seq_len(ncol(dos))) {
    focal_dos <- dos[focal, locus]
    outside_dos <- dos[outside, locus]
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
    counts <- c(REF = focal_ref, ALT = focal_alt)
    outside_counts <- c(REF = outside_ref, ALT = outside_alt)

    for (allele in names(counts)) {
      if (counts[allele] > 0 && outside_counts[allele] == 0) {
        expected <- rarefaction_probability(
          focal_copies, counts[allele], min_gene_copies
        )
        private_contribution[locus, p] <-
          private_contribution[locus, p] + expected
        private_row_index <- private_row_index + 1L
        private_rows[[private_row_index]] <- data.frame(
          population = population_levels[p], locus = colnames(dos)[locus],
          allele = allele, focal_allele_copies = counts[allele],
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
    population = character(), locus = character(), allele = character(),
    focal_allele_copies = integer(), focal_called_gene_copies = integer(),
    outside_called_samples = integer(), outside_call_rate = numeric(),
    rarefied_detection_probability = numeric()
  )
}
private_summary <- data.frame(
  population = population_levels,
  eligible_loci = colSums(!is.na(private_contribution)),
  observed_private_alleles = vapply(
    population_levels,
    function(z) sum(private_detail$population == z), integer(1)
  ),
  rarefied_private_alleles = colSums(private_contribution, na.rm = TRUE),
  private_allelic_richness = mean_na(private_contribution),
  rarefied_gene_copies = min_gene_copies,
  minimum_outside_call_rate = min_outside_call_rate,
  stringsAsFactors = FALSE
)
point$private_allelic_richness <- private_summary$private_allelic_richness[
  match(point$population, private_summary$population)
]

boot_private <- bootstrap_means(private_contribution)
private_ci <- do.call(rbind, lapply(seq_along(population_levels), function(i) {
  values <- boot_private[i, ]
  values <- values[is.finite(values)]
  interval <- if (length(values)) {
    quantile(values, c(0.025, 0.975), names = FALSE)
  } else c(NA_real_, NA_real_)
  data.frame(
    population = population_levels[i], metric = "private_allelic_richness",
    estimate = point$private_allelic_richness[i],
    lower_95 = interval[1], upper_95 = interval[2], nboot = nboot,
    stringsAsFactors = FALSE
  )
}))
confidence_intervals <- rbind(confidence_intervals, private_ci)

write.table(
  point, file.path(output_dir, "population_snp_diversity.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  confidence_intervals,
  file.path(output_dir, "population_snp_diversity_bootstrap_ci.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  private_summary, file.path(output_dir, "private_alleles_by_population.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  private_detail, file.path(output_dir, "private_alleles_by_locus.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

plot_data <- merge(
  confidence_intervals,
  point[c("population", "pop_gen")],
  by = "population", all.x = TRUE, sort = FALSE
)
metric_labels <- c(
  Ho = "Observed heterozygosity",
  Hs = "Gene diversity",
  Fis_ratio = "FIS",
  allelic_richness = "Rarefied allelic richness",
  private_allelic_richness = "Private allelic richness"
)
plot_data$metric <- factor(plot_data$metric, levels = names(metric_labels))
p <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = pop_gen, y = estimate, color = pop_gen)
) +
  ggplot2::geom_pointrange(
    ggplot2::aes(ymin = lower_95, ymax = upper_95),
    position = ggplot2::position_jitter(width = 0.12, height = 0), alpha = 0.75
  ) +
  ggplot2::facet_wrap(
    ~metric, scales = "free_y", ncol = 2,
    labeller = ggplot2::as_labeller(metric_labels)
  ) +
  ggplot2::scale_color_manual(values = c(uniform = "#2166AC", admixed = "#B2182B")) +
  ggplot2::labs(x = NULL, y = NULL, color = "Population class") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "top")
ggplot2::ggsave(
  file.path(output_dir, "population_snp_diversity_diagnostic.pdf"),
  p, width = 8.2, height = 6.5
)

saveRDS(
  list(basic_stats = basic, allelic_richness = allelic,
       private_alleles = private_summary,
       private_alleles_by_locus = private_detail,
       population_table = point, bootstrap_ci = confidence_intervals),
  file.path(output_dir, "population_snp_diversity_objects.rds")
)
message("Population diversity outputs written to ", output_dir)

