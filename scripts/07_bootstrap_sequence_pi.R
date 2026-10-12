#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 07_bootstrap_sequence_pi.R PER_LOCUS_PI.tsv.gz OUTPUT.tsv")
}
nboot <- as.integer(Sys.getenv("NBOOT", "1000"))
seed <- as.integer(Sys.getenv("BOOTSTRAP_SEED", "20261012"))
dat <- read.delim(args[1], check.names = FALSE, stringsAsFactors = FALSE)
required <- c(
  "locus", "population", "pop_gen", "different_allele_pairs",
  "callable_allele_pairs"
)
absent <- setdiff(required, names(dat))
if (length(absent)) stop("Columns not found: ", paste(absent, collapse = ", "))

set.seed(seed)
result <- do.call(rbind, lapply(split(dat, dat$population), function(x) {
  n <- nrow(x)
  counts <- rmultinom(nboot, n, rep(1 / n, n))
  numerator <- as.numeric(crossprod(x$different_allele_pairs, counts))
  denominator <- as.numeric(crossprod(x$callable_allele_pairs, counts))
  values <- numerator / denominator
  values <- values[is.finite(values)]
  interval <- quantile(values, c(0.025, 0.975), names = FALSE)
  data.frame(
    population = x$population[1], pop_gen = x$pop_gen[1], n_loci = n,
    pi = sum(x$different_allele_pairs) / sum(x$callable_allele_pairs),
    lower_95 = interval[1], upper_95 = interval[2], nboot = nboot,
    stringsAsFactors = FALSE
  )
}))
write.table(result, args[2], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Sequence nucleotide-diversity bootstrap table: ", args[2])

