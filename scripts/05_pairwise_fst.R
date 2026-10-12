#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 05_pairwise_fst.R ONE_SNP_PER_LOCUS.vcf.gz SAMPLE_MAP.tsv OUTPUT_DIR")
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
source(file.path(script_dir, "popgen_helpers.R"))

required_packages(c("StAMPP", "dartR", "adegenet", "vcfR", "ggplot2"))
output_dir <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
nboot <- as.integer(Sys.getenv("NBOOT", "1000"))
ncores_default <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
ncores <- as.integer(Sys.getenv("NCORES", as.character(ncores_default)))
if (!is.finite(ncores) || ncores < 1) ncores <- 1L
seed <- as.integer(Sys.getenv("BOOTSTRAP_SEED", "20261012"))

loaded <- load_dosage(args[1], args[2])
population_levels <- unique(as.character(loaded$map$pop_gen_population))

message(
  "Calculating pairwise Weir-Cockerham FST with ", ncores,
  " cores and ", nboot, " locus bootstraps..."
)
set.seed(seed)
result <- dartR::gl.fst.pop(
  loaded$gl, nboots = nboot, percent = 95,
  nclusters = ncores, verbose = 2
)

if (nboot < 2) {
  fst <- as.matrix(result)
  lower <- upper <- matrix(
    NA_real_, nrow(fst), ncol(fst), dimnames = dimnames(fst)
  )
  boot <- NULL
} else {
  fst <- result$Fsts
  boot <- as.data.frame(result$Bootstraps, stringsAsFactors = FALSE)
  lower <- upper <- matrix(
    NA_real_, nrow(fst), ncol(fst), dimnames = dimnames(fst)
  )
  for (i in seq_len(nrow(boot))) {
    p1 <- as.character(boot$Population1[i])
    p2 <- as.character(boot$Population2[i])
    lo <- as.numeric(boot[["Lower bound CI limit"]][i])
    hi <- as.numeric(boot[["Upper bound CI limit"]][i])
    lower[p1, p2] <- lower[p2, p1] <- lo
    upper[p1, p2] <- upper[p2, p1] <- hi
  }
}

# StAMPP returns the informative triangle. Mirror it for tables and plotting.
fst[upper.tri(fst)] <- t(fst)[upper.tri(fst)]
diag(fst) <- 0
fst <- fst[population_levels, population_levels, drop = FALSE]
lower <- lower[population_levels, population_levels, drop = FALSE]
upper <- upper[population_levels, population_levels, drop = FALSE]

write.table(
  data.frame(population = rownames(fst), fst, check.names = FALSE),
  file.path(output_dir, "pairwise_wc_fst_matrix.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  data.frame(population = rownames(lower), lower, check.names = FALSE),
  file.path(output_dir, "pairwise_wc_fst_lower95_matrix.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  data.frame(population = rownames(upper), upper, check.names = FALSE),
  file.path(output_dir, "pairwise_wc_fst_upper95_matrix.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

fst_long <- matrix_to_long(fst, "Fst")
lower_long <- matrix_to_long(lower, "lower_95")
upper_long <- matrix_to_long(upper, "upper_95")
for (object_name in c("fst_long", "lower_long", "upper_long")) {
  object <- get(object_name)
  object$population_1 <- as.character(object$population_1)
  object$population_2 <- as.character(object$population_2)
  assign(object_name, object)
}
long <- merge(fst_long, lower_long, by = c("population_1", "population_2"), all = TRUE)
long <- merge(long, upper_long, by = c("population_1", "population_2"), all = TRUE)
long <- long[long$population_1 < long$population_2, , drop = FALSE]
long$nboot <- nboot
classes <- unique(loaded$map[c("pop_gen_population", "pop_gen")])
class_lookup <- setNames(classes$pop_gen, classes$pop_gen_population)
long$class_1 <- unname(class_lookup[long$population_1])
long$class_2 <- unname(class_lookup[long$population_2])
long$comparison_class <- ifelse(
  long$class_1 == long$class_2,
  paste(long$class_1, long$class_2, sep = "-"),
  "admixed-uniform"
)
write.table(
  long, file.path(output_dir, "pairwise_wc_fst_long.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

heat <- matrix_to_long(fst, "Fst")
heat$population_1 <- factor(heat$population_1, levels = population_levels)
heat$population_2 <- factor(heat$population_2, levels = rev(population_levels))
p <- ggplot2::ggplot(heat, ggplot2::aes(population_1, population_2, fill = Fst)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(option = "C", na.value = "white") +
  ggplot2::coord_fixed() +
  ggplot2::labs(x = NULL, y = NULL, fill = expression(F[ST])) +
  ggplot2::theme_minimal(base_size = 7) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top"
  )
ggplot2::ggsave(
  file.path(output_dir, "pairwise_wc_fst_heatmap.pdf"),
  p, width = 13, height = 12, limitsize = FALSE
)

saveRDS(
  list(fst = fst, bootstrap = boot, stampp_result = result, long = long),
  file.path(output_dir, "pairwise_wc_fst_objects.rds")
)
message("Pairwise FST outputs written to ", output_dir)

