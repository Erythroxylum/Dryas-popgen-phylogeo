#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 05_pairwise_fst.R ONE_SNP_PER_LOCUS.vcf.gz SAMPLE_MAP.tsv OUTPUT_DIR")
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
source(file.path(script_dir, "popgen_helpers.R"))

required_packages(c("hierfstat", "dartR", "adegenet", "vcfR", "ggplot2"))
output_dir <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
nboot <- as.integer(Sys.getenv("NBOOT", "1000"))
seed <- as.integer(Sys.getenv("BOOTSTRAP_SEED", "20261012"))

loaded <- load_dosage(args[1], args[2])
hf_object <- make_hierfstat(loaded$dos, loaded$map$pop_gen_population)
hf <- hf_object$data
population_levels <- hf_object$levels

message("Calculating pairwise Weir-Cockerham FST...")
fst <- hierfstat::pairwise.WCfst(hf, diploid = TRUE)
fst <- rename_population_matrix(fst, population_levels)

message("Bootstrapping ", nboot, " times over independent ipyrad loci...")
set.seed(seed)
boot <- hierfstat::boot.ppfst(
  dat = hf, nboot = nboot, quant = c(0.025, 0.975), diploid = TRUE
)
lower <- rename_population_matrix(boot$ll, population_levels)
upper <- rename_population_matrix(boot$ul, population_levels)

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
  list(fst = fst, bootstrap = boot, long = long),
  file.path(output_dir, "pairwise_wc_fst_objects.rds")
)
message("Pairwise FST outputs written to ", output_dir)

