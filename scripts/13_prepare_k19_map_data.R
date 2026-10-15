#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste(
    "Usage: 13_prepare_k19_map_data.R",
    "ADMIXTURE_INPUT_DIR METADATA.csv OUTPUT_DIR"
  ))
}

suppressPackageStartupMessages(library(sf))

admixture_dir <- args[1]
metadata_file <- args[2]
output_dir <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

jitter_radius_km <- as.numeric(Sys.getenv("MAP_JITTER_RADIUS_KM", "45"))
if (!is.finite(jitter_radius_km) || jitter_radius_km < 0) {
  stop("MAP_JITTER_RADIUS_KM must be a non-negative number")
}

selected <- data.frame(
  K = 16:19,
  replicate = c(4L, 1L, 5L, 10L),
  stringsAsFactors = FALSE
)
selected$key <- paste0("K", selected$K, "_R", selected$replicate)

safe_correlation <- function(x, y) {
  sx <- sd(x, na.rm = TRUE)
  sy <- sd(y, na.rm = TRUE)
  if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) return(-Inf)
  value <- suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
  if (is.finite(value)) value else -Inf
}

align_to_previous <- function(previous, current) {
  scores <- outer(
    seq_len(ncol(previous)), seq_len(ncol(current)),
    Vectorize(function(i, j) safe_correlation(previous[, i], current[, j]))
  )
  previous_used <- rep(FALSE, ncol(previous))
  current_used <- rep(FALSE, ncol(current))
  matched <- matrix(integer(), nrow = 0, ncol = 2)

  for (step in seq_len(min(ncol(previous), ncol(current)))) {
    available <- scores
    available[previous_used, ] <- -Inf
    available[, current_used] <- -Inf
    if (!any(is.finite(available))) {
      stop("No finite component correlations remain during alignment")
    }
    chosen <- which(available == max(available), arr.ind = TRUE)[1, ]
    matched <- rbind(matched, chosen)
    previous_used[chosen[1]] <- TRUE
    current_used[chosen[2]] <- TRUE
  }
  matched_current <- matched[order(matched[, 1]), 2]
  new_order <- c(matched_current, which(!current_used))
  current[, new_order, drop = FALSE]
}

find_q_file <- function(k, replicate) {
  candidates <- c(
    file.path(admixture_dir, paste0("K", k, ".rep", replicate, ".Q")),
    file.path(
      admixture_dir, paste0("rep", replicate),
      paste0("K", k, ".rep", replicate, ".Q")
    )
  )
  found <- candidates[file.exists(candidates)]
  if (!length(found)) {
    stop("Q file not found. Checked: ", paste(candidates, collapse = ", "))
  }
  found[1]
}

sample_order_file <- file.path(admixture_dir, "sample_order.txt")
tree_order_file <- file.path(admixture_dir, "tree_plot_order.txt")
if (!file.exists(sample_order_file)) stop("Missing: ", sample_order_file)
if (!file.exists(tree_order_file)) {
  stop(
    "Missing tree_plot_order.txt. This is required to reproduce the K16 ",
    "component ordering used in the selected-run ADMIXTURE figure."
  )
}

q_sample_order <- readLines(sample_order_file)
tree_plot_order <- readLines(tree_order_file)
if (anyDuplicated(q_sample_order) || anyDuplicated(tree_plot_order)) {
  stop("An ADMIXTURE order file contains duplicate sample IDs")
}
if (!setequal(q_sample_order, tree_plot_order)) {
  stop("sample_order.txt and tree_plot_order.txt contain different samples")
}

metadata <- read.csv(
  metadata_file, check.names = FALSE, stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)
required_metadata <- c(
  "sampleID", "phyloID", "spID", "geoID", "Latitude", "Longitude",
  "pop_gen", "pop_gen_population"
)
absent <- setdiff(required_metadata, names(metadata))
if (length(absent)) {
  stop("Metadata columns missing: ", paste(absent, collapse = ", "))
}
if (anyDuplicated(metadata$sampleID)) stop("Metadata sampleID values are duplicated")

map_metadata <- metadata[match(q_sample_order, metadata$sampleID), , drop = FALSE]
if (any(is.na(map_metadata$sampleID))) {
  stop("Some ADMIXTURE samples are absent from metadata")
}
map_metadata$Latitude <- suppressWarnings(as.numeric(map_metadata$Latitude))
map_metadata$Longitude <- suppressWarnings(as.numeric(map_metadata$Longitude))
if (any(!is.finite(map_metadata$Latitude) | !is.finite(map_metadata$Longitude))) {
  bad <- map_metadata$sampleID[
    !is.finite(map_metadata$Latitude) | !is.finite(map_metadata$Longitude)
  ]
  stop("Missing coordinates for ADMIXTURE samples: ", paste(bad, collapse = ", "))
}

q_matrices <- list()
for (i in seq_len(nrow(selected))) {
  q <- as.matrix(read.table(
    find_q_file(selected$K[i], selected$replicate[i]), header = FALSE
  ))
  storage.mode(q) <- "double"
  if (nrow(q) != length(q_sample_order) || ncol(q) != selected$K[i]) {
    stop("Unexpected dimensions for ", selected$key[i])
  }
  if (any(!is.finite(q)) || any(q < 0)) stop("Invalid Q values in ", selected$key[i])
  if (any(abs(rowSums(q) - 1) > 0.002)) {
    stop("Q rows do not sum to one in ", selected$key[i])
  }
  rownames(q) <- q_sample_order
  colnames(q) <- paste0("raw", seq_len(ncol(q)))
  q_matrices[[selected$key[i]]] <- q
}

# Reproduce the component ordering used by scripts/12_plot_selected_runs_cex.R:
# K16 is ordered by the position of each component maximum in active-tree
# order, and each subsequent K is aligned to the previous selected solution.
first_key <- selected$key[1]
first_q <- q_matrices[[first_key]]
ordered_first <- first_q[match(tree_plot_order, rownames(first_q)), , drop = FALSE]
first_order <- order(apply(ordered_first, 2, which.max))
q_matrices[[first_key]] <- first_q[, first_order, drop = FALSE]
for (i in 2:nrow(selected)) {
  q_matrices[[selected$key[i]]] <- align_to_previous(
    q_matrices[[selected$key[i - 1]]], q_matrices[[selected$key[i]]]
  )
}

q19 <- q_matrices[["K19_R10"]]
q19 <- q19[match(q_sample_order, rownames(q19)), , drop = FALSE]
raw_component <- colnames(q19)
colnames(q19) <- sprintf("Q%02d", seq_len(ncol(q19)))

population_class <- ifelse(
  map_metadata$pop_gen %in% c("uniform", "admixed"),
  map_metadata$pop_gen, "not_selected"
)
map_data <- data.frame(
  sampleID = map_metadata$sampleID,
  phyloID = map_metadata$phyloID,
  spID = map_metadata$spID,
  geoID = map_metadata$geoID,
  Latitude = map_metadata$Latitude,
  Longitude = map_metadata$Longitude,
  pop_gen = population_class,
  pop_gen_population = map_metadata$pop_gen_population,
  total_minor_ancestry = 1 - apply(q19, 1, max),
  dominant_component = sprintf("Q%02d", max.col(q19, ties.method = "first")),
  q19,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Project first and apply deterministic sunflower offsets to samples sharing a
# geoID. This avoids stochastic map layouts while preserving regional centers.
points <- st_as_sf(
  map_data, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE
)
points <- st_transform(points, 3995)
xy <- st_coordinates(points)
map_data$x_3995 <- xy[, 1]
map_data$y_3995 <- xy[, 2]
map_data$x_plot <- map_data$x_3995
map_data$y_plot <- map_data$y_3995

group_id <- trimws(as.character(map_data$geoID))
missing_group <- is.na(group_id) | !nzchar(group_id)
group_id[missing_group] <- paste0("sample_", map_data$sampleID[missing_group])
golden_angle <- pi * (3 - sqrt(5))

for (idx in split(seq_len(nrow(map_data)), group_id)) {
  idx <- idx[order(map_data$sampleID[idx])]
  n <- length(idx)
  if (n < 2 || jitter_radius_km == 0) next
  center_x <- mean(map_data$x_3995[idx])
  center_y <- mean(map_data$y_3995[idx])
  rank <- seq_len(n) - 1
  radius <- jitter_radius_km * 1000 * sqrt(rank / max(rank))
  angle <- rank * golden_angle
  map_data$x_plot[idx] <- center_x + radius * cos(angle)
  map_data$y_plot[idx] <- center_y + radius * sin(angle)
}

q_columns <- colnames(q19)
selected_rows <- map_data$pop_gen %in% c("uniform", "admixed") &
  !is.na(map_data$pop_gen_population) & nzchar(map_data$pop_gen_population)
population_indices <- split(
  which(selected_rows), map_data$pop_gen_population[selected_rows]
)
population_rows <- lapply(names(population_indices), function(population) {
  idx <- population_indices[[population]]
  classes <- unique(map_data$pop_gen[idx])
  if (length(classes) != 1) stop("Population has mixed classes: ", population)
  lineages <- sort(unique(map_data$spID[idx]))
  values <- colMeans(map_data[idx, q_columns, drop = FALSE])
  data.frame(
    population = population,
    pop_gen = classes,
    species_lineage = paste(lineages, collapse = ";"),
    n_samples = length(idx),
    mean_latitude = mean(map_data$Latitude[idx]),
    mean_longitude = mean(map_data$Longitude[idx]),
    mean_total_minor_ancestry = mean(map_data$total_minor_ancestry[idx]),
    as.list(values), check.names = FALSE, stringsAsFactors = FALSE
  )
})
population_data <- do.call(rbind, population_rows)
population_data <- population_data[order(
  population_data$pop_gen, population_data$population
), , drop = FALSE]

component_label_file <- file.path(
  dirname(admixture_dir), "K19_component_labels.tsv"
)
if (!file.exists(component_label_file)) {
  stop("Missing K19 component label table: ", component_label_file)
}
component_key <- read.delim(
  component_label_file, check.names = FALSE, stringsAsFactors = FALSE
)
needed_component <- c("component", "label", "color")
if (length(setdiff(needed_component, names(component_key)))) {
  stop("K19_component_labels.tsv requires component, label, and color")
}
if (!setequal(component_key$component, q_columns)) {
  stop("K19_component_labels.tsv must contain Q01--Q19 exactly once")
}
component_key <- component_key[match(q_columns, component_key$component), ]
component_table <- data.frame(
  displayed_component = q_columns,
  original_K19_component = raw_component,
  label = component_key$label,
  color = component_key$color,
  stringsAsFactors = FALSE
)

write.table(
  map_data, file.path(output_dir, "aligned_K19_sample_map_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  population_data, file.path(output_dir, "aligned_K19_population_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  component_table, file.path(output_dir, "aligned_K19_components.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Samples exported: ", nrow(map_data))
message("Population summaries exported: ", nrow(population_data))
message("K19 map tables: ", output_dir)

