#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(paste(
    "Usage: 14_plot_phylogeography.R SAMPLE_MAP.tsv POPULATION_K19.tsv",
    "LGM_BOUNDARY.shp DIVERSITY_DIR REFUGIA_METRICS.tsv REGIONS.tsv OUTPUT.pdf"
  ))
}

required <- c(
  "ggplot2", "sf", "scatterpie", "rnaturalearth", "rnaturalearthdata",
  "ggrepel", "patchwork", "viridisLite"
)
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(scatterpie))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(patchwork))

sample_file <- args[1]
population_q_file <- args[2]
lgm_file <- args[3]
diversity_dir <- args[4]
refugia_file <- args[5]
regions_file <- args[6]
output_file <- args[7]
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

pie_radius_km <- as.numeric(Sys.getenv("MAP_PIE_RADIUS_KM", "60"))
map_min_latitude <- as.numeric(Sys.getenv("MAP_MIN_LATITUDE", "35"))
map_width <- as.numeric(Sys.getenv("MAP_WIDTH", "10"))
map_height <- as.numeric(Sys.getenv("MAP_HEIGHT", "10"))
figure_width <- as.numeric(Sys.getenv("FIGURE_WIDTH", "17"))
figure_height <- as.numeric(Sys.getenv("FIGURE_HEIGHT", "11"))
base_size <- as.numeric(Sys.getenv("FIGURE_BASE_SIZE", "9"))

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

samples <- read_tsv(sample_file)
population_q <- read_tsv(population_q_file)
refugia <- read_tsv(refugia_file)
regions <- read_tsv(regions_file)
pi_point <- read_tsv(file.path(diversity_dir, "filtered_sequence_population_pi.tsv"))
pi_ci <- read_tsv(file.path(diversity_dir, "filtered_sequence_pi_bootstrap_ci.tsv"))

q_columns <- grep("^Q[0-9]{2}$", names(samples), value = TRUE)
if (length(q_columns) != 19) stop("Expected Q01--Q19 in sample map data")
needed_sample <- c("sampleID", "pop_gen", "x_plot", "y_plot", q_columns)
if (length(setdiff(needed_sample, names(samples)))) {
  stop("Sample map table lacks required plotting columns")
}

class_colors <- c(
  uniform = "#2166AC", admixed = "#B2182B", not_selected = "#666666"
)
region_colors <- c(
  refugium = "#1B7837", contact = "#D95F02", mixed = "#7B3294"
)
component_colors <- setNames(viridisLite::turbo(22)[seq_len(19)], q_columns)

land <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")
land <- st_transform(st_make_valid(land), 3995)
lgm <- st_read(lgm_file, quiet = TRUE)
if (is.na(st_crs(lgm))) stop("LGM boundary has no coordinate reference system")
lgm <- st_transform(st_make_valid(lgm), 3995)

show_region <- tolower(trimws(as.character(regions$show))) %in% c("true", "t", "1", "yes")
regions <- regions[show_region, , drop = FALSE]
if (nrow(regions)) {
  required_regions <- c(
    "label", "region_type", "longitude", "latitude", "radius_km"
  )
  absent <- setdiff(required_regions, names(regions))
  if (length(absent)) stop("Region columns missing: ", paste(absent, collapse = ", "))
  unknown <- setdiff(unique(regions$region_type), names(region_colors))
  if (length(unknown)) stop("Unknown region types: ", paste(unknown, collapse = ", "))
  region_points <- st_as_sf(
    regions, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE
  )
  region_points <- st_transform(region_points, 3995)
  region_buffers <- st_buffer(region_points, dist = regions$radius_km * 1000)
  region_xy <- st_coordinates(region_points)
  region_labels <- data.frame(
    label = regions$label, region_type = regions$region_type,
    x = region_xy[, 1], y = region_xy[, 2], stringsAsFactors = FALSE
  )
}

# A projected latitude circle gives a symmetric polar crop.
crop_points <- st_as_sf(
  data.frame(lon = seq(-180, 180, by = 2), lat = map_min_latitude),
  coords = c("lon", "lat"), crs = 4326
)
crop_xy <- st_coordinates(st_transform(crop_points, 3995))
map_limit <- max(sqrt(rowSums(crop_xy^2)))

map_plot <- ggplot() +
  geom_sf(data = land, fill = "#D9D9D9", color = "#B8B8B8", linewidth = 0.12) +
  geom_sf(
    data = lgm, fill = "#CFE8F3", color = "#6BAED6",
    linewidth = 0.25, alpha = 0.30
  )

if (nrow(regions)) {
  for (region_type in names(region_colors)) {
    keep <- regions$region_type == region_type
    if (!any(keep)) next
    map_plot <- map_plot + geom_sf(
      data = region_buffers[keep, ],
      fill = grDevices::adjustcolor(region_colors[region_type], alpha.f = 0.07),
      color = region_colors[region_type], linewidth = 0.45,
      linetype = if (region_type == "contact") "dashed" else "solid"
    )
    map_plot <- map_plot + geom_text_repel(
      data = region_labels[region_labels$region_type == region_type, ],
      aes(x = x, y = y, label = label),
      color = region_colors[region_type],
      size = 2.6, seed = 20261014, min.segment.length = 0,
      box.padding = 0.25, point.padding = 0.1, max.overlaps = Inf,
      show.legend = FALSE
    )
  }
}

# Invisible points create an outline-color legend; the pie borders themselves
# are drawn in separate fixed-color layers below.
map_plot <- map_plot + geom_point(
  data = samples,
  aes(x = x_plot, y = y_plot, color = pop_gen),
  alpha = 0, size = 0.01
)

for (population_class in names(class_colors)) {
  z <- samples[samples$pop_gen == population_class, , drop = FALSE]
  if (!nrow(z)) next
  map_plot <- map_plot + scatterpie::geom_scatterpie(
    data = z,
    aes(x = x_plot, y = y_plot, r = pie_radius_km * 1000),
    cols = q_columns,
    color = class_colors[population_class], linewidth = 0.28,
    alpha = 0.98
  )
}

map_plot <- map_plot +
  scale_fill_manual(
    values = component_colors,
    breaks = q_columns,
    labels = paste("Component", seq_along(q_columns)),
    name = "K19 ancestry component"
  ) +
  scale_color_manual(
    values = class_colors,
    breaks = names(class_colors),
    labels = c("Uniform population", "Admixed population", "Not selected"),
    name = "Population assignment"
  ) +
  guides(
    fill = guide_legend(ncol = 2, order = 1, override.aes = list(alpha = 1)),
    color = guide_legend(order = 2, override.aes = list(alpha = 1, size = 2.8))
  ) +
  coord_sf(
    crs = st_crs(3995),
    xlim = c(-map_limit, map_limit), ylim = c(-map_limit, map_limit),
    expand = FALSE
  ) +
  labs(
    title = "Geographic structure and secondary contact in Dryas",
    subtitle = "Sample-level K19 ancestry; pale blue indicates the supplied LGM boundary layer",
    x = NULL, y = NULL
  ) +
  theme_void(base_size = base_size) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.height = grid::unit(3.2, "mm"),
    legend.key.width = grid::unit(4.2, "mm")
  )

# Population pi with its locus-bootstrap intervals.
pi_ci_small <- pi_ci[c("population", "lower_95", "upper_95")]
names(pi_ci_small)[-1] <- c("pi_lower_95", "pi_upper_95")
pi_data <- merge(
  pi_point[c("population", "pop_gen", "pi")], pi_ci_small,
  by = "population", all.x = TRUE, sort = FALSE
)
pi_data$pop_gen <- factor(pi_data$pop_gen, levels = c("uniform", "admixed"))

pi_class_plot <- ggplot(pi_data, aes(x = pop_gen, y = pi, color = pop_gen)) +
  geom_boxplot(width = 0.48, outlier.shape = NA, color = "#555555", fill = "white") +
  geom_jitter(width = 0.12, height = 0, size = 1.7, alpha = 0.78) +
  scale_color_manual(values = class_colors[c("uniform", "admixed")], guide = "none") +
  scale_x_discrete(labels = c(uniform = "Uniform", admixed = "Admixed")) +
  labs(title = "Population nucleotide diversity", x = NULL, y = expression(pi)) +
  theme_classic(base_size = base_size)

# Within-lineage private allelic richness is defined only for lineages with at
# least two uniform reference populations.
private_data <- refugia[
  refugia$status == "estimated" &
    is.finite(refugia$within_lineage_private_allelic_richness),,
  drop = FALSE
]
private_data$display_population <- gsub("_", " ", private_data$population)
private_data$display_population <- reorder(
  private_data$display_population,
  private_data$within_lineage_private_allelic_richness
)
lineages <- sort(unique(private_data$species_lineage))
lineage_colors <- setNames(hcl.colors(length(lineages), "Dark 3"), lineages)

private_plot <- ggplot(
  private_data,
  aes(
    x = within_lineage_private_allelic_richness,
    y = display_population, color = species_lineage
  )
) +
  geom_segment(
    aes(
      x = private_lower_95, xend = private_upper_95,
      y = display_population, yend = display_population
    ),
    linewidth = 0.55
  ) +
  geom_point(size = 1.8) +
  scale_color_manual(values = lineage_colors, guide = "none") +
  labs(
    title = "Locally unique variation",
    subtitle = "Uniform populations compared within lineage",
    x = "Within-lineage private allelic richness", y = NULL
  ) +
  theme_classic(base_size = base_size) +
  theme(axis.text.y = element_text(size = base_size * 0.66))

# Link population diversity directly to non-dominant K19 ancestry. This is a
# descriptive ADMIXTURE-derived index, not an independent test of gene flow.
relationship <- merge(
  pi_data, population_q[c("population", "mean_total_minor_ancestry")],
  by = "population", all = FALSE, sort = FALSE
)
label_order <- order(
  rank(-relationship$pi) + rank(-relationship$mean_total_minor_ancestry)
)
relationship$label <- ""
relationship$label[label_order[seq_len(min(8, nrow(relationship)))]] <-
  gsub("_", " ", relationship$population[label_order[seq_len(min(8, nrow(relationship)))]] )

relationship_plot <- ggplot(
  relationship,
  aes(x = mean_total_minor_ancestry, y = pi, color = pop_gen)
) +
  geom_point(size = 2, alpha = 0.82) +
  geom_text_repel(
    aes(label = label), size = 2.4, seed = 20261014,
    min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_color_manual(
    values = class_colors[c("uniform", "admixed")],
    labels = c(uniform = "Uniform", admixed = "Admixed"),
    name = NULL
  ) +
  labs(
    title = "Diversity and mixed ancestry",
    x = "Mean non-dominant K19 ancestry", y = expression(pi)
  ) +
  theme_classic(base_size = base_size) +
  theme(legend.position = "top")

map_only <- sub("\\.pdf$", "_map.pdf", output_file, ignore.case = TRUE)
ggsave(map_only, map_plot, width = map_width, height = map_height, units = "in")

side_panels <- pi_class_plot / private_plot / relationship_plot +
  plot_layout(heights = c(0.78, 1.45, 1.0))
composite <- map_plot | side_panels
composite <- composite + plot_layout(widths = c(2.25, 1.0)) +
  plot_annotation(tag_levels = "A")
ggsave(
  output_file, composite,
  width = figure_width, height = figure_height, units = "in"
)

message("Map-only figure: ", map_only)
message("Composite figure: ", output_file)

