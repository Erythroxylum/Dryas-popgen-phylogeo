#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(paste(
    "Usage: 14_plot_phylogeography.R SAMPLE_MAP.tsv POPULATION_K19.tsv",
    "LGM_BOUNDARY.shp DIVERSITY_DIR REFUGIA_METRICS.tsv REGIONS.tsv OUTPUT.pdf"
  ))
}

required <- c(
  "ggplot2", "sf", "scatterpie", "rnaturalearthdata",
  "ggrepel", "ggnewscale", "patchwork"
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
map_min_latitude <- as.numeric(Sys.getenv("MAP_MIN_LATITUDE", "25"))
map_width <- as.numeric(Sys.getenv("MAP_WIDTH", "12"))
map_height <- as.numeric(Sys.getenv("MAP_HEIGHT", "10"))
figure_width <- as.numeric(Sys.getenv("FIGURE_WIDTH", "20"))
figure_height <- as.numeric(Sys.getenv("FIGURE_HEIGHT", "12.5"))
base_size <- as.numeric(Sys.getenv("FIGURE_BASE_SIZE", "9"))
map_land_color <- Sys.getenv("MAP_LAND_COLOR", "#F2F2F2")
map_sea_color <- Sys.getenv("MAP_SEA_COLOR", "#FFFFFF")
map_border_color <- Sys.getenv("MAP_BORDER_COLOR", "#777777")

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

samples <- read_tsv(sample_file)
population_q <- read_tsv(population_q_file)
refugia <- read_tsv(refugia_file)
regions <- read_tsv(regions_file)
component_label_file <- file.path(dirname(regions_file), "K19_component_labels.tsv")
component_key <- read_tsv(component_label_file)
pi_point <- read_tsv(file.path(diversity_dir, "filtered_sequence_population_pi.tsv"))
pi_ci <- read_tsv(file.path(diversity_dir, "filtered_sequence_pi_bootstrap_ci.tsv"))
snp_point <- read_tsv(file.path(diversity_dir, "population_snp_diversity.tsv"))

q_columns <- grep("^Q[0-9]{2}$", names(samples), value = TRUE)
if (length(q_columns) != 19) stop("Expected Q01--Q19 in sample map data")
needed_sample <- c(
  "sampleID", "pop_gen", "pop_gen_population", "x_plot", "y_plot", q_columns
)
if (length(setdiff(needed_sample, names(samples)))) {
  stop("Sample map table lacks required plotting columns")
}

needed_snp <- c("population", "pop_gen", "private_allelic_richness")
if (length(setdiff(needed_snp, names(snp_point)))) {
  stop("population_snp_diversity.tsv lacks private-allele columns")
}

# Classify diversity signatures using the across-population medians. These are
# descriptive quadrants rather than demographic-model assignments. Explicit
# cutoffs can be supplied for sensitivity analyses through environment values.
signature_data <- merge(
  pi_point[c("population", "pop_gen", "pi")],
  snp_point[c("population", "private_allelic_richness")],
  by = "population", all = FALSE, sort = FALSE
)
if (!nrow(signature_data) ||
    any(!is.finite(signature_data$pi)) ||
    any(!is.finite(signature_data$private_allelic_richness))) {
  stop("Non-finite pi or private allelic richness in signature data")
}

environment_cutoff <- function(name, values) {
  raw <- trimws(Sys.getenv(name, ""))
  if (!nzchar(raw)) return(median(values, na.rm = TRUE))
  value <- suppressWarnings(as.numeric(raw))
  if (!is.finite(value)) stop(name, " must be numeric")
  value
}

is_drummondii <- grepl("drummond", signature_data$population, ignore.case = TRUE)
classification_reference <- signature_data[!is_drummondii, , drop = FALSE]
if (!nrow(classification_reference)) stop("No non-drummondii populations remain")

pi_cutoff <- environment_cutoff(
  "SIGNATURE_PI_CUTOFF", classification_reference$pi
)
private_cutoff <- environment_cutoff(
  "SIGNATURE_PRIVATE_CUTOFF",
  classification_reference$private_allelic_richness
)

signature_data$signature <- ifelse(
  signature_data$pi < pi_cutoff &
    signature_data$private_allelic_richness >= private_cutoff,
  "isolated_relict",
  ifelse(
    signature_data$pi >= pi_cutoff &
      signature_data$private_allelic_richness >= private_cutoff,
    "persistent_refugium",
    ifelse(
      signature_data$pi < pi_cutoff &
        signature_data$private_allelic_richness < private_cutoff,
      "founder_population", "secondary_contact"
    )
  )
)
signature_data$signature[is_drummondii] <- "not_applicable"

signature_levels <- c(
  "isolated_relict", "persistent_refugium",
  "founder_population", "secondary_contact"
)
signature_labels <- c(
  isolated_relict = "Isolated relict",
  persistent_refugium = "Persistent refugium",
  founder_population = "Founder population",
  secondary_contact = "Secondary contact"
)
signature_colors <- c(
  isolated_relict = "#CC79A7",
  persistent_refugium = "#009E73",
  founder_population = "#0072B2",
  secondary_contact = "#D55E00"
)
plot_signature_colors <- c(signature_colors, not_applicable = "#777777")
plot_signature_labels <- c(
  signature_labels, not_applicable = "Drummondii (excluded)"
)
outline_colors <- c(plot_signature_colors, not_selected = "#B3B3B3")
outline_labels <- c(
  plot_signature_labels,
  not_selected = "Outside focal population comparison"
)

signature_data$signature <- factor(
  signature_data$signature, levels = c(signature_levels, "not_applicable")
)
samples$signature <- as.character(signature_data$signature[
  match(samples$pop_gen_population, signature_data$population)
])
samples$signature[samples$pop_gen == "not_selected"] <- "not_selected"
missing_signatures <- samples$sampleID[
  samples$pop_gen != "not_selected" & is.na(samples$signature)
]
if (length(missing_signatures)) {
  stop(
    "Selected samples lack a diversity-signature assignment: ",
    paste(missing_signatures, collapse = ", ")
  )
}

signature_output <- file.path(
  dirname(sample_file), "population_refugial_signatures.tsv"
)
signature_threshold_output <- file.path(
  dirname(sample_file), "population_refugial_signature_thresholds.tsv"
)
write.table(
  signature_data, signature_output,
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  data.frame(
    metric = c("pi", "private_allelic_richness"),
    cutoff = c(pi_cutoff, private_cutoff),
    default_rule = "median_across_selected_non_drummondii_populations",
    stringsAsFactors = FALSE
  ),
  signature_threshold_output,
  sep = "\t", quote = FALSE, row.names = FALSE
)
region_colors <- c(
  refugial_region = "#1B7837", admixture_region = "#D95F02"
)
region_labels <- c(
  refugial_region = "Refugial region",
  admixture_region = "Admixture region"
)
needed_component <- c("component", "label", "color")
if (length(setdiff(needed_component, names(component_key)))) {
  stop("K19_component_labels.tsv requires component, label, and color")
}
if (!setequal(component_key$component, q_columns)) {
  stop("K19_component_labels.tsv must contain Q01--Q19 exactly once")
}
component_key <- component_key[match(q_columns, component_key$component), ]
component_colors <- setNames(component_key$color, component_key$component)
component_labels <- setNames(component_key$label, component_key$component)

# Match the mapmixture basemap: Natural Earth 1:50m country polygons with
# visible national boundaries. The plot remains in Arctic Polar Stereographic
# coordinates so the full circum-Arctic distribution is represented evenly.
land <- rnaturalearthdata::countries50[, "geometry"]
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
  unknown <- setdiff(unique(regions$region_type), c("refugium", "contact", "mixed"))
  if (length(unknown)) stop("Unknown region types: ", paste(unknown, collapse = ", "))
  regions$region_class <- ifelse(
    regions$region_type == "refugium", "refugial_region", "admixture_region"
  )
  region_points <- st_as_sf(
    regions, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE
  )
  region_points <- st_transform(region_points, 3995)
  region_buffers <- st_buffer(region_points, dist = regions$radius_km * 1000)
}

# A projected latitude circle gives a symmetric polar crop.
crop_points <- st_as_sf(
  data.frame(lon = seq(-180, 180, by = 2), lat = map_min_latitude),
  coords = c("lon", "lat"), crs = 4326
)
crop_xy <- st_coordinates(st_transform(crop_points, 3995))
map_limit <- max(sqrt(rowSums(crop_xy^2)))

map_plot <- ggplot() +
  geom_sf(data = land, fill = map_land_color, color = NA) +
  geom_sf(
    data = lgm, fill = "#CFE8F3", color = "#6BAED6",
    linewidth = 0.25, alpha = 0.30
  ) +
  # Redraw present-day national and coastal boundaries above the LGM layer.
  geom_sf(
    data = land, fill = NA, color = map_border_color, linewidth = 0.18
  )

if (nrow(regions)) {
  map_plot <- map_plot +
    geom_sf(
      data = region_buffers,
      aes(color = region_class), fill = NA, linewidth = 0.65
    ) +
    scale_color_manual(
      values = region_colors,
      breaks = names(region_colors),
      labels = region_labels[names(region_colors)],
      name = "Geographic regions",
      guide = guide_legend(ncol = 1, order = 1)
    ) +
    ggnewscale::new_scale_color()
}

# Invisible points create an outline-color legend; the pie borders themselves
# are drawn in separate fixed-color layers below.
map_plot <- map_plot + geom_point(
  data = samples,
  aes(x = x_plot, y = y_plot, color = signature),
  alpha = 0, size = 0.01
)

for (signature_class in names(outline_colors)) {
  z <- samples[samples$signature == signature_class, , drop = FALSE]
  if (!nrow(z)) next
  map_plot <- map_plot + scatterpie::geom_scatterpie(
    data = z,
    aes(x = x_plot, y = y_plot, r = pie_radius_km * 1000),
    cols = q_columns,
    color = outline_colors[signature_class], linewidth = 0.38,
    alpha = 0.98
  )
}

map_plot <- map_plot +
  scale_fill_manual(
    values = component_colors,
    breaks = q_columns,
    labels = component_labels[q_columns],
    name = "K19 ancestry component"
  ) +
  scale_color_manual(
    values = outline_colors,
    breaks = signature_levels,
    labels = signature_labels[signature_levels],
    name = "Diversity signature"
  ) +
  guides(
    fill = guide_legend(ncol = 1, order = 3, override.aes = list(alpha = 1)),
    color = guide_legend(ncol = 1, order = 2,
                         override.aes = list(alpha = 1, size = 2.8))
  ) +
  coord_sf(
    crs = st_crs(3995),
    xlim = c(-map_limit, map_limit), ylim = c(-map_limit, map_limit),
    expand = FALSE
  ) +
  labs(
    title = "Geographic structure and secondary contact in Dryas",
    subtitle = "Sample-level K19 ancestry; pale blue indicates the supplied LGM boundary layer",
    x = "Longitude", y = "Latitude",
    caption = "Gray pie outlines indicate samples outside the focal population comparison."
  ) +
  theme_minimal(base_size = base_size) +
  theme(
    panel.background = element_rect(fill = map_sea_color, color = NA),
    panel.grid.major = element_line(color = "#D2D2D2", linewidth = 0.22),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.32),
    plot.background = element_rect(fill = map_sea_color, color = NA),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "left",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.text = element_text(size = base_size * 0.72),
    legend.title = element_text(size = base_size * 0.80, face = "bold"),
    legend.key.height = grid::unit(2.8, "mm"),
    legend.key.width = grid::unit(4.2, "mm"),
    plot.margin = margin(4, 1, 4, 4)
  )

# Population pi with its locus-bootstrap intervals.
pi_ci_small <- pi_ci[c("population", "lower_95", "upper_95")]
names(pi_ci_small)[-1] <- c("pi_lower_95", "pi_upper_95")
pi_data <- merge(
  pi_point[c("population", "pop_gen", "pi")], pi_ci_small,
  by = "population", all.x = TRUE, sort = FALSE
)
pi_data$signature <- signature_data$signature[
  match(pi_data$population, signature_data$population)
]
pi_data$pop_gen <- factor(pi_data$pop_gen, levels = c("uniform", "admixed"))

pi_class_plot <- ggplot(pi_data, aes(x = pop_gen, y = pi)) +
  geom_boxplot(width = 0.48, outlier.shape = NA, color = "#555555", fill = "white") +
  geom_jitter(
    aes(color = signature), width = 0.12, height = 0,
    size = 1.8, alpha = 0.82
  ) +
  scale_color_manual(values = plot_signature_colors, guide = "none") +
  scale_x_discrete(labels = c(uniform = "Uniform", admixed = "Admixed")) +
  labs(
    title = "Population nucleotide diversity",
    x = "Population type", y = expression(pi)
  ) +
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
private_data$signature <- signature_data$signature[
  match(private_data$population, signature_data$population)
]

private_plot <- ggplot(
  private_data,
  aes(
    x = within_lineage_private_allelic_richness, y = display_population,
    color = signature
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
  scale_color_manual(values = plot_signature_colors, guide = "none") +
  labs(
    title = "Locally unique variation",
    subtitle = "Uniform populations compared within lineage",
    x = "Within-lineage private allelic richness", y = "Population"
  ) +
  theme_classic(base_size = base_size) +
  theme(axis.text.y = element_text(size = base_size * 0.66))

# Link population diversity directly to non-dominant K19 ancestry. This is a
# descriptive ADMIXTURE-derived index, not an independent test of gene flow.
relationship <- merge(
  pi_data, population_q[c("population", "mean_total_minor_ancestry")],
  by = "population", all = FALSE, sort = FALSE
)
admixture_geography <- function(x) {
  x <- sub("^x_?", "", x)
  x <- sub("^[^_]+_", "", x)
  gsub("_", " ", x)
}
relationship$label <- ifelse(
  relationship$pop_gen == "admixed",
  admixture_geography(relationship$population), ""
)

relationship_plot <- ggplot(
  relationship,
  aes(x = mean_total_minor_ancestry, y = pi, color = signature)
) +
  geom_point(size = 2, alpha = 0.82) +
  geom_text_repel(
    aes(label = label), size = 2.15, seed = 20261014,
    min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE,
    box.padding = 0.28, point.padding = 0.12, force = 1.4
  ) +
  scale_color_manual(
    values = plot_signature_colors, guide = "none"
  ) +
  labs(
    title = "Diversity and mixed ancestry",
    x = "Mean non-dominant K19 ancestry", y = expression(pi)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.12))) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = base_size) +
  theme(
    legend.position = "none",
    plot.margin = margin(4, 20, 4, 1)
  )

# Standalone four-quadrant summary. Drummondii is excluded because it is not
# equivalent to the ingroup populations in the refugial comparison.
signature_plot_data <- signature_data[!is_drummondii, , drop = FALSE]
percentile_rank <- function(x) {
  if (length(x) <= 1L || length(unique(x)) <= 1L) return(rep(50, length(x)))
  100 * (rank(x, ties.method = "average") - 1) / (length(x) - 1)
}
signature_plot_data$pi_percentile <- percentile_rank(signature_plot_data$pi)
signature_plot_data$private_percentile <- percentile_rank(
  signature_plot_data$private_allelic_richness
)
quadrant_text <- data.frame(
  signature = factor(signature_levels, levels = signature_levels),
  x = c(3, 97, 3, 97),
  y = c(97, 97, 3, 3),
  hjust = c(0, 1, 0, 1),
  vjust = c(1, 1, 0, 0),
  label = c(
    "Isolated relict:\nlow pi, high private alleles",
    "Persistent refugium:\nhigh pi, high private alleles",
    "Founder population:\nlow pi, low private alleles",
    "Secondary contact:\nhigh pi, low private alleles"
  )
)

signature_plot <- ggplot(
  signature_plot_data,
  aes(x = pi_percentile, y = private_percentile)
) +
  geom_vline(xintercept = 50, color = "#4D4D4D", linewidth = 0.45) +
  geom_hline(yintercept = 50, color = "#4D4D4D", linewidth = 0.45) +
  geom_point(
    aes(color = signature), shape = 21, fill = "white",
    size = 2.6, stroke = 0.95, alpha = 1
  ) +
  geom_text(
    data = quadrant_text,
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    color = "black",
    fontface = "bold", size = 2.6, lineheight = 0.95,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_color_manual(
    values = signature_colors,
    breaks = signature_levels,
    labels = signature_labels[signature_levels],
    name = "Diversity signature"
  ) +
  labs(
    title = "Population diversity signatures",
    subtitle = "Rank-standardized diversity space; drummondii excluded",
    x = expression("Nucleotide diversity (" * pi * ") percentile"),
    y = "Private allelic richness percentile",
    caption = paste0(
      "Descriptive hypotheses; raw-value median cutoffs exclude drummondii. ",
      "Private alleles are rarefied."
    )
  ) +
  scale_x_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 25),
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 25),
    expand = expansion(mult = 0.02)
  ) +
  theme_classic(base_size = base_size + 1) +
  guides(
    color = guide_legend(
      ncol = 2,
      override.aes = list(shape = 21, fill = "white", size = 2.6)
    )
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 10.5),
    plot.subtitle = element_text(size = 7.5),
    plot.caption = element_text(hjust = 0, size = 5.8),
    axis.title = element_text(size = 8.2),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 7.2),
    legend.text = element_text(size = 6.8),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  )

signature_pdf <- sub(
  "\\.pdf$", "_refugial_signatures.pdf", output_file, ignore.case = TRUE
)
ggsave(signature_pdf, signature_plot, width = 4, height = 5, units = "in")

map_only <- sub("\\.pdf$", "_map.pdf", output_file, ignore.case = TRUE)
ggsave(map_only, map_plot, width = map_width, height = map_height, units = "in")

pi_class_plot <- pi_class_plot + theme(plot.margin = margin(4, 8, 4, 1))
private_plot <- private_plot + theme(plot.margin = margin(4, 8, 4, 1))
side_panels <- pi_class_plot / private_plot / relationship_plot +
  plot_layout(heights = c(0.75, 1.25, 1.45))
composite <- map_plot | side_panels
composite <- composite + plot_layout(widths = c(2.0, 1.08)) +
  plot_annotation(tag_levels = "A")
ggsave(
  output_file, composite,
  width = figure_width, height = figure_height, units = "in"
)

message("Map-only figure: ", map_only)
message("Composite figure: ", output_file)
message("Refugial-signature figure: ", signature_pdf)
message("Refugial-signature assignments: ", signature_output)
message("Refugial-signature thresholds: ", signature_threshold_output)

