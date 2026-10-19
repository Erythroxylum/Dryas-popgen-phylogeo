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
  "ggrepel", "ggnewscale", "patchwork", "scales"
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
unassigned_pie_scale <- as.numeric(Sys.getenv(
  "MAP_UNASSIGNED_PIE_SCALE", "0.70"
))
unassigned_pie_alpha <- as.numeric(Sys.getenv(
  "MAP_UNASSIGNED_PIE_ALPHA", "0.35"
))
population_halo_min_radius_km <- as.numeric(Sys.getenv(
  "MAP_POPULATION_HALO_MIN_RADIUS_KM", "60"
))
population_halo_padding_km <- as.numeric(Sys.getenv(
  "MAP_POPULATION_HALO_PADDING_KM", "20"
))
population_halo_linewidth <- as.numeric(Sys.getenv(
  "MAP_POPULATION_HALO_LINEWIDTH", "0.75"
))
population_halo_alpha <- as.numeric(Sys.getenv(
  "MAP_POPULATION_HALO_ALPHA", "0.90"
))
show_hypothesis_regions <- tolower(trimws(Sys.getenv(
  "MAP_SHOW_HYPOTHESIS_REGIONS", "FALSE"
))) %in% c("true", "t", "1", "yes")
merge_same_signature_halos <- tolower(trimws(Sys.getenv(
  "MAP_MERGE_SAME_SIGNATURE_HALOS", "TRUE"
))) %in% c("true", "t", "1", "yes")
show_leader_lines <- tolower(trimws(Sys.getenv(
  "MAP_SHOW_LEADER_LINES", "TRUE"
))) %in% c("true", "t", "1", "yes")
leader_line_color <- Sys.getenv("MAP_LEADER_LINE_COLOR", "#595959")
leader_line_width <- as.numeric(Sys.getenv("MAP_LEADER_LINE_WIDTH", "0.24"))
leader_line_alpha <- as.numeric(Sys.getenv("MAP_LEADER_LINE_ALPHA", "0.70"))
coordinate_point_size <- as.numeric(Sys.getenv(
  "MAP_COORDINATE_POINT_SIZE", "0.45"
))
beringia_pie_radius_km <- as.numeric(Sys.getenv(
  "BERINGIA_PIE_RADIUS_KM", as.character(pie_radius_km)
))
map_min_latitude <- as.numeric(Sys.getenv("MAP_MIN_LATITUDE", "25"))
map_center_lon <- as.numeric(Sys.getenv("MAP_CENTER_LON", "0"))
map_center_lat <- as.numeric(Sys.getenv("MAP_CENTER_LAT", "90"))
map_half_width_x_raw <- trimws(Sys.getenv("MAP_HALF_WIDTH_X_KM", ""))
map_half_height_raw <- trimws(Sys.getenv("MAP_HALF_HEIGHT_KM", ""))
map_width <- as.numeric(Sys.getenv("MAP_WIDTH", "12"))
map_height <- as.numeric(Sys.getenv("MAP_HEIGHT", "10"))
figure_width <- as.numeric(Sys.getenv("FIGURE_WIDTH", "20"))
figure_height <- as.numeric(Sys.getenv("FIGURE_HEIGHT", "12.5"))
base_size <- as.numeric(Sys.getenv("FIGURE_BASE_SIZE", "9"))
axis_title_size <- as.numeric(Sys.getenv(
  "FIGURE_AXIS_TITLE_SIZE", as.character(base_size * 1.05)
))
axis_text_size <- as.numeric(Sys.getenv(
  "FIGURE_AXIS_TEXT_SIZE", as.character(base_size * 0.90)
))
signature_width <- as.numeric(Sys.getenv("SIGNATURE_WIDTH", "4"))
signature_height <- as.numeric(Sys.getenv("SIGNATURE_HEIGHT", "5"))
signature_x_expand_left <- as.numeric(Sys.getenv(
  "SIGNATURE_X_EXPAND_LEFT", "0.04"
))
signature_x_expand_right <- as.numeric(Sys.getenv(
  "SIGNATURE_X_EXPAND_RIGHT", "0.18"
))
signature_y_expand_bottom <- as.numeric(Sys.getenv(
  "SIGNATURE_Y_EXPAND_BOTTOM", "0.18"
))
signature_y_expand_top <- as.numeric(Sys.getenv(
  "SIGNATURE_Y_EXPAND_TOP", "0.04"
))
signature_text_x_inset <- as.numeric(Sys.getenv(
  "SIGNATURE_TEXT_X_INSET", "0.02"
))
signature_text_y_inset <- as.numeric(Sys.getenv(
  "SIGNATURE_TEXT_Y_INSET", "0.02"
))
signature_text_right_hjust <- as.numeric(Sys.getenv(
  "SIGNATURE_TEXT_RIGHT_HJUST", "0.80"
))
signature_si_width <- as.numeric(Sys.getenv("SIGNATURE_SI_WIDTH", "8"))
signature_si_height <- as.numeric(Sys.getenv("SIGNATURE_SI_HEIGHT", "7"))
signature_si_label_size <- as.numeric(Sys.getenv(
  "SIGNATURE_SI_LABEL_SIZE", "2.8"
))
heterozygosity_point_size <- as.numeric(Sys.getenv(
  "HETEROZYGOSITY_POINT_SIZE", "1.35"
))
map_land_color <- Sys.getenv("MAP_LAND_COLOR", "#F2F2F2")
map_sea_color <- Sys.getenv("MAP_SEA_COLOR", "#FFFFFF")
map_border_color <- Sys.getenv("MAP_BORDER_COLOR", "#777777")
lgm_fill_color <- Sys.getenv("LGM_FILL_COLOR", "#BDD7E7")
lgm_fill_alpha <- as.numeric(Sys.getenv("LGM_FILL_ALPHA", "0.70"))
lgm_min_latitude <- as.numeric(Sys.getenv("LGM_MIN_LATITUDE", "0"))
lgm_line_color <- Sys.getenv("LGM_LINE_COLOR", "#6BAED6")
lgm_line_width <- as.numeric(Sys.getenv("LGM_LINE_WIDTH", "0.35"))
lgm_line_alpha <- as.numeric(Sys.getenv("LGM_LINE_ALPHA", "0.65"))
beringia_center_lon <- as.numeric(Sys.getenv("BERINGIA_CENTER_LON", "180"))
beringia_center_lat <- as.numeric(Sys.getenv("BERINGIA_CENTER_LAT", "65"))
beringia_half_width_km <- as.numeric(Sys.getenv("BERINGIA_HALF_WIDTH_KM", "2250"))
beringia_half_width_x_km <- as.numeric(Sys.getenv(
  "BERINGIA_HALF_WIDTH_X_KM", as.character(beringia_half_width_km * 1.40)
))
beringia_half_height_km <- as.numeric(Sys.getenv(
  "BERINGIA_HALF_HEIGHT_KM", as.character(beringia_half_width_km * 0.80)
))
beringia_jitter_radius_km <- as.numeric(Sys.getenv("BERINGIA_JITTER_RADIUS_KM", "45"))
beringia_map_width <- as.numeric(Sys.getenv("BERINGIA_MAP_WIDTH", "12"))
beringia_map_height <- as.numeric(Sys.getenv("BERINGIA_MAP_HEIGHT", "9"))
beringia_legend_ncol <- as.integer(Sys.getenv("BERINGIA_LEGEND_NCOL", "4"))

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
individual_heterozygosity <- read_tsv(file.path(
  diversity_dir, "individual_snp_heterozygosity.tsv"
))

q_columns <- grep("^Q[0-9]{2}$", names(samples), value = TRUE)
if (length(q_columns) != 19) stop("Expected Q01--Q19 in sample map data")
needed_sample <- c(
  "sampleID", "geoID", "Latitude", "Longitude", "pop_gen",
  "pop_gen_population", "x_3995", "y_3995", "x_plot", "y_plot", q_columns
)
if (length(setdiff(needed_sample, names(samples)))) {
  stop("Sample map table lacks required plotting columns")
}
if (!is.finite(lgm_fill_alpha) || lgm_fill_alpha < 0 || lgm_fill_alpha > 1) {
  stop("LGM_FILL_ALPHA must be between zero and one")
}
if (!is.finite(unassigned_pie_scale) || unassigned_pie_scale <= 0 ||
    unassigned_pie_scale > 1) {
  stop("MAP_UNASSIGNED_PIE_SCALE must be greater than 0 and no greater than 1")
}
if (!is.finite(unassigned_pie_alpha) || unassigned_pie_alpha < 0 ||
    unassigned_pie_alpha > 1) {
  stop("MAP_UNASSIGNED_PIE_ALPHA must be between 0 and 1")
}
if (!is.finite(population_halo_min_radius_km) ||
    population_halo_min_radius_km <= 0) {
  stop("MAP_POPULATION_HALO_MIN_RADIUS_KM must be positive")
}
if (!is.finite(population_halo_padding_km) ||
    population_halo_padding_km < 0) {
  stop("MAP_POPULATION_HALO_PADDING_KM must be non-negative")
}
if (!is.finite(population_halo_linewidth) || population_halo_linewidth <= 0) {
  stop("MAP_POPULATION_HALO_LINEWIDTH must be positive")
}
if (!is.finite(population_halo_alpha) || population_halo_alpha < 0 ||
    population_halo_alpha > 1) {
  stop("MAP_POPULATION_HALO_ALPHA must be between 0 and 1")
}
if (!is.finite(leader_line_width) || leader_line_width <= 0) {
  stop("MAP_LEADER_LINE_WIDTH must be positive")
}
if (!is.finite(leader_line_alpha) || leader_line_alpha < 0 ||
    leader_line_alpha > 1) {
  stop("MAP_LEADER_LINE_ALPHA must be between 0 and 1")
}
if (!is.finite(coordinate_point_size) || coordinate_point_size < 0) {
  stop("MAP_COORDINATE_POINT_SIZE must be non-negative")
}
if (!is.finite(lgm_min_latitude) || lgm_min_latitude < -90 ||
    lgm_min_latitude > 90) {
  stop("LGM_MIN_LATITUDE must be between -90 and 90")
}
if (!is.finite(map_center_lon) || map_center_lon < -180 ||
    map_center_lon > 180) {
  stop("MAP_CENTER_LON must be between -180 and 180")
}
if (!is.finite(map_center_lat) || map_center_lat < -90 ||
    map_center_lat > 90) {
  stop("MAP_CENTER_LAT must be between -90 and 90")
}
if (!is.finite(lgm_line_width) || lgm_line_width <= 0) {
  stop("LGM_LINE_WIDTH must be a positive number")
}
if (!is.finite(lgm_line_alpha) || lgm_line_alpha < 0 || lgm_line_alpha > 1) {
  stop("LGM_LINE_ALPHA must be between zero and one")
}
signature_controls <- c(
  SIGNATURE_WIDTH = signature_width,
  SIGNATURE_HEIGHT = signature_height,
  SIGNATURE_X_EXPAND_LEFT = signature_x_expand_left,
  SIGNATURE_X_EXPAND_RIGHT = signature_x_expand_right,
  SIGNATURE_Y_EXPAND_BOTTOM = signature_y_expand_bottom,
  SIGNATURE_Y_EXPAND_TOP = signature_y_expand_top,
  SIGNATURE_TEXT_X_INSET = signature_text_x_inset,
  SIGNATURE_TEXT_Y_INSET = signature_text_y_inset,
  SIGNATURE_TEXT_RIGHT_HJUST = signature_text_right_hjust
)
if (any(!is.finite(signature_controls)) ||
    any(signature_controls[c("SIGNATURE_WIDTH", "SIGNATURE_HEIGHT")] <= 0) ||
    any(signature_controls[grep("EXPAND|INSET", names(signature_controls))] < 0) ||
    any(signature_controls[grep("TEXT_.*INSET", names(signature_controls))] >= 0.5)) {
  stop(paste(
    "Signature dimensions must be positive; axis expansions and text insets",
    "must be non-negative, and text insets must be below 0.5"
  ))
}
if (signature_text_right_hjust < 0 || signature_text_right_hjust > 1) {
  stop("SIGNATURE_TEXT_RIGHT_HJUST must be between 0 and 1")
}
if (!is.finite(signature_si_width) || signature_si_width <= 0 ||
    !is.finite(signature_si_height) || signature_si_height <= 0 ||
    !is.finite(signature_si_label_size) || signature_si_label_size <= 0) {
  stop("SIGNATURE_SI_WIDTH, SIGNATURE_SI_HEIGHT, and SIGNATURE_SI_LABEL_SIZE must be positive")
}
if (!is.finite(heterozygosity_point_size) || heterozygosity_point_size <= 0) {
  stop("HETEROZYGOSITY_POINT_SIZE must be positive")
}
if (!is.finite(axis_title_size) || axis_title_size <= 0 ||
    !is.finite(axis_text_size) || axis_text_size <= 0) {
  stop("FIGURE_AXIS_TITLE_SIZE and FIGURE_AXIS_TEXT_SIZE must be positive")
}

needed_individual <- c("sampleID", "snp_heterozygosity", "call_rate")
if (length(setdiff(needed_individual, names(individual_heterozygosity)))) {
  stop("individual_snp_heterozygosity.tsv lacks required columns")
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
  signature_data$private_allelic_richness >= private_cutoff,
  "candidate_refugium",
  ifelse(
    signature_data$pi < pi_cutoff,
    "founder_population", "secondary_contact"
  )
)
signature_data$signature[is_drummondii] <- "not_applicable"

signature_levels <- c(
  "candidate_refugium", "founder_population", "secondary_contact"
)
signature_labels <- c(
  candidate_refugium = "Candidate refugium",
  founder_population = "Founder population",
  secondary_contact = "Secondary contact"
)
signature_colors <- c(
  candidate_refugium = "#009E73",
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
lgm_wgs84 <- st_transform(lgm, 4326)
lgm_ymax <- vapply(
  st_geometry(lgm_wgs84),
  function(g) unname(st_bbox(g)["ymax"]),
  numeric(1)
)
keep_lgm <- is.finite(lgm_ymax) & lgm_ymax >= lgm_min_latitude
lgm_wgs84 <- lgm_wgs84[keep_lgm, , drop = FALSE]
if (!nrow(lgm_wgs84)) {
  stop("No LGM features remain north of LGM_MIN_LATITUDE")
}
lgm_wgs84 <- st_make_valid(lgm_wgs84)
lgm_geometry_types <- sort(unique(as.character(st_geometry_type(lgm_wgs84))))
lgm_has_polygon <- any(
  lgm_geometry_types %in% c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION")
)
lgm <- st_transform(lgm_wgs84, 3995)
message(
  "Northern LGM features retained: ", nrow(lgm_wgs84),
  "; geometry types: ",
  paste(lgm_geometry_types, collapse = ", ")
)
if (!lgm_has_polygon) {
  message(
    "LGM input is line geometry: drawing a blue boundary layer; ",
    "LGM_FILL_COLOR and LGM_FILL_ALPHA do not apply."
  )
}

show_region <- show_hypothesis_regions &
  tolower(trimws(as.character(regions$show))) %in% c("true", "t", "1", "yes")
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

build_population_halos <- function(data, x_column, y_column, crs, pie_radius_km) {
  population_name <- trimws(as.character(data$pop_gen_population))
  keep <- data$pop_gen != "not_selected" &
    !is.na(population_name) & nzchar(population_name) &
    is.finite(data[[x_column]]) & is.finite(data[[y_column]])
  index_groups <- split(which(keep), population_name[keep])
  if (!length(index_groups)) {
    return(st_sf(
      population = character(), signature = character(),
      radius_km = numeric(), geometry = st_sfc(crs = crs)
    ))
  }
  halo_list <- lapply(names(index_groups), function(population) {
    idx <- index_groups[[population]]
    xy <- cbind(data[[x_column]][idx], data[[y_column]][idx])
    center <- colMeans(xy)
    distances <- sqrt(rowSums(sweep(xy, 2, center, "-")^2))
    radius_m <- max(
      population_halo_min_radius_km * 1000,
      max(distances) +
        (pie_radius_km + population_halo_padding_km) * 1000
    )
    population_signature <- unique(na.omit(as.character(data$signature[idx])))
    if (length(population_signature) != 1) {
      stop("Population has inconsistent diversity signatures: ", population)
    }
    st_sf(
      population = population,
      signature = population_signature,
      radius_km = radius_m / 1000,
      geometry = st_buffer(
        st_sfc(st_point(center), crs = crs), dist = radius_m
      )
    )
  })
  do.call(rbind, halo_list)
}

dissolve_population_halos <- function(halos) {
  # Drummondii is retained as a sample pie but has no refugial-signature halo.
  halos <- halos[halos$signature %in% signature_levels, , drop = FALSE]
  if (!nrow(halos) || !merge_same_signature_halos) return(halos)
  dissolved <- lapply(signature_levels, function(signature_class) {
    z <- halos[halos$signature == signature_class, , drop = FALSE]
    if (!nrow(z)) return(NULL)
    st_sf(
      signature = signature_class,
      population_count = nrow(z),
      populations = paste(z$population, collapse = ";"),
      geometry = st_union(st_geometry(z))
    )
  })
  do.call(rbind, dissolved[!vapply(dissolved, is.null, logical(1))])
}

map_population_halos <- build_population_halos(
  samples, "x_plot", "y_plot", st_crs(3995), pie_radius_km
)
map_population_halos <- dissolve_population_halos(map_population_halos)

build_leader_segments <- function(data, origin_x, origin_y, display_x, display_y) {
  distance <- sqrt(
    (data[[display_x]] - data[[origin_x]])^2 +
      (data[[display_y]] - data[[origin_y]])^2
  )
  keep <- is.finite(distance) & distance > 1
  data.frame(
    x = data[[origin_x]][keep], y = data[[origin_y]][keep],
    xend = data[[display_x]][keep], yend = data[[display_y]][keep],
    stringsAsFactors = FALSE
  )
}

map_leaders <- build_leader_segments(
  samples, "x_3995", "y_3995", "x_plot", "y_plot"
)

# A projected latitude circle gives a symmetric polar crop.
crop_points <- st_as_sf(
  data.frame(lon = seq(-180, 180, by = 2), lat = map_min_latitude),
  coords = c("lon", "lat"), crs = 4326
)
crop_xy <- st_coordinates(st_transform(crop_points, 3995))
map_limit <- max(sqrt(rowSums(crop_xy^2)))

# Keep the Arctic Polar Stereographic projection but optionally shift and crop
# its plotting window. This changes the composition, not the map projection or
# the coordinates used for the ancestry pies.
map_center_point <- st_as_sf(
  data.frame(lon = map_center_lon, lat = map_center_lat),
  coords = c("lon", "lat"), crs = 4326
)
map_center_xy <- st_coordinates(st_transform(map_center_point, 3995))[1, ]
map_half_width_x <- if (nzchar(map_half_width_x_raw)) {
  as.numeric(map_half_width_x_raw) * 1000
} else {
  map_limit
}
map_half_height <- if (nzchar(map_half_height_raw)) {
  as.numeric(map_half_height_raw) * 1000
} else {
  map_limit
}
if (!is.finite(map_half_width_x) || map_half_width_x <= 0) {
  stop("MAP_HALF_WIDTH_X_KM must be a positive number")
}
if (!is.finite(map_half_height) || map_half_height <= 0) {
  stop("MAP_HALF_HEIGHT_KM must be a positive number")
}
map_xlim <- map_center_xy[1] + c(-map_half_width_x, map_half_width_x)
map_ylim <- map_center_xy[2] + c(-map_half_height, map_half_height)

map_plot <- ggplot() +
  geom_sf(data = land, fill = map_land_color, color = NA)
if (lgm_has_polygon) {
  map_plot <- map_plot + geom_sf(
    data = lgm, fill = lgm_fill_color, color = NA, alpha = lgm_fill_alpha
  )
} else {
  map_plot <- map_plot + geom_sf(
    data = lgm, fill = NA, color = lgm_line_color,
    linewidth = lgm_line_width, alpha = lgm_line_alpha
  )
}
# Redraw present-day national and coastal boundaries above the LGM layer.
map_plot <- map_plot + geom_sf(
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

# Population-derived halos replace broad, manually asserted refuge/contact
# circles. Their transparent interiors preserve the ancestry pies and LGM map.
for (signature_class in signature_levels) {
  z <- map_population_halos[
    map_population_halos$signature == signature_class, , drop = FALSE
  ]
  if (!nrow(z)) next
  map_plot <- map_plot + geom_sf(
    data = z, fill = NA, color = outline_colors[signature_class],
    linewidth = population_halo_linewidth, alpha = population_halo_alpha
  )
}

# Invisible points create an outline-color legend; the pie borders themselves
# are drawn in separate fixed-color layers below.
map_plot <- map_plot + geom_point(
  data = samples,
  aes(x = x_plot, y = y_plot, color = signature),
  alpha = 0, size = 0.01
)

pie_draw_order <- c("not_selected", "not_applicable", signature_levels)
for (signature_class in pie_draw_order) {
  z <- samples[samples$signature == signature_class, , drop = FALSE]
  if (!nrow(z)) next
  radius_multiplier <- if (signature_class == "not_selected") {
    unassigned_pie_scale
  } else {
    1
  }
  layer_alpha <- if (signature_class == "not_selected") {
    unassigned_pie_alpha
  } else {
    0.98
  }
  z$pie_radius_m <- pie_radius_km * radius_multiplier * 1000
  map_plot <- map_plot + scatterpie::geom_scatterpie(
    data = z,
    aes(x = x_plot, y = y_plot, r = pie_radius_m),
    cols = q_columns,
    color = outline_colors[signature_class], linewidth = 0.38,
    alpha = layer_alpha
  )
}

# Draw leaders above the pies so they remain visible when the jitter distance
# is smaller than the pie radius. The point at the line origin is the true
# collection coordinate; the other endpoint is the displayed pie center.
if (show_leader_lines && nrow(map_leaders)) {
  map_plot <- map_plot +
    geom_segment(
      data = map_leaders,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE, color = leader_line_color,
      linewidth = leader_line_width, alpha = leader_line_alpha,
      lineend = "round"
    ) +
    geom_point(
      data = map_leaders, aes(x = x, y = y), inherit.aes = FALSE,
      color = leader_line_color, size = coordinate_point_size,
      alpha = leader_line_alpha
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
    fill = guide_legend(
      nrow = 2, byrow = TRUE, order = 3,
      override.aes = list(alpha = 1, color = NA, linewidth = 0)
    ),
    color = guide_legend(
      nrow = 1, byrow = TRUE, order = 2,
      override.aes = list(alpha = 1, size = 2.8)
    )
  ) +
  coord_sf(
    crs = st_crs(3995),
    xlim = map_xlim, ylim = map_ylim,
    expand = FALSE
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal(base_size = base_size) +
  theme(
    panel.background = element_rect(fill = map_sea_color, color = NA),
    panel.grid.major = element_line(color = "#D2D2D2", linewidth = 0.22),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.32),
    plot.background = element_rect(fill = map_sea_color, color = NA),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(size = axis_title_size, color = "black"),
    axis.text = element_text(size = axis_text_size, color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.text = element_text(size = base_size * 0.82),
    legend.title = element_text(size = base_size * 0.92, face = "bold"),
    legend.key = element_blank(),
    legend.key.height = grid::unit(2.8, "mm"),
    legend.key.width = grid::unit(4.2, "mm"),
    plot.margin = margin(4, 1, 4, 4)
  )

# A second map isolates Beringia in a wide Lambert azimuthal equal-area
# projection. Coordinates and deterministic within-locality offsets are
# recalculated in this projection rather than reusing the polar-map offsets.
beringia_crs <- st_crs(sprintf(
  paste0(
    "+proj=laea +lat_0=%s +lon_0=%s +datum=WGS84 ",
    "+units=m +no_defs"
  ),
  beringia_center_lat, beringia_center_lon
))
if (!is.finite(beringia_half_width_km) || beringia_half_width_km <= 0) {
  stop("BERINGIA_HALF_WIDTH_KM must be a positive number")
}
if (!is.finite(beringia_half_width_x_km) || beringia_half_width_x_km <= 0) {
  stop("BERINGIA_HALF_WIDTH_X_KM must be a positive number")
}
if (!is.finite(beringia_half_height_km) || beringia_half_height_km <= 0) {
  stop("BERINGIA_HALF_HEIGHT_KM must be a positive number")
}
if (!is.finite(beringia_jitter_radius_km) || beringia_jitter_radius_km < 0) {
  stop("BERINGIA_JITTER_RADIUS_KM must be a non-negative number")
}
if (!is.finite(beringia_pie_radius_km) || beringia_pie_radius_km <= 0) {
  stop("BERINGIA_PIE_RADIUS_KM must be a positive number")
}
if (is.na(beringia_legend_ncol) || beringia_legend_ncol < 1) {
  stop("BERINGIA_LEGEND_NCOL must be a positive integer")
}

beringia_samples_sf <- st_as_sf(
  samples, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE
)
beringia_samples_sf <- st_transform(beringia_samples_sf, beringia_crs)
beringia_xy <- st_coordinates(beringia_samples_sf)
samples$beringia_origin_x <- beringia_xy[, 1]
samples$beringia_origin_y <- beringia_xy[, 2]
samples$beringia_x <- samples$beringia_origin_x
samples$beringia_y <- samples$beringia_origin_y

beringia_group <- trimws(as.character(samples$geoID))
missing_beringia_group <- is.na(beringia_group) | !nzchar(beringia_group)
beringia_group[missing_beringia_group] <- paste0(
  "sample_", samples$sampleID[missing_beringia_group]
)
golden_angle <- pi * (3 - sqrt(5))
for (idx in split(seq_len(nrow(samples)), beringia_group)) {
  idx <- idx[order(samples$sampleID[idx])]
  n <- length(idx)
  if (n < 2 || beringia_jitter_radius_km == 0) next
  center_x <- mean(samples$beringia_x[idx])
  center_y <- mean(samples$beringia_y[idx])
  point_rank <- seq_len(n) - 1
  radius <- beringia_jitter_radius_km * 1000 * sqrt(point_rank / max(point_rank))
  angle <- point_rank * golden_angle
  samples$beringia_x[idx] <- center_x + radius * cos(angle)
  samples$beringia_y[idx] <- center_y + radius * sin(angle)
}

beringia_leaders <- build_leader_segments(
  samples,
  "beringia_origin_x", "beringia_origin_y", "beringia_x", "beringia_y"
)

beringia_land <- st_transform(land, beringia_crs)
beringia_lgm <- st_transform(lgm, beringia_crs)
beringia_population_halos <- build_population_halos(
  samples, "beringia_x", "beringia_y", beringia_crs,
  beringia_pie_radius_km
)
beringia_population_halos <- dissolve_population_halos(
  beringia_population_halos
)
if (nrow(regions)) {
  beringia_region_points <- st_transform(region_points, beringia_crs)
  beringia_region_buffers <- st_buffer(
    beringia_region_points, dist = regions$radius_km * 1000
  )
}

beringia_plot <- ggplot() +
  geom_sf(data = beringia_land, fill = map_land_color, color = NA)
if (lgm_has_polygon) {
  beringia_plot <- beringia_plot + geom_sf(
    data = beringia_lgm,
    fill = lgm_fill_color, color = NA, alpha = lgm_fill_alpha
  )
} else {
  beringia_plot <- beringia_plot + geom_sf(
    data = beringia_lgm, fill = NA, color = lgm_line_color,
    linewidth = lgm_line_width, alpha = lgm_line_alpha
  )
}
beringia_plot <- beringia_plot + geom_sf(
  data = beringia_land, fill = NA,
  color = map_border_color, linewidth = 0.20
)

if (nrow(regions)) {
  beringia_plot <- beringia_plot +
    geom_sf(
      data = beringia_region_buffers,
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

for (signature_class in signature_levels) {
  z <- beringia_population_halos[
    beringia_population_halos$signature == signature_class, , drop = FALSE
  ]
  if (!nrow(z)) next
  beringia_plot <- beringia_plot + geom_sf(
    data = z, fill = NA, color = outline_colors[signature_class],
    linewidth = population_halo_linewidth, alpha = population_halo_alpha
  )
}

beringia_plot <- beringia_plot + geom_point(
  data = samples,
  aes(x = beringia_x, y = beringia_y, color = signature),
  alpha = 0, size = 0.01
)

for (signature_class in pie_draw_order) {
  z <- samples[samples$signature == signature_class, , drop = FALSE]
  if (!nrow(z)) next
  radius_multiplier <- if (signature_class == "not_selected") {
    unassigned_pie_scale
  } else {
    1
  }
  layer_alpha <- if (signature_class == "not_selected") {
    unassigned_pie_alpha
  } else {
    0.98
  }
  z$pie_radius_m <- beringia_pie_radius_km * radius_multiplier * 1000
  beringia_plot <- beringia_plot + scatterpie::geom_scatterpie(
    data = z,
    aes(x = beringia_x, y = beringia_y, r = pie_radius_m),
    cols = q_columns,
    color = outline_colors[signature_class], linewidth = 0.38,
    alpha = layer_alpha
  )
}

if (show_leader_lines && nrow(beringia_leaders)) {
  beringia_plot <- beringia_plot +
    geom_segment(
      data = beringia_leaders,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE, color = leader_line_color,
      linewidth = leader_line_width, alpha = leader_line_alpha,
      lineend = "round"
    ) +
    geom_point(
      data = beringia_leaders, aes(x = x, y = y), inherit.aes = FALSE,
      color = leader_line_color, size = coordinate_point_size,
      alpha = leader_line_alpha
    )
}

beringia_x_limit <- beringia_half_width_x_km * 1000
beringia_y_limit <- beringia_half_height_km * 1000
beringia_plot <- beringia_plot +
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
    fill = guide_legend(
      ncol = beringia_legend_ncol, order = 3,
      override.aes = list(alpha = 1, color = NA, linewidth = 0)
    ),
    color = guide_legend(
      ncol = min(beringia_legend_ncol, length(signature_levels)), order = 2,
      override.aes = list(alpha = 1, size = 2.8)
    )
  ) +
  coord_sf(
    crs = beringia_crs,
    xlim = c(-beringia_x_limit, beringia_x_limit),
    ylim = c(-beringia_y_limit, beringia_y_limit),
    expand = FALSE
  ) +
  labs(
    title = "Beringia: ancestry and diversity signatures",
    subtitle = paste0(
      "Lambert azimuthal equal-area projection centered at ",
      beringia_center_lat, " N, ", beringia_center_lon, " E"
    ),
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal(base_size = base_size) +
  theme(
    aspect.ratio = beringia_half_height_km / beringia_half_width_x_km,
    panel.background = element_rect(fill = map_sea_color, color = NA),
    panel.grid.major = element_line(color = "#D2D2D2", linewidth = 0.22),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.32),
    plot.background = element_rect(fill = map_sea_color, color = NA),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(size = axis_title_size, color = "black"),
    axis.text = element_text(size = axis_text_size, color = "black"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.text = element_text(size = base_size * 0.72),
    legend.title = element_text(size = base_size * 0.80, face = "bold"),
    legend.key = element_blank(),
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

# Descriptive per-individual heterozygosity. Population class comes from the
# same metadata assignments used for the map; unassigned samples are retained
# as a visual reference rather than treated as a demographic population.
individual_data <- merge(
  individual_heterozygosity[
    c("sampleID", "snp_heterozygosity", "call_rate")
  ],
  samples[c("sampleID", "pop_gen")],
  by = "sampleID", all = FALSE, sort = FALSE
)
individual_data <- individual_data[
  is.finite(individual_data$snp_heterozygosity), , drop = FALSE
]
if (!nrow(individual_data)) {
  stop("No mapped samples have finite individual SNP heterozygosity")
}
individual_data$sample_class <- factor(
  individual_data$pop_gen,
  levels = c("uniform", "admixed", "not_selected"),
  labels = c("Uniform", "Admixed", "Unassigned")
)
if (any(is.na(individual_data$sample_class))) {
  stop("Unexpected population class in individual heterozygosity data")
}
individual_class_colors <- c(
  Uniform = "#0072B2", Admixed = "#D55E00", Unassigned = "#999999"
)
class_counts <- table(individual_data$sample_class)
class_axis_labels <- setNames(
  paste0(names(class_counts), "\n(n = ", as.integer(class_counts), ")"),
  names(class_counts)
)

individual_heterozygosity_plot <- ggplot(
  individual_data,
  aes(x = sample_class, y = snp_heterozygosity, fill = sample_class)
) +
  geom_violin(
    trim = FALSE, scale = "width", color = "#4D4D4D",
    linewidth = 0.35, alpha = 0.42
  ) +
  geom_boxplot(
    width = 0.14, outlier.shape = NA, fill = "white",
    color = "#333333", linewidth = 0.38
  ) +
  geom_jitter(
    aes(color = sample_class), width = 0.10, height = 0,
    shape = 16, size = heterozygosity_point_size,
    alpha = 0.58, show.legend = FALSE
  ) +
  scale_fill_manual(values = individual_class_colors, guide = "none") +
  scale_color_manual(values = individual_class_colors, guide = "none") +
  scale_x_discrete(labels = class_axis_labels) +
  labs(
    x = "Population assignment",
    y = "Individual observed heterozygosity"
  ) +
  theme_classic(base_size = base_size + 1) +
  theme(
    axis.title = element_text(size = axis_title_size),
    axis.text.x = element_text(size = axis_text_size, color = "black"),
    axis.text.y = element_text(size = axis_text_size, color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Refugial-signature summary. Drummondii is excluded because it is not
# equivalent to the ingroup populations in the refugial comparison.
signature_plot_data <- signature_data[!is_drummondii, , drop = FALSE]
if (any(signature_plot_data$pi <= 0)) {
  stop("The refugial-signature log pi axis requires pi > 0")
}
if (any(signature_plot_data$private_allelic_richness < 0)) {
  stop("Private allelic richness cannot be negative")
}
x_log_range <- range(log10(signature_plot_data$pi))
y_sqrt_range <- range(sqrt(signature_plot_data$private_allelic_richness))
x_log_span <- diff(x_log_range)
y_sqrt_span <- diff(y_sqrt_range)
if (x_log_span == 0) x_log_span <- 1
if (y_sqrt_span == 0) y_sqrt_span <- 1
x_panel_range <- c(
  x_log_range[1] - signature_x_expand_left * x_log_span,
  x_log_range[2] + signature_x_expand_right * x_log_span
)
y_panel_range <- c(
  max(0, y_sqrt_range[1] - signature_y_expand_bottom * y_sqrt_span),
  y_sqrt_range[2] + signature_y_expand_top * y_sqrt_span
)
x_text <- 10^c(
  x_panel_range[1] + signature_text_x_inset * diff(x_panel_range),
  x_panel_range[2] - signature_text_x_inset * diff(x_panel_range)
)
y_text <- c(
  y_panel_range[2] - signature_text_y_inset * diff(y_panel_range),
  y_panel_range[1] + signature_text_y_inset * diff(y_panel_range)
)^2
quadrant_text <- data.frame(
  signature = factor(signature_levels, levels = signature_levels),
  x = c(sqrt(x_text[1] * x_text[2]), x_text[1], x_text[2]),
  y = c(y_text[1], y_text[2], y_text[2]),
  hjust = c(
    0.5, 0, signature_text_right_hjust
  ),
  vjust = c(1, 0, 0),
  label = c(
    "Candidate refugium",
    "Founder population",
    "Secondary contact"
  )
)

signature_plot <- ggplot(
  signature_plot_data,
  aes(x = pi, y = private_allelic_richness)
) +
  # Draw the vertical separator to the lower panel boundary, then mask its
  # upper portion so it separates only founder and secondary-contact space.
  geom_vline(
    xintercept = pi_cutoff, color = "#4D4D4D", linewidth = 0.45
  ) +
  annotate(
    "segment", x = pi_cutoff, xend = pi_cutoff,
    y = private_cutoff, yend = Inf,
    color = "white", linewidth = 0.95
  ) +
  geom_hline(yintercept = private_cutoff, color = "#4D4D4D", linewidth = 0.45) +
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
    x = expression("Nucleotide diversity (" * pi * "; log scale)"),
    y = "Private allelic richness (square-root scale)",
    caption = paste0(
      "Axes show original values on transformed scales;\n",
      "assignments use raw-value medians. ",
      "Private alleles are rarefied."
    )
  ) +
  scale_x_log10(
    labels = scales::label_scientific(digits = 2),
    expand = expansion(
      mult = c(signature_x_expand_left, signature_x_expand_right)
    )
  ) +
  scale_y_sqrt(
    labels = scales::label_number(accuracy = 0.001),
    expand = expansion(
      mult = c(signature_y_expand_bottom, signature_y_expand_top)
    )
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
    plot.caption = element_text(hjust = 0, size = 5.8),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    legend.title = element_text(size = 7.2),
    legend.text = element_text(size = 6.8),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  )

# Large supplementary version with every analyzed population labeled. The
# main composite uses the unlabeled form so its structure remains legible.
signature_label_data <- signature_plot_data
signature_label_data$display_population <- gsub(
  "_", " ", sub("^x_?", "", signature_label_data$population)
)
signature_si_plot <- signature_plot +
  geom_text_repel(
    data = signature_label_data,
    aes(label = display_population),
    size = signature_si_label_size, seed = 20261018,
    min.segment.length = 0, max.overlaps = Inf,
    box.padding = 0.38, point.padding = 0.20,
    force = 1.3, color = "black", show.legend = FALSE
  ) +
  theme(
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8)
  )

signature_pdf <- sub(
  "\\.pdf$", "_refugial_signatures.pdf", output_file, ignore.case = TRUE
)
ggsave(
  signature_pdf, signature_plot,
  width = signature_width, height = signature_height, units = "in"
)

signature_si_pdf <- sub(
  "\\.pdf$", "_refugial_signatures_labeled.pdf",
  output_file, ignore.case = TRUE
)
ggsave(
  signature_si_pdf, signature_si_plot,
  width = signature_si_width, height = signature_si_height, units = "in"
)

map_only <- sub("\\.pdf$", "_map.pdf", output_file, ignore.case = TRUE)
ggsave(map_only, map_plot, width = map_width, height = map_height, units = "in")

beringia_map <- sub(
  "\\.pdf$", "_beringia_map.pdf", output_file, ignore.case = TRUE
)
ggsave(
  beringia_map, beringia_plot,
  width = beringia_map_width, height = beringia_map_height, units = "in"
)

# The main figure contains exactly three panels: the circumpolar ancestry map,
# the compact refugial-signature summary, and descriptive individual Ho.
signature_panel <- signature_plot +
  theme(
    legend.position = "none",
    plot.caption = element_blank(),
    plot.margin = margin(5, 7, 5, 5)
  )
individual_heterozygosity_plot <- individual_heterozygosity_plot +
  theme(plot.margin = margin(5, 7, 5, 5))
side_panels <- signature_panel / individual_heterozygosity_plot +
  plot_layout(heights = c(1.08, 0.92))
composite <- map_plot | side_panels
composite <- composite + plot_layout(widths = c(2.05, 1.0)) +
  plot_annotation(tag_levels = "A")
ggsave(
  output_file, composite,
  width = figure_width, height = figure_height, units = "in"
)

message("Map-only figure: ", map_only)
message("Beringia equal-area map: ", beringia_map)
message("Composite figure: ", output_file)
message("Refugial-signature figure: ", signature_pdf)
message("Labeled SI refugial-signature figure: ", signature_si_pdf)
message("Refugial-signature assignments: ", signature_output)
message("Refugial-signature thresholds: ", signature_threshold_output)

