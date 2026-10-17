#!/usr/bin/env Rscript

# Convert a polyline LGM boundary layer into fillable ice polygons.
#
# Example:
# MIN_AREA_KM2=1000 Rscript scripts/15_polygonize_lgm_boundary.R \
#   data/LGM_boundary/LGM_boundary.shp \
#   data/LGM_boundary/LGM_ice_polygons.gpkg

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(paste(
    "Usage: 15_polygonize_lgm_boundary.R",
    "INPUT_BOUNDARY.shp OUTPUT_POLYGONS.gpkg"
  ))
}

required <- "sf"
missing_packages <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(sf))

input_file <- args[1]
output_file <- args[2]
minimum_area_km2 <- as.numeric(Sys.getenv("MIN_AREA_KM2", "1000"))
minimum_latitude <- as.numeric(Sys.getenv("LGM_MIN_LATITUDE", "0"))
simplify_tolerance_m <- as.numeric(Sys.getenv("LGM_SIMPLIFY_TOLERANCE_M", "0"))
equal_area_crs <- as.integer(Sys.getenv("LGM_EQUAL_AREA_EPSG", "6931"))

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
if (!is.finite(minimum_area_km2) || minimum_area_km2 < 0) {
  stop("MIN_AREA_KM2 must be a non-negative number")
}
if (!is.finite(minimum_latitude) || minimum_latitude < -90 ||
    minimum_latitude > 90) {
  stop("LGM_MIN_LATITUDE must be between -90 and 90")
}
if (!is.finite(simplify_tolerance_m) || simplify_tolerance_m < 0) {
  stop("LGM_SIMPLIFY_TOLERANCE_M must be a non-negative number")
}
if (is.na(equal_area_crs) || equal_area_crs <= 0) {
  stop("LGM_EQUAL_AREA_EPSG must be a valid EPSG integer")
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

boundary <- st_read(input_file, quiet = TRUE)
if (is.na(st_crs(boundary))) {
  stop("The input LGM boundary has no coordinate reference system")
}
if (!nrow(boundary)) stop("The input LGM boundary contains no features")

# Remove features wholly south of the requested latitude before applying the
# Arctic equal-area projection. This also avoids projection failures caused by
# unrelated Southern Hemisphere features in global source layers.
boundary_wgs84 <- st_transform(boundary, 4326)
feature_ymax <- vapply(
  st_geometry(boundary_wgs84),
  function(g) unname(st_bbox(g)["ymax"]),
  numeric(1)
)
keep_north <- is.finite(feature_ymax) & feature_ymax >= minimum_latitude
boundary_wgs84 <- boundary_wgs84[keep_north, , drop = FALSE]
if (!nrow(boundary_wgs84)) {
  stop("No LGM boundary features remain north of LGM_MIN_LATITUDE")
}

boundary_wgs84 <- st_make_valid(boundary_wgs84)
input_types <- sort(unique(as.character(st_geometry_type(boundary_wgs84))))
if (!any(input_types %in% c("LINESTRING", "MULTILINESTRING",
                            "GEOMETRYCOLLECTION"))) {
  stop(
    "Expected line geometry, but retained features are: ",
    paste(input_types, collapse = ", ")
  )
}

if (all(input_types %in% c("LINESTRING", "MULTILINESTRING"))) {
  linework <- boundary_wgs84
} else {
  linework <- suppressWarnings(
    st_collection_extract(boundary_wgs84, "LINESTRING", warn = FALSE)
  )
}
linework <- linework[!st_is_empty(linework), , drop = FALSE]
if (!nrow(linework)) stop("No usable linework could be extracted")

# EPSG:6931 (by default) is a Northern Hemisphere equal-area projection. It
# makes the polygon filter an actual square-kilometer threshold rather than a
# threshold distorted by latitude or Web Mercator coordinates.
linework_equal_area <- st_transform(linework, equal_area_crs)
merged_linework <- st_union(st_geometry(linework_equal_area))
noded_linework <- st_node(merged_linework)
polygonized <- st_polygonize(noded_linework)
polygon_geometry <- suppressWarnings(
  st_collection_extract(polygonized, "POLYGON", warn = FALSE)
)
if (!length(polygon_geometry)) {
  stop(paste(
    "Polygonization produced no closed polygons.",
    "The source line endpoints may not form closed rings."
  ))
}

polygons <- st_sf(geometry = polygon_geometry)
polygons <- st_make_valid(polygons)
polygons <- polygons[!st_is_empty(polygons), , drop = FALSE]
area_before_km2 <- as.numeric(st_area(polygons)) / 1e6
keep_area <- is.finite(area_before_km2) & area_before_km2 >= minimum_area_km2
retained <- polygons[keep_area, , drop = FALSE]
if (!nrow(retained)) {
  stop(
    "No polygons meet MIN_AREA_KM2=", minimum_area_km2,
    ". Rerun with a smaller threshold."
  )
}

# Dissolve touching and overlapping pieces after filtering slivers, then cast
# back to individual polygon parts so the output retains useful area fields.
dissolved <- st_make_valid(st_union(st_geometry(retained)))
retained_geometry <- suppressWarnings(st_cast(dissolved, "POLYGON"))
retained <- st_sf(geometry = retained_geometry)
retained <- retained[!st_is_empty(retained), , drop = FALSE]

if (simplify_tolerance_m > 0) {
  retained$geometry <- st_simplify(
    retained$geometry,
    dTolerance = simplify_tolerance_m,
    preserveTopology = TRUE
  )
  retained <- st_make_valid(retained)
  retained <- retained[!st_is_empty(retained), , drop = FALSE]
}

retained$area_km2 <- as.numeric(st_area(retained)) / 1e6
retained <- retained[
  is.finite(retained$area_km2) & retained$area_km2 >= minimum_area_km2,
  , drop = FALSE
]
retained <- retained[order(retained$area_km2, decreasing = TRUE), , drop = FALSE]
retained$polygon_id <- seq_len(nrow(retained))
retained$source <- basename(input_file)
retained$min_area_km2 <- minimum_area_km2
retained$simplify_m <- simplify_tolerance_m
retained <- retained[c(
  "polygon_id", "area_km2", "source", "min_area_km2",
  "simplify_m", "geometry"
)]

# WGS84 makes the derived file portable; the plotting script reprojects it to
# the polar and Beringia map coordinate systems.
output <- st_transform(retained, 4326)
st_write(output, output_file, delete_layer = TRUE, quiet = TRUE)

summary_file <- paste0(tools::file_path_sans_ext(output_file), "_summary.tsv")
summary_table <- data.frame(
  input_file = normalizePath(input_file),
  output_file = normalizePath(output_file),
  retained_input_features = nrow(boundary_wgs84),
  candidate_polygons = length(area_before_km2),
  polygons_below_threshold = sum(!keep_area),
  output_polygons = nrow(output),
  minimum_area_km2 = minimum_area_km2,
  simplify_tolerance_m = simplify_tolerance_m,
  output_area_km2 = sum(output$area_km2),
  stringsAsFactors = FALSE
)
write.table(
  summary_table, summary_file,
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Input geometry types: ", paste(input_types, collapse = ", "))
message("Candidate polygons: ", length(area_before_km2))
message("Removed below ", minimum_area_km2, " km2: ", sum(!keep_area))
message("Output polygons: ", nrow(output))
message("Polygon layer: ", output_file)
message("Summary: ", summary_file)
