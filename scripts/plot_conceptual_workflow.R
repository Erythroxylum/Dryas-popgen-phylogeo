#!/usr/bin/env Rscript

# ==============================================================================
# Conceptual workflow: resolving evolutionary lineages in Dryas
#
# Panel A: Genome-wide lineage discovery
#          ADMIXTURE -> s170 lineages -> BPP sliding-window population trees
#
# Panel B: Demographic model construction
#          Pairwise BPP IM screens -> hierarchical MSC-M models
#
# Panel C: Evolutionary interpretation
#          tau -> divergence time
#          theta + tau + W -> aab/abb simulations -> gdi
#
# IMPORTANT:
# Map coordinates, ADMIXTURE proportions, and phylogenies in this figure are
# SCHEMATIC. They are intended to communicate the analysis workflow and should
# not be interpreted as empirical results.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------------------

pkgs <- c(
  "ggplot2",
  "patchwork",
  "sf",
  "rnaturalearth",
  "rnaturalearthdata",
  "dplyr",
  "tidyr",
  "grid",
  "scales"
)

missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Install missing packages first:\n",
    paste0(
      "install.packages(c(",
      paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
      "))"
    )
  )
}

library(ggplot2)
library(patchwork)
library(sf)
library(rnaturalearth)
library(dplyr)
library(tidyr)
library(grid)
library(scales)


# ------------------------------------------------------------------------------
# 1. Output
# ------------------------------------------------------------------------------

outdir <- "figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

outfile_pdf <- file.path(outdir, "FIG_conceptual_Dryas_workflow.pdf")
outfile_png <- file.path(outdir, "FIG_conceptual_Dryas_workflow.png")


# ------------------------------------------------------------------------------
# 2. Theme and colors
# ------------------------------------------------------------------------------

theme_concept <- theme_void(base_size = 10) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0,
      margin = margin(b = 7)
    ),
    plot.subtitle = element_text(
      size = 9,
      hjust = 0,
      margin = margin(b = 5)
    ),
    plot.margin = margin(5, 5, 5, 5)
  )

# Conceptual lineage palette.
# These colors are illustrative and can later be replaced with manuscript colors.
lineage_cols <- c(
  "#D73027",
  "#FC8D59",
  "#FEE08B",
  "#91CF60",
  "#1A9850",
  "#66C2A5",
  "#3288BD",
  "#5E4FA2",
  "#984EA3"
)

names(lineage_cols) <- paste0("L", 1:9)

grey_dark  <- "grey20"
grey_mid   <- "grey55"
grey_light <- "grey85"
green_mig  <- "#55B748"


# ------------------------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------------------------

arrow_segment <- function(x, y, xend, yend,
                          linewidth = 0.6,
                          colour = grey_dark,
                          type = "closed") {

  geom_segment(
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    linewidth = linewidth,
    colour = colour,
    arrow = arrow(
      length = unit(0.12, "inches"),
      type = type
    )
  )
}


# Simple rectangular node/label
box_label <- function(label,
                      x,
                      y,
                      width = 0.28,
                      height = 0.10,
                      fill = "white",
                      colour = "black",
                      fontsize = 3.4,
                      linewidth = 0.5) {

  list(
    annotate(
      "rect",
      xmin = x - width / 2,
      xmax = x + width / 2,
      ymin = y - height / 2,
      ymax = y + height / 2,
      fill = fill,
      colour = colour,
      linewidth = linewidth
    ),
    annotate(
      "text",
      x = x,
      y = y,
      label = label,
      size = fontsize,
      fontface = "bold"
    )
  )
}


# Draw a tiny rooted tree from manually supplied segments.
# This keeps the conceptual trees dependency-free.
tiny_tree <- function(x0 = 0,
                      y0 = 0,
                      sx = 1,
                      sy = 1,
                      topology = 1,
                      colour = "grey25",
                      linewidth = 0.45,
                      alpha = 1) {

  if (topology == 1) {
    seg <- data.frame(
      x = c(0.00, 0.28, 0.28, 0.55, 0.55, 0.78, 0.78),
      y = c(0.50, 0.50, 0.50, 0.28, 0.72, 0.17, 0.39),
      xend = c(0.28, 0.28, 0.55, 0.78, 0.78, 1.00, 1.00),
      yend = c(0.50, 0.28, 0.72, 0.17, 0.83, 0.17, 0.39)
    )
    # Add vertical connectors
    seg <- bind_rows(
      seg,
      data.frame(
        x = c(0.28, 0.55, 0.78),
        y = c(0.28, 0.28, 0.17),
        xend = c(0.28, 0.55, 0.78),
        yend = c(0.72, 0.72, 0.39)
      )
    )
  } else if (topology == 2) {
    seg <- data.frame(
      x = c(0.00, 0.30, 0.30, 0.55, 0.55, 0.80, 0.80),
      y = c(0.50, 0.50, 0.50, 0.65, 0.35, 0.18, 0.52),
      xend = c(0.30, 0.30, 0.55, 0.80, 0.80, 1.00, 1.00),
      yend = c(0.50, 0.25, 0.65, 0.82, 0.18, 0.18, 0.52)
    )
    seg <- bind_rows(
      seg,
      data.frame(
        x = c(0.30, 0.55, 0.80),
        y = c(0.25, 0.35, 0.18),
        xend = c(0.30, 0.55, 0.80),
        yend = c(0.65, 0.65, 0.52)
      )
    )
  } else {
    seg <- data.frame(
      x = c(0.00, 0.25, 0.25, 0.50, 0.50, 0.75, 0.75),
      y = c(0.50, 0.50, 0.50, 0.20, 0.80, 0.62, 0.90),
      xend = c(0.25, 0.25, 0.50, 0.75, 0.75, 1.00, 1.00),
      yend = c(0.50, 0.20, 0.80, 0.10, 0.62, 0.62, 0.90)
    )
    seg <- bind_rows(
      seg,
      data.frame(
        x = c(0.25, 0.50, 0.75),
        y = c(0.20, 0.20, 0.62),
        xend = c(0.25, 0.50, 0.75),
        yend = c(0.80, 0.80, 0.90)
      )
    )
  }

  seg <- seg |>
    mutate(
      x = x0 + x * sx,
      xend = x0 + xend * sx,
      y = y0 + y * sy,
      yend = y0 + yend * sy
    )

  geom_segment(
    data = seg,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    linewidth = linewidth,
    colour = colour,
    alpha = alpha,
    lineend = "round"
  )
}


# A small four-tip demographic tree for Panel B
demographic_tree <- function(x0,
                             y0,
                             width = 1,
                             height = 1,
                             migration = TRUE,
                             migration_y = 0.45) {

  tree <- data.frame(
    x = c(
      0.00, 0.28,
      0.28, 0.28,
      0.28, 0.55,
      0.55, 0.55,
      0.55, 0.82,
      0.55, 0.82
    ),
    y = c(
      0.50, 0.50,
      0.20, 0.20,
      0.80, 0.80,
      0.80, 0.80,
      0.62, 0.62,
      0.95, 0.95
    ),
    xend = c(
      0.28, 0.28,
      1.00, 0.28,
      0.55, 0.55,
      0.82, 0.82,
      1.00, 0.82,
      1.00, 0.82
    ),
    yend = c(
      0.50, 0.20,
      0.20, 0.80,
      0.80, 0.62,
      0.62, 0.95,
      0.62, 0.95,
      0.95, 0.62
    )
  )

  tree <- tree |>
    mutate(
      x = x0 + x * width,
      xend = x0 + xend * width,
      y = y0 + y * height,
      yend = y0 + yend * height
    )

  layers <- list(
    geom_segment(
      data = tree,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.65,
      colour = "black"
    )
  )

  if (migration) {
    layers <- c(
      layers,
      list(
        geom_segment(
          aes(
            x = x0 + 0.57 * width,
            y = y0 + migration_y * height,
            xend = x0 + 0.92 * width,
            yend = y0 + migration_y * height
          ),
          inherit.aes = FALSE,
          linewidth = 0.65,
          colour = green_mig,
          arrow = arrow(
            length = unit(0.09, "inches"),
            type = "closed"
          )
        )
      )
    )
  }

  layers
}


# ------------------------------------------------------------------------------
# 4. PANEL A1: Arctic map with schematic samples
# ------------------------------------------------------------------------------

world <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
)

# Keep Northern Hemisphere; projection will crop later.
world_north <- world |>
  st_crop(
    xmin = -180,
    xmax = 180,
    ymin = 25,
    ymax = 90
  )

# Lambert azimuthal equal-area centered on North Pole.
crs_arctic <- paste0(
  "+proj=laea +lat_0=90 +lon_0=0 ",
  "+datum=WGS84 +units=km +no_defs"
)

world_arctic <- st_transform(world_north, crs_arctic)


# Schematic geographic sampling centers.
# THESE ARE NOT EMPIRICAL SAMPLE LOCATIONS.
set.seed(12)

centers <- tibble(
  region = c(
    "Alaska", "W Canada", "Greenland",
    "Iceland", "Scandinavia", "Alps",
    "Siberia", "Far East", "Arctic Canada"
  ),
  lon = c(-150, -120, -42, -19, 18, 10, 95, 145, -90),
  lat = c(67, 58, 70, 65, 68, 46, 68, 63, 72),
  lineage = paste0("L", 1:9)
)

dummy_samples <- centers |>
  slice(rep(1:n(), each = 5)) |>
  mutate(
    lon = lon + rnorm(n(), 0, 5),
    lat = lat + rnorm(n(), 0, 2.0)
  )

dummy_sf <- st_as_sf(
  dummy_samples,
  coords = c("lon", "lat"),
  crs = 4326
) |>
  st_transform(crs_arctic)

dummy_xy <- bind_cols(
  dummy_samples,
  as.data.frame(st_coordinates(dummy_sf))
)


p_map <- ggplot() +
  geom_sf(
    data = world_arctic,
    fill = "grey96",
    colour = "grey72",
    linewidth = 0.18
  ) +
  geom_point(
    data = dummy_xy,
    aes(x = X, y = Y, colour = lineage),
    size = 1.6,
    alpha = 0.9
  ) +
  scale_colour_manual(values = lineage_cols, guide = "none") +
  coord_sf(
    xlim = c(-6000, 6000),
    ylim = c(-6000, 6000),
    expand = FALSE
  ) +
  labs(
    title = "Genome-wide sampling",
    subtitle = "384 individuals across the circumpolar range"
  ) +
  theme_concept +
  theme(
    panel.border = element_rect(
      colour = "grey60",
      fill = NA,
      linewidth = 0.4
    )
  )


# ------------------------------------------------------------------------------
# 5. PANEL A2: miniature ADMIXTURE plot
# ------------------------------------------------------------------------------

set.seed(5)

n_ind <- 42
K <- 9

admix <- matrix(0, nrow = n_ind, ncol = K)

groups <- rep(1:9, length.out = n_ind)

for (i in seq_len(n_ind)) {

  g <- groups[i]

  x <- rep(0.005, K)
  x[g] <- runif(1, 0.70, 0.96)

  # Add small admixture with neighboring lineage.
  ng <- ifelse(g == K, 1, g + 1)
  x[ng] <- runif(1, 0.02, 0.20)

  x <- x / sum(x)
  admix[i, ] <- x
}

admix_df <- as.data.frame(admix)
names(admix_df) <- paste0("L", 1:K)

admix_df$individual <- seq_len(n_ind)

admix_long <- admix_df |>
  pivot_longer(
    cols = starts_with("L"),
    names_to = "lineage",
    values_to = "ancestry"
  )

p_admix <- ggplot(
  admix_long,
  aes(
    x = individual,
    y = ancestry,
    fill = lineage
  )
) +
  geom_col(width = 1) +
  scale_fill_manual(values = lineage_cols, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "ADMIXTURE",
    subtitle = "Genome-wide ancestry defines s170 lineages"
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5),
    panel.border = element_rect(
      colour = "grey60",
      fill = NA,
      linewidth = 0.4
    ),
    plot.margin = margin(5, 5, 5, 5)
  )


# ------------------------------------------------------------------------------
# 6. PANEL A3: sliding-window BPP population trees
# ------------------------------------------------------------------------------

p_trees <- ggplot() +
  xlim(0, 5.3) +
  ylim(0, 4.0) +

  # faded background trees
  tiny_tree(
    x0 = 0.25, y0 = 2.2,
    sx = 2.0, sy = 1.25,
    topology = 1,
    colour = "grey75",
    linewidth = 0.55
  ) +
  tiny_tree(
    x0 = 0.55, y0 = 1.75,
    sx = 2.0, sy = 1.25,
    topology = 2,
    colour = "grey65",
    linewidth = 0.55
  ) +
  tiny_tree(
    x0 = 0.85, y0 = 1.30,
    sx = 2.0, sy = 1.25,
    topology = 3,
    colour = "grey55",
    linewidth = 0.55
  ) +

  # foreground tree
  tiny_tree(
    x0 = 1.15, y0 = 0.85,
    sx = 2.0, sy = 1.25,
    topology = 1,
    colour = "grey15",
    linewidth = 0.75
  ) +

  annotate(
    "text",
    x = 2.6,
    y = 3.75,
    label = "BPP sliding-window population trees",
    fontface = "bold",
    size = 3.7
  ) +
  annotate(
    "text",
    x = 2.6,
    y = 3.40,
    label = "28 populations across chromosomes/windows",
    size = 3.0
  ) +
  annotate(
    "text",
    x = 2.6,
    y = 0.42,
    label = "recurrent relationships + topological discordance",
    size = 3.0,
    colour = grey_mid
  ) +
  theme_void() +
  theme(
    panel.border = element_rect(
      colour = "grey60",
      fill = NA,
      linewidth = 0.4
    )
  )


# ------------------------------------------------------------------------------
# 7. Combine Panel A
# ------------------------------------------------------------------------------

panel_A <- (
  p_map |
    p_admix |
    p_trees
) +
  plot_layout(widths = c(1.15, 0.85, 1.0)) +
  plot_annotation(
    title = "A   Discover and define evolutionary lineages",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0
      )
    )
  )


# ------------------------------------------------------------------------------
# 8. PANEL B1: reduction to nine major lineages
# ------------------------------------------------------------------------------

lineage_nodes <- tibble(
  lineage = paste0("L", 1:9),
  x = rep(1:3, each = 3),
  y = rep(c(3, 2, 1), 3)
)

p_lineages <- ggplot() +
  geom_point(
    data = lineage_nodes,
    aes(x = x, y = y, fill = lineage),
    shape = 21,
    size = 8,
    stroke = 0.5
  ) +
  geom_text(
    data = lineage_nodes,
    aes(x = x, y = y, label = lineage),
    size = 3,
    fontface = "bold"
  ) +
  scale_fill_manual(values = lineage_cols, guide = "none") +
  xlim(0.4, 3.6) +
  ylim(0.45, 3.75) +
  annotate(
    "text",
    x = 2,
    y = 3.55,
    label = "Nine major lineages",
    fontface = "bold",
    size = 4
  ) +
  annotate(
    "text",
    x = 2,
    y = 0.62,
    label = "hierarchical reduction of population structure",
    size = 2.8,
    colour = grey_mid
  ) +
  theme_void() +
  theme(
    panel.border = element_rect(
      colour = "grey60",
      fill = NA,
      linewidth = 0.4
    )
  )


# ------------------------------------------------------------------------------
# 9. PANEL B2: pairwise BPP IM screening matrix
# ------------------------------------------------------------------------------

pair_grid <- expand_grid(
  source = paste0("L", 1:9),
  target = paste0("L", 1:9)
) |>
  mutate(
    i = match(source, paste0("L", 1:9)),
    j = match(target, paste0("L", 1:9)),
    signal = "none"
  )

# Schematic candidate signals only.
# NOT empirical migration estimates.
major_pairs <- tibble(
  i = c(1, 2, 4, 5, 7, 8),
  j = c(2, 1, 5, 4, 8, 7)
)

minor_pairs <- tibble(
  i = c(2, 3, 5, 6, 8),
  j = c(3, 2, 6, 5, 9)
)

pair_grid <- pair_grid |>
  mutate(
    signal = case_when(
      paste(i, j) %in% paste(major_pairs$i, major_pairs$j) ~ "candidate",
      paste(i, j) %in% paste(minor_pairs$i, minor_pairs$j) ~ "weak",
      i == j ~ "self",
      TRUE ~ "none"
    )
  )

matrix_cols <- c(
  self = "grey30",
  candidate = green_mig,
  weak = "grey65",
  none = "white"
)

p_pairwise <- ggplot(
  pair_grid,
  aes(x = j, y = 10 - i, fill = signal)
) +
  geom_tile(
    colour = "grey75",
    linewidth = 0.25
  ) +
  scale_fill_manual(values = matrix_cols, guide = "none") +
  scale_x_continuous(
    breaks = 1:9,
    labels = paste0("L", 1:9),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = 1:9,
    labels = rev(paste0("L", 1:9)),
    expand = c(0, 0)
  ) +
  coord_equal() +
  labs(
    title = "Pairwise BPP IM screening",
    subtitle = "Candidate migration edges among major lineages"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 6.5),
    axis.text.y = element_text(size = 6.5),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8),
    panel.border = element_rect(
      colour = "grey60",
      fill = NA,
      linewidth = 0.4
    )
  )


# ------------------------------------------------------------------------------
# 10. PANEL B3: three hierarchical MSC-M models
# ------------------------------------------------------------------------------

p_models <- ggplot() +
  xlim(0, 10) +
  ylim(0, 5.0) +

  # H4d
  demographic_tree(
    x0 = 0.5,
    y0 = 1.25,
    width = 2.2,
    height = 2.2,
    migration = TRUE,
    migration_y = 0.43
  ) +

  # AIH
  demographic_tree(
    x0 = 3.9,
    y0 = 1.25,
    width = 2.2,
    height = 2.2,
    migration = TRUE,
    migration_y = 0.67
  ) +

  # ajan-alas
  demographic_tree(
    x0 = 7.3,
    y0 = 1.25,
    width = 2.2,
    height = 2.2,
    migration = TRUE,
    migration_y = 0.55
  ) +

  annotate(
    "text",
    x = 1.6, y = 4.35,
    label = "H4d",
    fontface = "bold",
    size = 4.2
  ) +
  annotate(
    "text",
    x = 5.0, y = 4.35,
    label = "AIH",
    fontface = "bold",
    size = 4.2
  ) +
  annotate(
    "text",
    x = 8.4, y = 4.35,
    label = "ajan–alas",
    fontface = "bold",
    size = 4.2
  ) +

  annotate(
    "text",
    x = 5,
    y = 0.55,
    label = expression(
      paste("chromosome-specific MSC-M estimates:  ", theta, "   ", tau, "   W")
    ),
    size = 3.6
  ) +

  annotate(
    "text",
    x = 5,
    y = 4.82,
    label = "Hierarchical demographic models",
    fontface = "bold",
    size = 4.3
  ) +

  theme_void() +
  theme(
    panel.border = element_rect(
      colour = "grey60",
      fill = NA,
      linewidth = 0.4
    )
  )


# ------------------------------------------------------------------------------
# 11. Combine Panel B
# ------------------------------------------------------------------------------

panel_B <- (
  p_lineages |
    p_pairwise |
    p_models
) +
  plot_layout(widths = c(0.65, 0.9, 1.65)) +
  plot_annotation(
    title = "B   Build migration-aware demographic models",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0
      )
    )
  )


# ------------------------------------------------------------------------------
# 12. PANEL C1: divergence-time branch
# ------------------------------------------------------------------------------

p_time <- ggplot() +
  xlim(0, 10) +
  ylim(0, 10) +

  # central model box
  annotate(
    "rect",
    xmin = 0.6, xmax = 4.1,
    ymin = 6.0, ymax = 8.5,
    fill = "grey96",
    colour = "black",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 2.35, y = 7.65,
    label = "Fitted MSC-M",
    fontface = "bold",
    size = 4.2
  ) +
  annotate(
    "text",
    x = 2.35, y = 6.75,
    label = expression(theta ~ "   " ~ tau ~ "   " ~ W),
    size = 4.2
  ) +

  # tau branch
  geom_segment(
    aes(x = 4.2, y = 7.25, xend = 6.1, yend = 7.25),
    linewidth = 0.65,
    arrow = arrow(
      length = unit(0.12, "inches"),
      type = "closed"
    )
  ) +

  annotate(
    "rect",
    xmin = 6.25, xmax = 9.5,
    ymin = 6.0, ymax = 8.5,
    fill = "white",
    colour = "black",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 7.875, y = 7.65,
    label = expression(tau ~ "\u2192 divergence time"),
    fontface = "bold",
    size = 3.8
  ) +
  annotate(
    "text",
    x = 7.875, y = 6.75,
    label = "Pleistocene diversification",
    size = 3.2
  ) +

  # gdi branch
  geom_segment(
    aes(x = 2.35, y = 5.9, xend = 2.35, yend = 4.55),
    linewidth = 0.65,
    arrow = arrow(
      length = unit(0.12, "inches"),
      type = "closed"
    )
  ) +

  annotate(
    "rect",
    xmin = 0.6, xmax = 4.1,
    ymin = 2.0, ymax = 4.35,
    fill = "white",
    colour = "black",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 2.35, y = 3.70,
    label = "aab / abb gene-tree simulations",
    fontface = "bold",
    size = 3.4
  ) +
  annotate(
    "text",
    x = 2.35, y = 2.85,
    label = expression(theta ~ "+" ~ tau ~ "+" ~ W),
    size = 3.4
  ) +
  annotate(
    "text",
    x = 2.35, y = 2.35,
    label = "1,000,000 trees / chromosome / direction",
    size = 2.6,
    colour = grey_mid
  ) +

  geom_segment(
    aes(x = 4.2, y = 3.15, xend = 6.1, yend = 3.15),
    linewidth = 0.65,
    arrow = arrow(
      length = unit(0.12, "inches"),
      type = "closed"
    )
  ) +

  annotate(
    "rect",
    xmin = 6.25, xmax = 9.5,
    ymin = 2.0, ymax = 4.35,
    fill = "white",
    colour = "black",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 7.875, y = 3.55,
    label = "genealogical divergence index",
    fontface = "bold",
    size = 3.4
  ) +
  annotate(
    "text",
    x = 7.875, y = 2.75,
    label = "How exclusive has the lineage become?",
    size = 3.0
  ) +

  labs(
    title = "Demographic models provide two complementary outputs"
  ) +
  theme_concept +
  theme(
    panel.border = element_rect(
      colour = "grey60",
      fill = NA,
      linewidth = 0.4
    )
  )


# ------------------------------------------------------------------------------
# 13. PANEL C2: gdi continuum
# ------------------------------------------------------------------------------

gdi_examples <- tibble(
  comparison = c(
    "ajan–alas",
    "integrifolia\npopulations",
    "EU–Carpathian/\nMacedonian octopetala"
  ),
  x = c(0.11, 0.33, 0.43),
  y = c(1.1, 1.85, 1.1)
)

p_gdi <- ggplot() +
  xlim(0, 1) +
  ylim(0, 3.5) +

  # heuristic regions
  annotate(
    "rect",
    xmin = 0, xmax = 0.2,
    ymin = 0.15, ymax = 3.1,
    fill = "grey94",
    colour = NA
  ) +
  annotate(
    "rect",
    xmin = 0.2, xmax = 0.7,
    ymin = 0.15, ymax = 3.1,
    fill = "grey98",
    colour = NA
  ) +
  annotate(
    "rect",
    xmin = 0.7, xmax = 1,
    ymin = 0.15, ymax = 3.1,
    fill = "grey94",
    colour = NA
  ) +

  geom_segment(
    aes(x = 0, y = 0.75, xend = 1, yend = 0.75),
    linewidth = 0.75,
    colour = "black"
  ) +

  geom_segment(
    aes(x = 0.2, y = 0.55, xend = 0.2, yend = 0.95),
    linewidth = 0.45,
    linetype = "dashed",
    colour = grey_mid
  ) +
  geom_segment(
    aes(x = 0.7, y = 0.55, xend = 0.7, yend = 0.95),
    linewidth = 0.45,
    linetype = "dashed",
    colour = grey_mid
  ) +

  geom_point(
    data = gdi_examples,
    aes(x = x, y = y),
    size = 3.3,
    shape = 21,
    fill = "white",
    stroke = 0.8
  ) +
  geom_segment(
    data = gdi_examples,
    aes(x = x, y = 0.78, xend = x, yend = y - 0.16),
    linewidth = 0.35,
    colour = grey_mid
  ) +
  geom_text(
    data = gdi_examples,
    aes(x = x, y = y + 0.32, label = comparison),
    size = 2.8,
    lineheight = 0.9
  ) +

  annotate(
    "text",
    x = 0.10, y = 2.85,
    label = "weak genealogical\nexclusivity",
    size = 3.1,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.45, y = 2.85,
    label = "intermediate / ambiguous",
    size = 3.1,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.85, y = 2.85,
    label = "strong genealogical\nexclusivity",
    size = 3.1,
    fontface = "bold"
  ) +

  annotate(
    "text",
    x = 0.2, y = 0.33,
    label = "0.2",
    size = 2.8
  ) +
  annotate(
    "text",
    x = 0.7, y = 0.33,
    label = "0.7",
    size = 2.8
  ) +
  annotate(
    "text",
    x = 0.5, y = 0.10,
    label = "Kornai et al. heuristic reference regions — not species thresholds",
    size = 2.7,
    colour = grey_mid
  ) +

  labs(
    title = "Genealogical exclusivity under the fitted MSC-M history",
    subtitle = "Selected Dryas lineage comparisons span a continuum of gdi"
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12
    ),
    plot.subtitle = element_text(
      size = 8.5
    ),
    panel.border = element_rect(
      colour = "grey60",
      fill = NA,
      linewidth = 0.4
    ),
    plot.margin = margin(5, 5, 5, 5)
  )


# ------------------------------------------------------------------------------
# 14. Combine Panel C
# ------------------------------------------------------------------------------

panel_C <- (
  p_time |
    p_gdi
) +
  plot_layout(widths = c(1.05, 1.0)) +
  plot_annotation(
    title = "C   Estimate divergence history and quantify lineage exclusivity",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0
      )
    )
  )


# ------------------------------------------------------------------------------
# 15. Full figure
# ------------------------------------------------------------------------------

full_figure <- (
  panel_A /
    panel_B /
    panel_C
) +
  plot_layout(
    heights = c(1.00, 0.92, 0.90)
  ) +
  plot_annotation(
    title = "From genome-wide population structure to migration-aware evolutionary inference in Dryas",
    subtitle = paste(
      "Conceptual workflow linking ancestry, population-tree variation,",
      "demographic modeling, divergence-time inference, and genealogical exclusivity"
    ),
    caption = paste(
      "Map locations, ancestry proportions, phylogenies, and pairwise migration signals",
      "shown here are schematic. Empirical analyses were performed chromosome-by-chromosome."
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 18,
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = 11,
        hjust = 0
      ),
      plot.caption = element_text(
        size = 8,
        colour = "grey40",
        hjust = 0
      ),
      plot.margin = margin(10, 10, 10, 10)
    )
  )


# ------------------------------------------------------------------------------
# 16. Export
# ------------------------------------------------------------------------------

ggsave(
  outfile_pdf,
  full_figure,
  width = 15,
  height = 12,
  units = "in",
  device = cairo_pdf
)

ggsave(
  outfile_png,
  full_figure,
  width = 15,
  height = 12,
  units = "in",
  dpi = 400,
  bg = "white"
)

message("Wrote:")
message("  ", outfile_pdf)
message("  ", outfile_png)
