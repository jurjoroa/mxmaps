suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(scales); library(patchwork)
  library(mxmaps)
})
sf_use_s2(FALSE)

# v3 vs v3.1 comparison: same data, same styling, four views.
#   1. Country-level choropleth (indigenous-language share)
#   2. Oaxaca zoom (largest state, 570 municipios — assignment quality)
#   3. BC Sur zoom (sparse state, 5-6 municipios — blob/absorption test)
#   4. Displacement map (km from true municipio centroid to hex centroid)
# Plus a quantitative table: displacement stats + tessellation coverage.

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

v3   <- readRDS("data/mxmunicipio_hex_sf_v3.rds")
v31  <- readRDS("data/mxmunicipio_hex_sf_v3_1.rds")
borders <- readRDS("data/mxstate_borders.rds")

# --- True municipio centroids (metric CRS for km distances) -----------------
data("mxmunicipio.map")
df <- mxmunicipio.map
regions <- unique(df$region)
polys <- lapply(regions, function(r) {
  sub <- df[df$region == r, ]
  groups <- unique(sub$group)
  rings <- lapply(groups, function(g) {
    pts <- sub[sub$group == g, ]
    pts <- pts[order(pts$order), ]
    coords <- cbind(pts$long, pts$lat)
    if (nrow(coords) < 4) return(NULL)
    if (!all(coords[1,] == coords[nrow(coords),])) coords <- rbind(coords, coords[1,])
    list(coords = coords, hole = any(pts$hole))
  })
  rings <- Filter(Negate(is.null), rings)
  if (length(rings) == 0) return(st_polygon())
  exterior <- lapply(rings[!sapply(rings, `[[`, "hole")], `[[`, "coords")
  holes    <- lapply(rings[sapply(rings,  `[[`, "hole")], `[[`, "coords")
  st_polygon(c(exterior, holes))
})
mxmun_sf <- st_sf(region = regions, geometry = st_sfc(polys, crs = 4326))
mxmun_sf$state_code <- substr(mxmun_sf$region, 1, 2)
mxmun_sf <- st_make_valid(mxmun_sf)

true_c <- mxmun_sf |>
  st_transform(6372) |>
  st_centroid(of_largest_polygon = TRUE) |>
  st_coordinates()

# Displacement (km): true centroid -> final hex centroid, matched by region
add_disp <- function(hex) {
  hex_c <- hex |>
    st_transform(6372) |>
    st_centroid(of_largest_polygon = TRUE) |>
    st_coordinates()
  idx <- match(hex$region, mxmun_sf$region)
  hex$disp_km <- sqrt((hex_c[, 1] - true_c[idx, 1])^2 +
                      (hex_c[, 2] - true_c[idx, 2])^2) / 1000
  hex
}
v3  <- add_disp(v3)
v31 <- add_disp(v31)

# --- Quantitative comparison -------------------------------------------------
disp_stats <- function(hex, label) {
  cat(sprintf("\n%s — displacement km: mean=%.2f median=%.2f p95=%.2f max=%.2f\n",
              label, mean(hex$disp_km), median(hex$disp_km),
              quantile(hex$disp_km, 0.95), max(hex$disp_km)))
  worst <- hex |> st_drop_geometry() |>
    group_by(state_code) |>
    summarise(mean_disp = mean(disp_km), .groups = "drop") |>
    arrange(desc(mean_disp)) |>
    head(5)
  cat("Worst 5 states by mean displacement:\n")
  print(as.data.frame(worst, row.names = NULL))
  invisible(hex)
}
disp_stats(v3,  "v3  ")
disp_stats(v31, "v3.1")

# Tessellation coverage: union of cells vs state polygon area
data("mxmunicipio.map")
state_sf <- mxmun_sf |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_make_valid()
state_sf_m <- st_transform(state_sf, 6372)

coverage <- function(hex, label) {
  hex_m <- st_transform(hex, 6372)
  cov <- vapply(sort(unique(hex$state_code)), function(s) {
    cells <- st_union(hex_m$geometry[hex_m$state_code == s])
    st_poly <- state_sf_m$geometry[state_sf_m$state_code == s]
    as.numeric(st_area(st_intersection(cells, st_poly)) / st_area(st_poly))
  }, numeric(1))
  cat(sprintf("%s tessellation coverage: min=%.1f%% mean=%.1f%%\n",
              label, 100 * min(cov), 100 * mean(cov)))
  invisible(cov)
}
cov3  <- coverage(v3,  "v3  ")
cov31 <- coverage(v31, "v3.1")

# --- Choropleth plots ---------------------------------------------------------
data("df_mxmunicipio_2020")
census <- df_mxmunicipio_2020 |>
  select(region, indigenous_language, pop) |>
  mutate(indigenous_share = indigenous_language / pop)

mk_map <- function(hex_data, title, subtitle) {
  ggplot() +
    geom_sf(data = hex_data, aes(fill = indigenous_share),
            color = "grey55", linewidth = 0.05) +
    scale_fill_distiller(palette = "YlGnBu", direction = 1,
                         labels = percent_format(accuracy = 1),
                         na.value = "grey90") +
    geom_sf(data = borders, fill = NA, color = "black", linewidth = 0.35) +
    coord_sf() +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(color = "grey40", size = 9),
          legend.position = "right") +
    labs(title = title, subtitle = subtitle)
}

hex_v3  <- v3  |> left_join(census, by = "region")
hex_v31 <- v31 |> left_join(census, by = "region")

p_full <- mk_map(hex_v3, "v3 — greedy outermost-first assignment",
                 "canonical (data/mxmunicipio_hex_sf_v3.rds)") +
  mk_map(hex_v31, "v3.1 — Hungarian assignment + balanced absorption",
         "proposed (data/mxmunicipio_hex_sf_v3_1.rds)") +
  plot_annotation(
    title = "Municipio hex cartogram — v3 vs v3.1",
    subtitle = "Share of population speaking an indigenous language, 2020 census",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.subtitle = element_text(color = "grey40", size = 11)))
ggsave("v3_vs_v3_1_full.png", p_full, width = 16, height = 8, dpi = 150, bg = "white")

# Zoom helper: side-by-side for one state
zoom_plot <- function(state_code, fname, note) {
  b <- st_bbox(state_sf$geometry[state_sf$state_code == state_code])
  xlim <- c(b["xmin"] - 0.15 * (b["xmax"] - b["xmin"]),
            b["xmax"] + 0.15 * (b["xmax"] - b["xmin"]))
  ylim <- c(b["ymin"] - 0.15 * (b["ymax"] - b["ymin"]),
            b["ymax"] + 0.15 * (b["ymax"] - b["ymin"]))
  mk_zoom <- function(hex_data, title) {
    ggplot() +
      geom_sf(data = hex_data, aes(fill = indigenous_share),
              color = "grey45", linewidth = 0.12) +
      scale_fill_distiller(palette = "YlGnBu", direction = 1,
                           labels = percent_format(accuracy = 1),
                           na.value = "grey90") +
      geom_sf(data = borders, fill = NA, color = "black", linewidth = 0.45) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      theme_void(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 13),
            legend.position = "none") +
      labs(title = title)
  }
  p <- mk_zoom(hex_v3 |> filter(state_code == state_code),
               "v3 — greedy") |
       mk_zoom(hex_v31 |> filter(state_code == state_code),
               "v3.1 — Hungarian + balanced") +
    plot_annotation(
      title = sprintf("%s zoom — %s", state_code, note),
      theme = theme(plot.title = element_text(face = "bold", size = 15)))
  ggsave(fname, p, width = 14, height = 8, dpi = 150, bg = "white")
  cat(sprintf("Wrote %s\n", fname))
}
zoom_plot("20", "v3_vs_v3_1_oaxaca.png", "Oaxaca, 570 municipios")
zoom_plot("03", "v3_vs_v3_1_bcs.png", "Baja California Sur, sparse — absorption blobs")

# --- Displacement map ---------------------------------------------------------
mk_disp <- function(hex_data, title) {
  ggplot() +
    geom_sf(data = hex_data, aes(fill = disp_km),
            color = "grey55", linewidth = 0.04) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1,
                         name = "km", na.value = "grey90") +
    geom_sf(data = borders, fill = NA, color = "black", linewidth = 0.3) +
    coord_sf() +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13)) +
    labs(title = title)
}
p_disp <- mk_disp(v3, "v3 — displacement (greedy)") |
          mk_disp(v31, "v3.1 — displacement (Hungarian)") +
  plot_annotation(
    title = "Distance from true municipio centroid to hex centroid",
    theme = theme(plot.title = element_text(face = "bold", size = 15)))
ggsave("v3_vs_v3_1_displacement.png", p_disp,
       width = 16, height = 8, dpi = 150, bg = "white")

cat("Wrote: v3_vs_v3_1_full.png, v3_vs_v3_1_oaxaca.png,",
    "v3_vs_v3_1_bcs.png, v3_vs_v3_1_displacement.png\n")
