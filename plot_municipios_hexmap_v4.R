suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(scales); library(patchwork)
  library(mxmaps)
})
sf_use_s2(FALSE)

# v3.1 vs v4 comparison: v4 = pure unclipped hexes, largest-fit sizing,
# patch approximates the state (borders may deviate — accepted trade-off).

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

v31 <- readRDS("data/mxmunicipio_hex_sf_v3_1.rds")
v4  <- readRDS("data/mxmunicipio_hex_sf_v4.rds")
borders <- readRDS("data/mxstate_borders.rds")

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
mxmun_sf <- st_make_valid(st_sf(region = regions,
                                geometry = st_sfc(polys, crs = 4326)))
mxmun_sf$state_code <- substr(mxmun_sf$region, 1, 2)

state_sf <- mxmun_sf |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_make_valid()
state_sf_m <- st_transform(state_sf, 6372)

true_c <- mxmun_sf |>
  st_transform(6372) |>
  st_centroid(of_largest_polygon = TRUE) |>
  st_coordinates()

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
v31 <- add_disp(v31)
v4  <- add_disp(v4)

disp_stats <- function(hex, label) {
  cat(sprintf("%s — disp km: mean=%.2f median=%.2f p95=%.2f max=%.2f\n",
              label, mean(hex$disp_km), median(hex$disp_km),
              quantile(hex$disp_km, 0.95), max(hex$disp_km)))
  invisible(hex)
}
disp_stats(v31, "v3.1")
disp_stats(v4,  "v4  ")

coverage <- function(hex, label) {
  hex_m <- st_transform(hex, 6372)
  cov <- vapply(sort(unique(hex$state_code)), function(s) {
    cells <- st_union(hex_m$geometry[hex_m$state_code == s])
    st_poly <- state_sf_m$geometry[state_sf_m$state_code == s]
    as.numeric(st_area(st_intersection(cells, st_poly)) / st_area(st_poly))
  }, numeric(1))
  cat(sprintf("%s state coverage by own cells: min=%.1f%% mean=%.1f%%\n",
              label, 100 * min(cov), 100 * mean(cov)))
  invisible(cov)
}
coverage(v31, "v3.1")
coverage(v4,  "v4  ")

# Cell purity: how close each cell's area is to its state's ideal hex area
purity <- function(hex, label) {
  hex_m <- st_transform(hex, 6372)
  areas <- as.numeric(st_area(hex_m))
  ideal <- tapply(areas, hex_m$state_code, median)
  ratio <- areas / ideal[hex_m$state_code]
  cat(sprintf("%s cell area ratio to state median: p5=%.2f p95=%.2f (1.0 = even)\n",
              label, quantile(ratio, 0.05), quantile(ratio, 0.95)))
  invisible(ratio)
}
purity(v31, "v3.1")
purity(v4,  "v4  ")

# --- Plots --------------------------------------------------------------------
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
          plot.subtitle = element_text(color = "grey40", size = 9)) +
    labs(title = title, subtitle = subtitle)
}

hex_v31 <- v31 |> left_join(census, by = "region")
hex_v4  <- v4  |> left_join(census, by = "region")

p_full <- mk_map(hex_v31, "v3.1 — clipped to true borders",
                 "canonical candidate (mxmunicipio_hex_sf_v3_1.rds)") +
  mk_map(hex_v4, "v4 — pure hexes, even space, borders approximate",
         "proposed (mxmunicipio_hex_sf_v4.rds)") +
  plot_annotation(
    title = "Municipio hex cartogram — v3.1 vs v4",
    subtitle = "Share of population speaking an indigenous language, 2020 census",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.subtitle = element_text(color = "grey40", size = 11)))
ggsave("v3_1_vs_v4_full.png", p_full, width = 16, height = 8, dpi = 150, bg = "white")

zoom_plot <- function(state_code, fname, note) {
  b <- st_bbox(state_sf$geometry[state_sf$state_code == state_code])
  pad_x <- 0.25 * (b["xmax"] - b["xmin"])
  pad_y <- 0.25 * (b["ymax"] - b["ymin"])
  xlim <- c(b["xmin"] - pad_x, b["xmax"] + pad_x)
  ylim <- c(b["ymin"] - pad_y, b["ymax"] + pad_y)
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
  p <- mk_zoom(hex_v31 |> filter(state_code == state_code), "v3.1 — clipped") |
       mk_zoom(hex_v4  |> filter(state_code == state_code), "v4 — pure hexes") +
    plot_annotation(
      title = sprintf("%s zoom — %s", state_code, note),
      theme = theme(plot.title = element_text(face = "bold", size = 15)))
  ggsave(fname, p, width = 14, height = 8, dpi = 150, bg = "white")
  cat(sprintf("Wrote %s\n", fname))
}
zoom_plot("03", "v3_1_vs_v4_bcs.png", "Baja California Sur — absorption blobs vs pure hexes")
zoom_plot("12", "v3_1_vs_v4_guerrero.png", "Guerrero — dense/sparse border seam")

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
p_disp <- mk_disp(v31, "v3.1 — displacement") |
          mk_disp(v4,  "v4 — displacement") +
  plot_annotation(
    title = "Distance from true municipio centroid to hex centroid",
    theme = theme(plot.title = element_text(face = "bold", size = 15)))
ggsave("v3_1_vs_v4_displacement.png", p_disp,
       width = 16, height = 8, dpi = 150, bg = "white")

cat("Wrote: v3_1_vs_v4_full.png, v3_1_vs_v4_bcs.png,",
    "v3_1_vs_v4_guerrero.png, v3_1_vs_v4_displacement.png\n")
