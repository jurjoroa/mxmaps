suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(scales); library(patchwork)
  library(mxmaps)
})
sf_use_s2(FALSE)

# v3.1 vs v5 comparison: v5 = v3.1 + area-balanced surplus absorption.
# Same containment (cells inside state boundaries, full tessellation),
# more even cell sizes inside each state.

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

v31 <- readRDS("data/mxmunicipio_hex_sf_v3_1.rds")
v5  <- readRDS("data/mxmunicipio_hex_sf_v5.rds")
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
v5  <- add_disp(v5)

# Area-ratio table (cell area vs state ideal share) — the evenness metric
ratio_table <- function(hex, label) {
  hex_m <- st_transform(hex, 6372)
  areas <- as.numeric(st_area(hex_m))
  share <- tapply(areas, hex_m$state_code, function(a) sum(a) / length(a))
  data.frame(version = label,
             state_code = hex_m$state_code,
             ratio = areas / share[hex_m$state_code])
}
rt31 <- ratio_table(v31, "v3.1 clipped borders")
rt5  <- ratio_table(v5,  "v5 balanced space")
rt <- rbind(rt31, rt5)
cat("Evenness (cell area / state ideal share):\n")
cat(sprintf("v3.1: p5=%.2f p25=%.2f p75=%.2f p95=%.2f\n",
            quantile(rt31$ratio, .05), quantile(rt31$ratio, .25),
            quantile(rt31$ratio, .75), quantile(rt31$ratio, .95)))
cat(sprintf("v5  : p5=%.2f p25=%.2f p75=%.2f p95=%.2f\n",
            quantile(rt5$ratio, .05), quantile(rt5$ratio, .25),
            quantile(rt5$ratio, .75), quantile(rt5$ratio, .95)))

disp_stats <- function(hex, label) {
  cat(sprintf("%s — disp km: mean=%.2f median=%.2f p95=%.2f max=%.2f\n",
              label, mean(hex$disp_km), median(hex$disp_km),
              quantile(hex$disp_km, 0.95), max(hex$disp_km)))
  invisible(hex)
}
disp_stats(v31, "v3.1")
disp_stats(v5,  "v5  ")

# --- Choropleth plots ----------------------------------------------------------
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
hex_v5  <- v5  |> left_join(census, by = "region")

p_full <- mk_map(hex_v31, "v3.1 — nearest absorption",
                 "mxmunicipio_hex_sf_v3_1.rds") +
  mk_map(hex_v5, "v5 — space balanced to the smallest cells",
         "mxmunicipio_hex_sf_v5.rds") +
  plot_annotation(
    title = "Municipio hex cartogram — v3.1 vs v5",
    subtitle = "Share of population speaking an indigenous language, 2020 census",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.subtitle = element_text(color = "grey40", size = 11)))
ggsave("v3_1_vs_v5_full.png", p_full, width = 16, height = 8, dpi = 150, bg = "white")

# Evenness distribution plot — the headline metric
p_even <- ggplot(rt, aes(x = version, y = pmin(ratio, 3), fill = version)) +
  geom_hline(yintercept = 1, color = "grey30", linetype = "dashed") +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("grey70", "#4DAF4A")) +
  coord_cartesian(ylim = c(0, 3)) +
  labs(title = "Cell size evenness inside each state",
       subtitle = "Cell area relative to its state's ideal equal share (1.0 = perfect; capped at 3)",
       x = NULL, y = "area ratio") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "grey40", size = 10))
ggsave("v3_1_vs_v5_evenness.png", p_even, width = 9, height = 7, dpi = 150, bg = "white")

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
  p <- mk_zoom(hex_v31 |> filter(state_code == state_code), "v3.1 — nearest") |
       mk_zoom(hex_v5  |> filter(state_code == state_code), "v5 — balanced") +
    plot_annotation(
      title = sprintf("%s zoom — %s", state_code, note),
      theme = theme(plot.title = element_text(face = "bold", size = 15)))
  ggsave(fname, p, width = 14, height = 8, dpi = 150, bg = "white")
  cat(sprintf("Wrote %s\n", fname))
}
zoom_plot("03", "v3_1_vs_v5_bcs.png", "Baja California Sur — sparse absorption")
zoom_plot("12", "v3_1_vs_v5_guerrero.png", "Guerrero — dense/sparse seam")

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
          mk_disp(v5,  "v5 — displacement") +
  plot_annotation(
    title = "Distance from true municipio centroid to hex centroid",
    theme = theme(plot.title = element_text(face = "bold", size = 15)))
ggsave("v3_1_vs_v5_displacement.png", p_disp,
       width = 16, height = 8, dpi = 150, bg = "white")

cat("Wrote: v3_1_vs_v5_full.png, v3_1_vs_v5_evenness.png,",
    "v3_1_vs_v5_bcs.png, v3_1_vs_v5_guerrero.png, v3_1_vs_v5_displacement.png\n")
