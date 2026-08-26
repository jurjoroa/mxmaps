# ---------------------------------------------------------------------------
# Purpose: render the v10 uniform-hexagon cartogram against true territory.
#          Everything is drawn in EPSG:6372 so the hexagons stay regular --
#          reprojecting to 4326 for display would shear them.
# Author:  Claude (Opus 5) for Jorge Roa
# Version: v10
# Date:    2026-08-25
# Usage:   Rscript plot_hexbin_v10.R [state_code ...]
# Output:  state_pages/hexbin_v10_national.png
#          state_pages/hexbin_v10_<code>_<abbr>.png
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(scales)
  library(patchwork); library(mxmaps)
})
sf_use_s2(FALSE)
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
dir.create("state_pages", showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
ONLY <- if (length(args)) sprintf("%02d", as.integer(args)) else NULL

hex <- readRDS("data/mxmunicipio_hexbin_v10.rds")
hex_st <- readRDS("data/mxstate_hexbin_v10.rds")

data("mxmunicipio.map")
df <- mxmunicipio.map
regions <- unique(df$region)
polys <- lapply(regions, function(r) {
  sub <- df[df$region == r, ]
  groups <- unique(sub$group)
  rings <- lapply(groups, function(g) {
    pts <- sub[sub$group == g, ]
    pts <- pts[order(pts$order), ]
    co <- cbind(pts$long, pts$lat)
    if (nrow(co) < 4) return(NULL)
    if (!all(co[1, ] == co[nrow(co), ])) co <- rbind(co, co[1, ])
    list(coords = co, hole = any(pts$hole))
  })
  rings <- Filter(Negate(is.null), rings)
  if (!length(rings)) return(st_polygon())
  st_polygon(c(lapply(rings[!sapply(rings, `[[`, "hole")], `[[`, "coords"),
               lapply(rings[ sapply(rings, `[[`, "hole")], `[[`, "coords")))
})
real <- st_transform(
  st_make_valid(st_sf(region = regions, geometry = st_sfc(polys, crs = 4326))),
  6372)
real$state_code <- substr(real$region, 1, 2)
real_st <- real |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop")

data("df_mxmunicipio_2020")
census <- df_mxmunicipio_2020 |>
  select(region, indigenous_language, pop) |>
  mutate(indigenous_share = indigenous_language / pop)
lookup <- df_mxmunicipio_2020 |>
  select(state_code, state_name, state_abbr) |>
  distinct()

hex <- hex |> left_join(census, by = "region")
real <- real |> left_join(census, by = "region")

panel <- function(cells, borders, title, lw = 0.08, blw = 0.5) {
  ggplot() +
    geom_sf(data = cells, aes(fill = indigenous_share),
            color = "grey30", linewidth = lw) +
    geom_sf(data = borders, fill = NA, color = "black", linewidth = blw) +
    scale_fill_distiller(palette = "YlGnBu", direction = 1,
                         labels = percent_format(accuracy = 1),
                         name = "Indigenous\nlanguage share",
                         na.value = "grey90") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "right") +
    labs(title = title)
}

# ---- national -------------------------------------------------------------
if (is.null(ONLY)) {
  p <- panel(hex, hex_st, "v10 — uniform hexbin (1 hex = 1 municipio)") +
    panel(real, real_st, "True municipio territory") +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Mexico — v10 uniform hexagon cartogram vs true territory",
      subtitle = paste("Every hexagon is identical (r = 17.5 km, 792 km2), so",
                       "each state's size tracks its MUNICIPIO COUNT, not its area.",
                       "\nShare of population speaking an indigenous language,",
                       "2020 census."),
      theme = theme(plot.title = element_text(face = "bold", size = 18),
                    plot.subtitle = element_text(color = "grey40", size = 11)))
  ggsave("state_pages/hexbin_v10_national.png", p,
         width = 20, height = 9, dpi = 110, bg = "white")
  cat("Wrote state_pages/hexbin_v10_national.png\n")
}

# ---- per state ------------------------------------------------------------
codes <- if (is.null(ONLY)) sort(unique(hex$state_code)) else ONLY
for (s in codes) {
  nm <- lookup$state_name[match(s, lookup$state_code)]
  ab <- lookup$state_abbr[match(s, lookup$state_code)]
  h <- hex[hex$state_code == s, ]
  r <- real[real$state_code == s, ]
  p <- (panel(h, hex_st[hex_st$state_code == s, ],
              sprintf("v10 hexbin — %d hexes", nrow(h)), lw = 0.15) +
        panel(r, real_st[real_st$state_code == s, ],
              "True territory", lw = 0.15)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("%s (%s) — uniform hexbin vs true territory", nm, ab),
      subtitle = "Cells are identical regular hexagons; only position and adjacency are optimised",
      theme = theme(plot.title = element_text(face = "bold", size = 17),
                    plot.subtitle = element_text(color = "grey40", size = 11)))
  f <- sprintf("state_pages/hexbin_v10_%s_%s.png", s, ab)
  ggsave(f, p, width = 17, height = 9, dpi = 110, bg = "white")
  cat("Wrote", f, "\n")
}
