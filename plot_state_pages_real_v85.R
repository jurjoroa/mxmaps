suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(scales); library(patchwork)
  library(mxmaps)
})
sf_use_s2(FALSE)

# True territory vs v7 hex cartogram, one page per state (plus a national
# overview): how much distortion does the final layout introduce?
# Output: state_pages/real_v85_<code>_<abbr>.png + real_v85_national.png

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
dir.create("state_pages", showWarnings = FALSE)

v7 <- readRDS("data/mxmunicipio_hex_sf_v8_5.rds")
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

data("df_mxmunicipio_2020")
census <- df_mxmunicipio_2020 |>
  select(region, indigenous_language, pop) |>
  mutate(indigenous_share = indigenous_language / pop)
state_lookup <- df_mxmunicipio_2020 |>
  select(state_code, state_name, state_abbr) |>
  distinct()

hex_v7 <- v7 |> left_join(census, by = "region")
real <- mxmun_sf |> left_join(census, by = "region")

mk_panel <- function(hex_data, real_data, xlim, ylim, title, show_real,
                     lw = 0.10) {
  p <- ggplot() +
    scale_fill_distiller(palette = "YlGnBu", direction = 1,
                         name = "Indigenous\nlanguage share",
                         labels = percent_format(accuracy = 1),
                         na.value = "grey90")
  if (show_real) {
    p <- p + geom_sf(data = real_data, aes(fill = indigenous_share),
                     color = "grey30", linewidth = 0.15)
  } else {
    p <- p + geom_sf(data = hex_data, aes(fill = indigenous_share),
                     color = "grey40", linewidth = lw) +
      geom_sf(data = borders, fill = NA, color = "black", linewidth = 0.5)
  }
  p + coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "right",
          legend.key.height = unit(1.2, "cm"),
          plot.margin = margin(4, 4, 4, 4)) +
    labs(title = title)
}

# --- National overview --------------------------------------------------------
p <- mk_panel(hex_v7, real, c(-118.5, -86.5), c(14, 32.8),
              "v8.5 — equal-area + cohesion", show_real = FALSE) |
     mk_panel(real, real, c(-118.5, -86.5), c(14, 32.8),
              "True municipio territory", show_real = TRUE) +
  plot_annotation(
    title = "Mexico — v8.5 equal-area municipio cartogram vs true territory",
    subtitle = "Share of population speaking an indigenous language, 2020 census",
    theme = theme(plot.title = element_text(face = "bold", size = 18),
                  plot.subtitle = element_text(color = "grey40", size = 12)))
ggsave("state_pages/real_v85_national.png", p,
       width = 16, height = 7, dpi = 150, bg = "white")
cat("Wrote state_pages/real_v85_national.png\n")

# --- Per state ----------------------------------------------------------------
for (s in sort(unique(mxmun_sf$state_code))) {
  b <- st_bbox(state_sf$geometry[state_sf$state_code == s])
  pad_x <- 0.18 * (b["xmax"] - b["xmin"])
  pad_y <- 0.18 * (b["ymax"] - b["ymin"])
  xlim <- c(b["xmin"] - pad_x, b["xmax"] + pad_x)
  ylim <- c(b["ymin"] - pad_y, b["ymax"] + pad_y)
  nm <- state_lookup$state_name[state_lookup$state_code == s][1]
  ab <- state_lookup$state_abbr[state_lookup$state_code == s][1]

  p7 <- mk_panel(hex_v7 |> filter(state_code == s), NULL, xlim, ylim,
                 "v8.5 — equal-area + cohesion", show_real = FALSE)
  real_s <- real |> filter(state_code == s)
  pr <- mk_panel(real_s, real_s, xlim, ylim,
                 "True territory", show_real = TRUE)
  p <- p7 + pr + plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("%s (%s) — distortion check", nm, ab),
      subtitle = "v8.5 equal-area municipio cartogram vs true municipio territory",
      theme = theme(plot.title = element_text(face = "bold", size = 18),
                    plot.subtitle = element_text(color = "grey40", size = 12)))
  fname <- sprintf("state_pages/real_v85_%s_%s.png", s, ab)
  ggsave(fname, p, width = 15, height = 8.5, dpi = 150, bg = "white")
  cat(sprintf("Wrote %s\n", fname))
}
