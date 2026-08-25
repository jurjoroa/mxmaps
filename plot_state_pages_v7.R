suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(scales); library(patchwork)
  library(mxmaps)
})
sf_use_s2(FALSE)

# One comparison page per state: v5 | v6, high resolution, shared legend.
# Output: state_pages/v6_<code>_<abbr>.png (32 pages)

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
dir.create("state_pages", showWarnings = FALSE)

v5  <- readRDS("data/mxmunicipio_hex_sf_v5.rds")
v6  <- readRDS("data/mxmunicipio_hex_sf_v6.rds")
v7  <- readRDS("data/mxmunicipio_hex_sf_v7.rds")
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

hex_v5 <- v5 |> left_join(census, by = "region")
hex_v6 <- v6 |> left_join(census, by = "region")
hex_v7 <- v7 |> left_join(census, by = "region")

n_by_state <- table(mxmun_sf$state_code)

mk_panel <- function(hex_data, xlim, ylim, title) {
  ggplot() +
    geom_sf(data = hex_data, aes(fill = indigenous_share),
            color = "grey40", linewidth = 0.15) +
    scale_fill_distiller(palette = "YlGnBu", direction = 1,
                         name = "Indigenous\nlanguage share",
                         labels = percent_format(accuracy = 1),
                         na.value = "grey90") +
    geom_sf(data = borders, fill = NA, color = "black", linewidth = 0.5) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "right",
          legend.key.height = unit(1.2, "cm"),
          plot.margin = margin(4, 4, 4, 4)) +
    labs(title = title)
}

for (s in names(n_by_state)) {
  b <- st_bbox(state_sf$geometry[state_sf$state_code == s])
  pad_x <- 0.18 * (b["xmax"] - b["xmin"])
  pad_y <- 0.18 * (b["ymax"] - b["ymin"])
  xlim <- c(b["xmin"] - pad_x, b["xmax"] + pad_x)
  ylim <- c(b["ymin"] - pad_y, b["ymax"] + pad_y)
  nm <- state_lookup$state_name[state_lookup$state_code == s][1]
  ab <- state_lookup$state_abbr[state_lookup$state_code == s][1]

  p5 <- mk_panel(hex_v6 |> filter(state_code == s), xlim, ylim,
                 "v6 — equal-area cells (CVT)")
  p6 <- mk_panel(hex_v7 |> filter(state_code == s), xlim, ylim,
                 "v7 — equal-area + cohesion")
  p <- p5 + p6 + plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("%s (%s) — %d municipios", nm, ab, n_by_state[s]),
      subtitle = "Municipio hex cartogram QA — black line = true state border",
      theme = theme(plot.title = element_text(face = "bold", size = 18),
                    plot.subtitle = element_text(color = "grey40", size = 12)))
  fname <- sprintf("state_pages/v7_%s_%s.png", s, ab)
  ggsave(fname, p, width = 15, height = 8.5, dpi = 150, bg = "white")
  cat(sprintf("Wrote %s\n", fname))
}
