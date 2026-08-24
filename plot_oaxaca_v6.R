suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(scales); library(patchwork)
  library(mxmaps)
})
sf_use_s2(FALSE)

# The Oaxaca example, reproduced: v3.1 | v5 | v6 (equal-area CVT cells)

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

v31 <- readRDS("data/mxmunicipio_hex_sf_v3_1.rds")
v5  <- readRDS("data/mxmunicipio_hex_sf_v5.rds")
v6  <- readRDS("/tmp/v6_state_20.rds")
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

hex_v31 <- v31 |> left_join(census, by = "region") |> filter(state_code == "20")
hex_v5  <- v5  |> left_join(census, by = "region") |> filter(state_code == "20")
hex_v6  <- v6  |> left_join(census, by = "region")

# evenness metric for the three versions (Oaxaca only)
ratio_of <- function(hex, label) {
  hex_m <- st_transform(hex, 6372)
  areas <- as.numeric(st_area(hex_m))
  share <- sum(areas) / length(areas)
  r <- areas / share
  cat(sprintf("%s Oaxaca area ratio: p5=%.2f p25=%.2f p75=%.2f p95=%.2f\n",
              label, quantile(r, .05), quantile(r, .25),
              quantile(r, .75), quantile(r, .95)))
  invisible(r)
}
ratio_of(v31 |> filter(state_code == "20"), "v3.1")
ratio_of(v5  |> filter(state_code == "20"), "v5  ")
ratio_of(v6,                                 "v6  ")

b <- st_bbox(state_sf$geometry[state_sf$state_code == "20"])
pad_x <- 0.18 * (b["xmax"] - b["xmin"]); pad_y <- 0.18 * (b["ymax"] - b["ymin"])
xlim <- c(b["xmin"] - pad_x, b["xmax"] + pad_x)
ylim <- c(b["ymin"] - pad_y, b["ymax"] + pad_y)

mk_panel <- function(hex_data, title) {
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

p <- mk_panel(hex_v31, "v3.1 — nearest absorption") +
     mk_panel(hex_v5,  "v5 — balanced absorption") +
     mk_panel(hex_v6,  "v6 — equal-area cells (CVT)") +
     plot_layout(guides = "collect") +
  plot_annotation(
    title = "Oaxaca (OAX) — 570 municipios",
    subtitle = "Municipio hex cartogram QA — black line = true state border",
    theme = theme(plot.title = element_text(face = "bold", size = 18),
                  plot.subtitle = element_text(color = "grey40", size = 12)))
ggsave("v3_1_v5_v6_oaxaca.png", p, width = 20, height = 8, dpi = 150, bg = "white")
cat("Wrote v3_1_v5_v6_oaxaca.png\n")
