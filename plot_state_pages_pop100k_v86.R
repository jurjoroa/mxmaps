suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(scales); library(patchwork)
  library(mxmaps)
})
sf_use_s2(FALSE)

# Distortion check with population: inhabitants per 100k, by municipio.
#   - state pages:   share of the STATE's population, per 100k
#   - national page: share of the NATIONAL population, per 100k
# True territory vs the v8.6 hex cartogram, same fill, same legend.
# Output: state_pages/pop100k_v86_<code>_<abbr>.png + pop100k_v86_national.png

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
dir.create("state_pages", showWarnings = FALSE)

v8 <- readRDS("data/mxmunicipio_hex_sf_v8_6.rds")
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
  select(region, state_code, pop)
nat_pop <- sum(census$pop)
state_pop <- census |> group_by(state_code) |>
  summarise(state_pop = sum(pop), .groups = "drop")
census <- census |>
  left_join(state_pop, by = "state_code") |>
  mutate(
    pop100k_state = pop / state_pop * 1e5,
    pop100k_nat = pop / nat_pop * 1e5
  ) |>
  select(region, pop100k_state, pop100k_nat)
state_lookup <- df_mxmunicipio_2020 |>
  select(state_code, state_name, state_abbr) |>
  distinct()

hex_v8 <- v8 |> left_join(census, by = "region")
real <- mxmun_sf |> left_join(census, by = "region")

mk_panel <- function(data, xlim, ylim, title, fill_var, show_real,
                     lw = 0.10) {
  p <- ggplot() +
    scale_fill_distiller(palette = "YlGnBu", direction = 1,
                         name = "Inhabitants\nper 100k",
                         labels = comma_format(accuracy = 1),
                         na.value = "grey90")
  if (show_real) {
    p <- p + geom_sf(data = data, aes(fill = .data[[fill_var]]),
                     color = "grey30", linewidth = 0.15)
  } else {
    p <- p + geom_sf(data = data, aes(fill = .data[[fill_var]]),
                     color = "grey40", linewidth = 0.10) +
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

# --- National overview (national rate) ----------------------------------------
p8 <- mk_panel(hex_v8, c(-118.5, -86.5), c(14, 32.8),
               "v8.6 — hex cartogram", "pop100k_nat", FALSE)
pr <- mk_panel(real, c(-118.5, -86.5), c(14, 32.8),
               "True municipio territory", "pop100k_nat", TRUE)
p <- p8 + pr + plot_layout(guides = "collect") +
  plot_annotation(
    title = "Mexico — inhabitants per 100k (national), v8.6 vs true territory",
    subtitle = "Each municipio's population as a share of the national population, per 100,000",
    theme = theme(plot.title = element_text(face = "bold", size = 18),
                  plot.subtitle = element_text(color = "grey40", size = 12)))
ggsave("state_pages/pop100k_v86_national.png", p,
       width = 16, height = 7, dpi = 150, bg = "white")
cat("Wrote state_pages/pop100k_v86_national.png\n")

# --- Per state (state rate) ----------------------------------------------------
for (s in sort(unique(mxmun_sf$state_code))) {
  b <- st_bbox(state_sf$geometry[state_sf$state_code == s])
  pad_x <- 0.18 * (b["xmax"] - b["xmin"])
  pad_y <- 0.18 * (b["ymax"] - b["ymin"])
  xlim <- c(b["xmin"] - pad_x, b["xmax"] + pad_x)
  ylim <- c(b["ymin"] - pad_y, b["ymax"] + pad_y)
  nm <- state_lookup$state_name[state_lookup$state_code == s][1]
  ab <- state_lookup$state_abbr[state_lookup$state_code == s][1]

  p8 <- mk_panel(hex_v8 |> filter(state_code == s), xlim, ylim,
                 "v8.6 — hex cartogram", "pop100k_state", FALSE)
  real_s <- real |> filter(state_code == s)
  pr <- mk_panel(real_s, xlim, ylim,
                 "True territory", "pop100k_state", TRUE)
  p <- p8 + pr + plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("%s (%s) — inhabitants per 100k (state)", nm, ab),
      subtitle = "Each municipio's population as a share of the state population, per 100,000",
      theme = theme(plot.title = element_text(face = "bold", size = 18),
                    plot.subtitle = element_text(color = "grey40", size = 12)))
  fname <- sprintf("state_pages/pop100k_v86_%s_%s.png", s, ab)
  ggsave(fname, p, width = 15, height = 8.5, dpi = 150, bg = "white")
  cat(sprintf("Wrote %s\n", fname))
}
