suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(scales); library(patchwork)
  library(mxmaps)
})
sf_use_s2(FALSE)

# Per-state visual QA: v3.1 vs v5 zoom for all 32 states, sorted by municipio
# count (dense states first — where balancing matters most). Four sheets of
# 8 states each: rows = states, columns = v3.1 | v5.

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

data("df_mxmunicipio_2020")
census <- df_mxmunicipio_2020 |>
  select(region, indigenous_language, pop) |>
  mutate(indigenous_share = indigenous_language / pop)
state_lookup <- df_mxmunicipio_2020 |>
  select(state_code, state_name, state_abbr) |>
  distinct()

hex_v31 <- v31 |> left_join(census, by = "region")
hex_v5  <- v5  |> left_join(census, by = "region")

n_by_state <- table(mxmun_sf$state_code)
state_order <- names(sort(n_by_state, decreasing = TRUE))

mk_zoom <- function(hex_data, borders_df, state_code, xlim, ylim, title) {
  ggplot() +
    geom_sf(data = hex_data, aes(fill = indigenous_share),
            color = "grey45", linewidth = 0.1) +
    scale_fill_distiller(palette = "YlGnBu", direction = 1,
                         labels = percent_format(accuracy = 1),
                         na.value = "grey90") +
    geom_sf(data = borders_df, fill = NA, color = "black", linewidth = 0.4) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          legend.position = "none",
          plot.margin = margin(2, 2, 2, 2)) +
    labs(title = title)
}

sheet_size <- 8
n_sheets <- ceiling(length(state_order) / sheet_size)
for (sh in seq_len(n_sheets)) {
  idx <- ((sh - 1) * sheet_size + 1):min(sh * sheet_size, length(state_order))
  states <- state_order[idx]
  row_list <- list()
  for (s in states) {
    b <- st_bbox(state_sf$geometry[state_sf$state_code == s])
    pad_x <- 0.18 * (b["xmax"] - b["xmin"])
    pad_y <- 0.18 * (b["ymax"] - b["ymin"])
    xlim <- c(b["xmin"] - pad_x, b["xmax"] + pad_x)
    ylim <- c(b["ymin"] - pad_y, b["ymax"] + pad_y)
    nm <- state_lookup$state_name[state_lookup$state_code == s][1]
    ab <- state_lookup$state_abbr[state_lookup$state_code == s][1]
    p31 <- mk_zoom(hex_v31 |> filter(state_code == s), borders, s,
                   xlim, ylim, sprintf("v3.1 — %s · %d", ab, n_by_state[s]))
    p5  <- mk_zoom(hex_v5  |> filter(state_code == s), borders, s,
                   xlim, ylim, sprintf("v5 — %s · %d", ab, n_by_state[s]))
    row_list[[length(row_list) + 1]] <- p31 | p5
  }
  p_sheet <- patchwork::wrap_plots(row_list, ncol = 1) +
    plot_annotation(
      title = sprintf("State-by-state QA — v3.1 vs v5 (sheet %d/%d, sorted by municipio count)",
                      sh, n_sheets),
      theme = theme(plot.title = element_text(face = "bold", size = 15)))
  fname <- sprintf("v3_1_vs_v5_states_sheet%d.png", sh)
  ggsave(fname, p_sheet, width = 13, height = 27, dpi = 110, bg = "white")
  cat(sprintf("Wrote %s (%s)\n", fname,
              paste(state_lookup$state_abbr[match(states, state_lookup$state_code)],
                    collapse = ", ")))
}
