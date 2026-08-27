# ------------------------------------------------------------------------------
# Purpose: Honest presentation of the municipio hex layout at national scale.
#          Real state outlines are fixed, so a state's TOTAL rendered area is
#          identical to its true territory: one BCS cell carries 222x the ink
#          of one Tlaxcala cell while standing for the same single data point.
#          Three figures answer that: (1) a national page with an explicit
#          km2-per-cell size reference, (2) the municipio layout beside the
#          state-level choropleth that is actually unbiased nationally, and
#          (3) a 32-panel small-multiples sheet, each state in its own frame,
#          which is the presentation where the layout is genuinely fair.
# Author : Jorge Roa (with Claude Code)
# Version: 1.0
# Date   : 2026-08-26
# Usage  : Rscript plot_national_honest.R
# Output : state_pages/honest_national_v88.png
#          state_pages/honest_municipio_vs_state_v88.png
#          state_pages/honest_small_multiples_v88.png
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(scales)
  library(patchwork); library(mxmaps)
})
sf_use_s2(FALSE)

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
dir.create("state_pages", showWarnings = FALSE)

hex <- readRDS("data/mxmunicipio_hex_sf_v8_8.rds")
borders <- readRDS("data/mxstate_borders.rds")

# --- true municipio polygons ---------------------------------------------------
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
  if (length(rings) == 0) return(st_polygon())
  st_polygon(c(lapply(rings[!sapply(rings, "[[", "hole")], "[[", "coords"),
               lapply(rings[ sapply(rings, "[[", "hole")], "[[", "coords")))
})
mx <- st_make_valid(st_sf(region = regions,
                          geometry = st_sfc(polys, crs = 4326)))
mx$state_code <- substr(mx$region, 1, 2)
state_sf <- mx |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_make_valid()

# --- variable: inhabitants per 100k -------------------------------------------
data("df_mxmunicipio_2020")
cen <- df_mxmunicipio_2020
nat_pop <- sum(cen$pop)
state_pop <- cen |>
  group_by(state_code) |>
  summarise(state_pop = sum(pop), .groups = "drop")
cen <- cen |>
  left_join(state_pop, by = "state_code") |>
  mutate(pop100k_nat = pop / nat_pop * 1e5,
         pop100k_state = pop / state_pop * 1e5)
lookup <- cen |> select(state_code, state_name, state_abbr) |> distinct()

hex_d <- hex |>
  left_join(cen |> select(region, pop100k_nat, pop100k_state), by = "region")

# --- cell areas, km2 -----------------------------------------------------------
hex_km2 <- as.numeric(st_area(st_transform(hex, 6372))) / 1e6
cell <- data.frame(region = hex$region, state_code = hex$state_code,
                   km2 = hex_km2)
per_state <- cell |>
  group_by(state_code) |>
  summarise(n = n(), km2_per_cell = mean(km2), .groups = "drop") |>
  left_join(lookup, by = "state_code")
nat_avg <- sum(cell$km2) / nrow(cell)

XL <- c(-118.5, -86.5)
YL <- c(14, 32.8)

# Population is extremely skewed (0.02 to 1,500 per 100k across municipios),
# so a linear ramp collapses 2,400 cells into the palest bin. Log10 keeps the
# palette and direction and makes the pattern readable.
fill_scale <- function(nm, acc = 1) {
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = nm,
                       trans = "log10",
                       labels = comma_format(accuracy = acc),
                       na.value = "grey90")
}

mk_hex_panel <- function(title, sub) {
  ggplot() +
    geom_sf(data = hex_d, aes(fill = pop100k_nat), color = "grey45",
            linewidth = 0.07) +
    geom_sf(data = borders, fill = NA, color = "black", linewidth = 0.45) +
    fill_scale("Inhabitants\nper 100k\n(national)\nlog scale", 0.1) +
    coord_sf(xlim = XL, ylim = YL, expand = FALSE) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(color = "grey35", size = 10),
          legend.position = "right",
          legend.key.height = unit(1.1, "cm")) +
    labs(title = title, subtitle = sub)
}

# ==============================================================================
# FIGURE 1 — national page with an explicit km2-per-cell size reference
# ==============================================================================
# Reference cells: the median-area cell of a sparse, a middling and a dense
# state, moved to a common origin and drawn at a common metric scale.
ref_codes <- c("03", "14", "29")           # BCS, Jalisco, Tlaxcala
hex_m <- st_transform(hex, 6372)
ref_list <- lapply(ref_codes, function(sc) {
  idx <- which(hex$state_code == sc)
  pick <- idx[which.min(abs(hex_km2[idx] - median(hex_km2[idx])))]
  g <- st_geometry(hex_m)[pick]
  ctr <- st_coordinates(st_centroid(g))[1, 1:2]
  st_sfc(st_geometry(g)[[1]] - ctr, crs = 6372)
})
ref_sf <- st_sf(
  state_code = ref_codes,
  abbr = per_state$state_abbr[match(ref_codes, per_state$state_code)],
  n = per_state$n[match(ref_codes, per_state$state_code)],
  km2 = per_state$km2_per_cell[match(ref_codes, per_state$state_code)],
  geometry = do.call(c, ref_list)
)
# lay them out left to right, spaced by the widest cell
w <- max(sapply(ref_list, function(g) diff(st_bbox(g)[c(1, 3)])))
offs <- (seq_along(ref_codes) - 1) * w * 3.5
ref_sf$geometry <- do.call(c, lapply(seq_along(ref_list), function(i)
  st_sfc(st_geometry(ref_list[[i]])[[1]] + c(offs[i], 0), crs = 6372)))
ref_sf$lab <- sprintf("%s\n1 cell = %s km2\n(%s x national mean)",
                      ref_sf$abbr, comma(round(ref_sf$km2)),
                      formatC(ref_sf$km2 / nat_avg, format = "f", digits = 2))
lab_xy <- st_coordinates(st_centroid(ref_sf))
ref_sf$lx <- lab_xy[, 1]
ref_sf$ly <- st_bbox(ref_sf)[2] - 0.20 * w

p_ref <- ggplot(ref_sf) +
  geom_sf(fill = "#41b6c4", color = "grey20", linewidth = 0.4) +
  geom_text(aes(x = lx, y = ly, label = lab), size = 3.4,
            lineheight = 1.1, vjust = 1) +
  coord_sf(datum = NA, clip = "off",
           ylim = c(st_bbox(ref_sf)[2] - 1.15 * w,
                    st_bbox(ref_sf)[4] + 0.20 * w)) +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "grey35", size = 9.5),
        plot.margin = margin(6, 20, 6, 20)) +
  labs(title = "SIZE REFERENCE — one cell = one municipio, drawn to scale",
       subtitle = paste("Every cell is one data point, but a Baja California",
                        "Sur cell covers 222x the area of a Tlaxcala cell.",
                        "\nBig cells are big states with few municipios, not",
                        "big values. Read colour, never size."))

p1 <- mk_hex_panel(
  "Mexico - inhabitants per 100k (national), v8.8 municipio layout",
  paste("2,469 cells, one per municipio, inside the real state outlines.",
        "Cell AREA carries no information.")) 

fig1 <- p1 / p_ref + plot_layout(heights = c(3, 1)) +
  plot_annotation(
    title = "The national municipio page, told honestly",
    subtitle = paste0("State outlines are exact and each state holds exactly ",
                      "its own municipio count, so cell size is forced to ",
                      "vary 222x between states.\nAt state granularity this ",
                      "page carries the same area bias as a true territory ",
                      "map (47.5% of ink misallocated, both maps)."),
    theme = theme(plot.title = element_text(face = "bold", size = 19),
                  plot.subtitle = element_text(color = "grey35", size = 11)))
ggsave("state_pages/honest_national_v88.png", fig1,
       width = 15, height = 12, dpi = 150, bg = "white")
cat("Wrote state_pages/honest_national_v88.png\n")

# ==============================================================================
# FIGURE 2 — municipio layout beside the state-level choropleth
# ==============================================================================
st_var <- state_pop |>
  mutate(pop100k_nat = state_pop / nat_pop * 1e5) |>
  left_join(lookup, by = "state_code")
state_d <- state_sf |> left_join(st_var, by = "state_code")

p_state <- ggplot() +
  geom_sf(data = state_d, aes(fill = pop100k_nat), color = "grey20",
          linewidth = 0.3) +
  fill_scale("Inhabitants\nper 100k\n(national)\nlog scale", 1) +
  coord_sf(xlim = XL, ylim = YL, expand = FALSE) +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "grey35", size = 10),
        legend.position = "right",
        legend.key.height = unit(1.1, "cm")) +
  labs(title = "State level - 32 areas, 32 data points",
       subtitle = paste("Unbiased at national scale: one polygon per data",
                        "point, no cell-size confound."))

p_mun <- mk_hex_panel(
  "Municipio level - 2,469 cells, 2,469 data points",
  paste("Fair WITHIN each state (within-state ink variance -80.5%),",
        "biased BETWEEN states (222x)."))

fig2 <- p_mun + p_state +
  plot_annotation(
    title = "Which national map is actually unbiased?",
    subtitle = paste0("Left: the municipio layout. Right: the same variable ",
                      "aggregated to states. Only the right-hand map gives ",
                      "every data point comparable ink at national scale;\nthe ",
                      "left-hand map should be read for within-state pattern, ",
                      "or on the per-state pages."),
    theme = theme(plot.title = element_text(face = "bold", size = 19),
                  plot.subtitle = element_text(color = "grey35", size = 11)))
ggsave("state_pages/honest_municipio_vs_state_v88.png", fig2,
       width = 18, height = 7.5, dpi = 150, bg = "white")
cat("Wrote state_pages/honest_municipio_vs_state_v88.png\n")

# ==============================================================================
# FIGURE 3 — small multiples, every state normalised to its own frame
# ==============================================================================
codes <- sort(unique(hex$state_code))
panels <- lapply(codes, function(s) {
  d <- hex_d |> filter(state_code == s)
  b <- st_bbox(d)
  # square the frame so all 32 panels share one aspect ratio
  cx <- mean(b[c(1, 3)]); cy <- mean(b[c(2, 4)])
  h <- max(diff(b[c(1, 3)]), diff(b[c(2, 4)])) / 2 * 1.06
  ab <- lookup$state_abbr[lookup$state_code == s][1]
  rng <- range(d$pop100k_state, na.rm = TRUE)
  ggplot() +
    geom_sf(data = d, aes(fill = pop100k_state), color = "grey50",
            linewidth = 0.06) +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, limits = rng,
                         trans = "log10", guide = "none",
                         na.value = "grey90") +
    coord_sf(xlim = c(cx - h, cx + h), ylim = c(cy - h, cy + h),
             expand = FALSE) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
          plot.subtitle = element_text(size = 6.6, hjust = 0.5,
                                       color = "grey40"),
          plot.margin = margin(2, 2, 2, 2)) +
    labs(title = sprintf("%s (n=%d)", ab, nrow(d)),
         subtitle = sprintf("%s-%s per 100k | %s km2/cell",
                            comma(round(rng[1])), comma(round(rng[2])),
                            comma(round(mean(
                              cell$km2[cell$state_code == s])))))
})
fig3 <- wrap_plots(panels, ncol = 6) +
  plot_annotation(
    title = "All 32 states, each in its own frame - the fair presentation",
    subtitle = paste0("Inhabitants per 100k of the STATE's population. Each ",
                      "panel has its own colour range and its own scale, ",
                      "printed under its title.\nWithin a state, cell areas ",
                      "vary only 3.0x (median) against 51.5x on the true ",
                      "territory map, so within a panel ink tracks the data ",
                      "and not the terrain.\nDo not compare colours or cell ",
                      "sizes across panels."),
    theme = theme(plot.title = element_text(face = "bold", size = 20),
                  plot.subtitle = element_text(color = "grey35", size = 11)))
ggsave("state_pages/honest_small_multiples_v88.png", fig3,
       width = 19, height = 22, dpi = 130, bg = "white", limitsize = FALSE)
cat("Wrote state_pages/honest_small_multiples_v88.png\n")
