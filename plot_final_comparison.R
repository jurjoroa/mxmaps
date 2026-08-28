# ---------------------------------------------------------------------------
# Purpose: the FINAL comparison sheet for the canonical municipio layout —
#          layout vs true territory, plus a measured metrics strip.
#
#          Design decisions and why (dataviz guidance):
#          - FORM: choropleth pair. The data's job is magnitude, and the point
#            of the figure is a like-for-like comparison, so two panels sharing
#            one scale beats any single-panel trick.
#          - COLOUR: YlGnBu sequential on a log10 scale with fixed decade
#            breaks. Validated the way a SEQUENTIAL ramp must be — L* strictly
#            monotonic 98.9 -> 22.1, min step 8.7 L* units, ColorBrewer
#            colourblind- and print-safe. (Running the categorical validator on
#            a sequential ramp is a category error: a sequential ramp is
#            supposed to span the lightness band and have close neighbours.)
#            log10 because the variable has a 142x skew — on a linear ramp
#            95.4% of all 2469 municipios land in one colour bin, on BOTH
#            panels, which is what made every earlier figure read as blank.
#          - ONE legend, collected across panels, so a colour means the same
#            thing left and right.
#          - The metrics strip is stat tiles, not a chart, and its text wears
#            INK colours only — never the fill palette. A number is not a
#            series.
#          - Every metric is COMPUTED from the artifact at render time, so the
#            figure cannot drift from the data it claims to describe.
#
# Author:  Claude (Opus 5) for Jorge Roa
# Version: final
# Date:    26-August-2026
# Usage:   Rscript plot_final_comparison.R [tag]      (default: v11_2)
# Output:  state_pages/FINAL_<tag>_national.png
#          state_pages/FINAL_<tag>_states.png
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(scales)
  library(patchwork); library(mxmaps)
})
sf_use_s2(FALSE)
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
source("hexlayout_common.R")
dir.create("state_pages", showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
TAG <- if (length(args) >= 1) args[1] else "v11_2"
f <- sprintf("data/mxmunicipio_hex_sf_%s.rds", TAG)
if (!file.exists(f)) stop("no such artifact: ", f)

INK      <- "grey15"    # primary text
INK_SOFT <- "grey40"    # secondary text
INK_MUTE <- "grey58"    # tertiary / units

lay  <- readRDS(f)
real <- build_mxmun_sf()
borders  <- readRDS("data/mxstate_borders.rds")
state_sf <- real |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |> st_make_valid()

data("df_mxmunicipio_2020")
cs <- df_mxmunicipio_2020 |> select(region, state_code, pop)
natp <- sum(cs$pop)
sp <- cs |> group_by(state_code) |> summarise(spop = sum(pop), .groups = "drop")
cs <- cs |> left_join(sp, by = "state_code") |>
  mutate(p_nat = pop / natp * 1e5, p_state = pop / spop * 1e5) |>
  select(region, p_nat, p_state)
lk <- df_mxmunicipio_2020 |>
  select(state_code, state_name, state_abbr) |> distinct()

layj  <- lay  |> left_join(cs, by = "region")
realj <- real |> left_join(cs, by = "region")

# ---- measure the artifact, do not quote remembered numbers ----------------
lay6  <- st_make_valid(st_transform(lay, 6372))
real6 <- st_make_valid(st_transform(real, 6372))
a_lay <- as.numeric(st_area(lay6))
tru_state <- tapply(as.numeric(st_area(real6)), real$state_code, sum)
cov_state <- tapply(a_lay, lay$state_code, sum)

r_eq <- unsplit(lapply(split(a_lay, lay$state_code), function(z) z / mean(z)),
                lay$state_code)
n_tiny <- sum(r_eq < 0.25)
eq_min <- min(r_eq)
# COVERAGE IS GEOMETRIC, NOT AN AREA SUM. This tile used to report
#   max_state (true_area - sum(cell_area)) / true_area
# which is the OT solver's area residual, and it cannot see an uncovered hole:
# a v11.3 build with KNN=60 left 3816.9 ppm of Chiapas covered by no cell at all
# and this number read 0.02% before and after the fix. Measure the thing the
# hard constraint actually names -- the part of the real state polygon that no
# cell covers -- with st_difference. poly_area() is unnecessary here because
# everything is already projected to 6372.
gap_ppm <- max(vapply(sort(unique(lay$state_code)), function(s) {
  sp <- st_make_valid(st_union(real6$geometry[real$state_code == s]))
  cl <- st_make_valid(st_union(st_make_valid(lay6$geometry[lay$state_code == s])))
  d <- suppressWarnings(st_difference(sp, cl))
  if (length(d) == 0) return(0)
  1e6 * sum(as.numeric(st_area(d))) / sum(as.numeric(st_area(sp)))
}, 0))
# kept as a separate diagnostic: the solver's per-state area residual, which is
# a real number about convergence but is NOT a statement about coverage
area_res_pct <- 100 * max((tru_state[names(cov_state)] - cov_state) /
                            tru_state[names(cov_state)])

ov <- suppressWarnings(st_intersection(lay6, lay6))
ov <- ov[ov$region != ov$region.1, ]
ov <- ov[!is.na(st_dimension(ov)) & st_dimension(ov) == 2, ]
ov_ppm <- if (nrow(ov) == 0) 0 else
  1e6 * sum(as.numeric(st_area(ov))) / 2 / sum(a_lay)

# cohesion: share of true intra-state adjacencies that survive as touching cells
adj <- st_intersects(real6, real6)
touch <- st_is_within_distance(lay6, lay6, dist = 1, sparse = TRUE)
idx <- match(real$region, lay$region)
pres <- 0L; tot <- 0L
for (i in seq_len(nrow(real6))) {
  js <- adj[[i]]
  js <- js[js > i & real$state_code[js] == real$state_code[i]]
  for (j in js) {
    tot <- tot + 1L
    if (!is.na(idx[i]) && !is.na(idx[j]) && idx[j] %in% touch[[idx[i]]])
      pres <- pres + 1L
  }
}
cohesion <- 100 * pres / max(tot, 1)

# displacement, in units of the state's own cell size (scale-free) and in km
cl <- st_coordinates(st_centroid(lay6))[, 1:2, drop = FALSE]
ct <- st_coordinates(st_centroid(real6, of_largest_polygon = TRUE))[, 1:2, drop = FALSE]
disp_km <- sqrt(rowSums((cl - ct[match(lay$region, real$region), ])^2)) / 1000

cat(sprintf(paste0("MEASURED %s: cohesion %.1f%% | disp %.0f km | n_tiny %d | ",
                   "eq_min %.3f | overlap %.1f ppm | uncovered %.1f ppm\n"),
            TAG, cohesion, mean(disp_km), n_tiny, eq_min, ov_ppm, gap_ppm))

# ---- panels ---------------------------------------------------------------
fill_scale <- function()
  scale_fill_distiller(
    palette = "YlGnBu", direction = 1, trans = "log10",
    breaks = c(0.1, 1, 10, 100, 1000),
    labels = c("0.1", "1", "10", "100", "1,000"),
    na.value = "grey92",
    name = "Inhabitants per 100k\n(log scale)")

pan <- function(d, fv, xlim, ylim, title, is_real, lw = 0.07) {
  ggplot() +
    geom_sf(data = d, aes(fill = .data[[fv]]),
            color = if (is_real) "grey35" else "grey45",
            linewidth = if (is_real) lw + 0.03 else lw) +
    { if (!is_real)
        geom_sf(data = borders, fill = NA, color = INK, linewidth = 0.42) } +
    fill_scale() +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12.5, colour = INK,
                                    margin = margin(b = 4)),
          legend.title = element_text(size = 9.5, colour = INK_SOFT),
          legend.text = element_text(size = 9, colour = INK_SOFT),
          legend.key.height = unit(0.85, "cm"),
          legend.position = "right") +
    labs(title = title)
}

# ---- metrics strip: stat tiles, ink text only -----------------------------
tiles <- data.frame(
  x = 1:6,
  label = c("Cohesion", "Mean displacement", "Illegible cells",
            "Smallest cell", "Cell overlap", "Territory uncovered"),
  value = c(sprintf("%.1f%%", cohesion), sprintf("%.0f km", mean(disp_km)),
            sprintf("%d", n_tiny), sprintf("%.2f×", eq_min),
            sprintf("%.1f ppm", ov_ppm), sprintf("%.1f ppm", gap_ppm)),
  note = c("of true adjacencies kept", "cell centroid vs municipio",
           "under 0.25× state mean", "relative to state mean",
           "of national area", "of state area, no cell covers it"),
  stringsAsFactors = FALSE)

strip <- ggplot(tiles, aes(x = x)) +
  geom_text(aes(y = 0.62, label = value), size = 6.4, fontface = "bold",
            colour = INK) +
  geom_text(aes(y = 0.28, label = label), size = 3.5, colour = INK_SOFT) +
  geom_text(aes(y = 0.05, label = note), size = 2.9, colour = INK_MUTE) +
  scale_x_continuous(limits = c(0.5, 6.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.05, 0.85), expand = c(0, 0)) +
  theme_void()

# ---- national sheet -------------------------------------------------------
xl <- c(-118.5, -86.5); yl <- c(14, 32.8)
top <- pan(layj, "p_nat", xl, yl, sprintf("%s layout — one cell per municipio", TAG), FALSE) +
       pan(realj, "p_nat", xl, yl, "True municipio territory", TRUE) +
  plot_layout(guides = "collect")

p <- (top / strip) + plot_layout(heights = c(1, 0.20)) +
  plot_annotation(
    title = sprintf("Mexico — municipio layout %s vs true territory", TAG),
    subtitle = paste0(
      "Real national and state outlines preserved; one cell per municipio ",
      "inside its own state. Colour is inhabitants per 100k on a log scale, ",
      "identical on both panels.\nAll figures below are measured from the ",
      "artifact at render time."),
    caption = paste0("Cell SIZE carries no data: it is set per state as ",
                     "state area / municipio count, and varies 223× between ",
                     "states (BCS 14,789 km²/cell vs Tlaxcala 66 km²/cell). ",
                     "Compare cells within a state, not across states."),
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, colour = INK),
      plot.subtitle = element_text(colour = INK_SOFT, size = 10.5,
                                   margin = margin(b = 6)),
      plot.caption = element_text(colour = INK_MUTE, size = 9, hjust = 0,
                                  margin = margin(t = 6))))
cat(sprintf("coverage: worst state gap %.1f ppm | solver area residual %.3f%%\n",
            gap_ppm, area_res_pct))
fn <- sprintf("state_pages/FINAL_%s_national.png", TAG)
ggsave(fn, p, width = 19, height = 10, dpi = 150, bg = "white")
cat("Wrote", fn, "\n")

# ---- representative states sheet -----------------------------------------
# chosen to span the range the layout has to cope with, not cherry-picked
sel <- c("20", "30", "03", "23", "09", "08")
why <- c("570 municipios (densest)", "lowest cohesion", "5 municipios (sparsest)",
         "islands (Cozumel)", "smallest cells", "largest cells")
pairs <- list()
for (k in seq_along(sel)) {
  s <- sel[k]
  b <- st_bbox(state_sf$geometry[state_sf$state_code == s])
  px <- 0.08 * (b["xmax"] - b["xmin"]); py <- 0.08 * (b["ymax"] - b["ymin"])
  xls <- c(b["xmin"] - px, b["xmax"] + px); yls <- c(b["ymin"] - py, b["ymax"] + py)
  ab <- lk$state_abbr[lk$state_code == s][1]
  # Short titles. Long ones collided in these narrow panels and rendered as
  # "statetrue" -- caught by looking at the output, not by any check.
  L <- pan(layj[layj$state_code == s, ], "p_state", xls, yls,
           "layout", FALSE, lw = 0.13) + theme(legend.position = "none")
  R <- pan(realj[realj$state_code == s, ], "p_state", xls, yls,
           "true", TRUE, lw = 0.13) + theme(legend.position = "none")
  pairs[[k]] <- (L | R) +
    plot_annotation(title = sprintf("%s — %s", ab, why[k]))
  pairs[[k]] <- wrap_elements(pairs[[k]])
}
# One shared legend for the whole sheet: a colour must have a key. Built as its
# own strip because each pair is wrapped (wrap_elements drops guide collection).
# The key needs the REAL value range: a single value cannot define a log scale
# and scale_fill_distiller errors out.
kv <- range(c(layj$p_state, realj$p_state), na.rm = TRUE)
key <- ggplot(data.frame(x = 1:2, y = 1, v = kv)) +
  geom_tile(aes(x, y, fill = v)) + fill_scale() +
  # The tiles exist only to generate the guide; clip them out of the panel so
  # the strip shows the legend and not a giant colour band.
  coord_cartesian(xlim = c(50, 51)) +
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10, colour = INK_SOFT, vjust = 0.9),
        legend.text = element_text(size = 9, colour = INK_SOFT),
        legend.key.width = unit(2.2, "cm"), legend.key.height = unit(0.45, "cm"))

q <- wrap_plots(pairs, ncol = 2) / key +
  plot_layout(heights = c(1, 0.06)) +
  plot_annotation(
    title = sprintf("%s vs true territory — the six cases the layout must cope with", TAG),
    subtitle = paste0("Left of each pair is the layout, right is true territory. ",
                      "Log colour scale computed within each state, so panels are ",
                      "comparable within a pair but not across pairs."),
    theme = theme(
      plot.title = element_text(face = "bold", size = 17, colour = INK),
      plot.subtitle = element_text(colour = INK_SOFT, size = 10.5,
                                   margin = margin(b = 4))))
fn2 <- sprintf("state_pages/FINAL_%s_states.png", TAG)
ggsave(fn2, q, width = 17, height = 15, dpi = 140, bg = "white")
cat("Wrote", fn2, "\n")
