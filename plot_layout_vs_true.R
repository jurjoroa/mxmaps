# ---------------------------------------------------------------------------
# Purpose: render ANY layout artifact against true territory, with a colour
#          scale that actually works. Replaces the plot_state_pages_pop100k_v8*
#          copies (one script, a version argument) and uses hexlayout_common.R
#          instead of re-transcribing the polygon builder a fourth time.
#
#          THE COLOUR FIX. Second iteration: log10, not quantile classes.
#          The previous pages used scale_fill_distiller() with no transform on
#          population share per 100k. That distribution has median 10.75 and
#          max 1525.6 -- a 142x skew -- so on a linear ramp **95.4% of all 2469
#          municipios fall in the lowest of nine equal-width colour bins**. The
#          pages were ~95% one pale colour, on BOTH panels, and that reading
#          failure had nothing whatsoever to do with the layout. Quantile
#          The first fix here used equal-COUNT quantile classes. That is wrong
#          for this page: these are distortion-check figures, so magnitude has
#          to survive, and equal-count bins force both panels to show the same
#          number of cells per colour -- they are made to look alike BY
#          CONSTRUCTION, which is precisely what the comparison must not do.
#          log10 with fixed decade breaks fixes the one-bin collapse, keeps
#          magnitude, and leaves the two panels honestly comparable.
#
# Author:  Claude (Opus 5) for Jorge Roa
# Version: shared
# Date:    26-August-2026
# Usage:   Rscript plot_layout_vs_true.R <tag> [state_code ...]
#          e.g. Rscript plot_layout_vs_true.R v11
#          tag names the artifact data/mxmunicipio_hex_sf_<tag>.rds
# Output:  state_pages/cls_<tag>_national.png
#          state_pages/cls_<tag>_<code>_<abbr>.png
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
if (length(args) < 1) stop("usage: Rscript plot_layout_vs_true.R <tag> [state ...]")
TAG <- args[1]
ONLY <- if (length(args) > 1) sprintf("%02d", as.integer(args[-1])) else NULL
N_CLASS <- 7

f <- sprintf("data/mxmunicipio_hex_sf_%s.rds", TAG)
if (!file.exists(f)) stop("no such artifact: ", f)
lay <- readRDS(f)

real <- build_mxmun_sf()
state_sf <- real |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_make_valid()
borders <- readRDS("data/mxstate_borders.rds")

data("df_mxmunicipio_2020")
cs <- df_mxmunicipio_2020 |> select(region, state_code, pop)
natp <- sum(cs$pop)
sp <- cs |> group_by(state_code) |> summarise(spop = sum(pop), .groups = "drop")
cs <- cs |>
  left_join(sp, by = "state_code") |>
  mutate(p_state = pop / spop * 1e5, p_nat = pop / natp * 1e5) |>
  select(region, p_state, p_nat)
lk <- df_mxmunicipio_2020 |>
  select(state_code, state_name, state_abbr) |> distinct()

lay  <- lay  |> left_join(cs, by = "region")
realj <- real |> left_join(cs, by = "region")

# Quantile classes over the pooled values of BOTH panels, so a colour means the
# same thing on the left and the right. Labels show the real value ranges,
# which a continuous ramp on a 142x-skewed variable cannot convey.
# Fixed decade breaks on a log10 scale, shared by both panels. Decades are
# data-independent, so a colour means the same absolute rate on every figure in
# the set -- unlike quantile classes, which are redefined per panel/state.
log_scale <- function() {
  scale_fill_distiller(
    palette = "YlGnBu", direction = 1, trans = "log10",
    breaks = c(0.1, 1, 10, 100, 1000),
    labels = c("0.1", "1", "10", "100", "1,000"),
    na.value = "grey90",
    name = "Inhabitants\nper 100k\n(log scale)")
}

pan <- function(d, fv, xlim, ylim, title, is_real, lw = 0.08) {
  ggplot() +
    { if (is_real)
        geom_sf(data = d, aes(fill = .data[[fv]]), color = "grey30",
                linewidth = lw + 0.04)
      else
        geom_sf(data = d, aes(fill = .data[[fv]]), color = "grey40",
                linewidth = lw) } +
    { if (!is_real)
        geom_sf(data = borders, fill = NA, color = "black", linewidth = 0.45) } +
    log_scale() +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "right",
          legend.key.height = unit(0.7, "cm")) +
    labs(title = title)
}

# ---- national -------------------------------------------------------------
if (is.null(ONLY)) {
  xl <- c(-118.5, -86.5); yl <- c(14, 32.8)
  p <- pan(lay, "p_nat", xl, yl, sprintf("%s layout", TAG), FALSE) +
       pan(realj, "p_nat", xl, yl, "True municipio territory", TRUE) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("Mexico — inhabitants per 100k, %s vs true territory", TAG),
      subtitle = paste0("log10 scale with fixed decade breaks. A linear ramp put 95.4% of ",
                        "municipios in one colour bin, which is why the earlier pages ",
                        "read as blank. Magnitude is preserved, so the two panels ",
                        "remain honestly comparable."),
      theme = theme(plot.title = element_text(face = "bold", size = 17),
                    plot.subtitle = element_text(color = "grey40", size = 10)))
  fn <- sprintf("state_pages/cls_%s_national.png", TAG)
  ggsave(fn, p, width = 19, height = 8, dpi = 130, bg = "white")
  cat("Wrote", fn, "\n")
}

# ---- per state ------------------------------------------------------------
codes <- if (is.null(ONLY)) sort(unique(real$state_code)) else ONLY
for (s in codes) {
  ls_ <- lay[lay$state_code == s, ]; rs_ <- realj[realj$state_code == s, ]
  b <- st_bbox(state_sf$geometry[state_sf$state_code == s])
  px <- 0.15 * (b["xmax"] - b["xmin"]); py <- 0.15 * (b["ymax"] - b["ymin"])
  xl <- c(b["xmin"] - px, b["xmax"] + px); yl <- c(b["ymin"] - py, b["ymax"] + py)
  nm <- lk$state_name[lk$state_code == s][1]
  ab <- lk$state_abbr[lk$state_code == s][1]
  p <- pan(ls_, "p_state", xl, yl, sprintf("%s layout", TAG), FALSE, lw = 0.14) +
       pan(rs_, "p_state", xl, yl, "True territory", TRUE, lw = 0.14) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("%s (%s) — inhabitants per 100k, %s vs true", nm, ab, TAG),
      subtitle = "log10 scale, fixed decade breaks shared by both panels",
      theme = theme(plot.title = element_text(face = "bold", size = 16),
                    plot.subtitle = element_text(color = "grey40", size = 10)))
  fn <- sprintf("state_pages/cls_%s_%s_%s.png", TAG, s, ab)
  ggsave(fn, p, width = 17, height = 8.5, dpi = 130, bg = "white")
  cat("Wrote", fn, "\n")
}
