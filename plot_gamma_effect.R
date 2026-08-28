# ---------------------------------------------------------------------------
# Purpose: show what re-keying the per-state gamma rule actually bought, which
#          a pair of national maps cannot show and a single number hides.
#
#          WHAT CHANGED AND WHY. The shipped rule set a state's area-variation
#          exponent from its CELL SIZE:
#            gamma = clamp(0.18 * (cs_km / cs_national - 0.75), 0, 0.45)
#          so a state got to vary its cell areas only if its cells were already
#          big. That was derived for legibility and never checked against
#          cohesion -- and it pins Veracruz (337 km2/cell) and Puebla (157) at
#          gamma = 0, which are two of the three worst states in the country.
#          The replacement picks each state's gamma by measured cohesion
#          subject to a hard legibility floor.
#
#          THREE PANELS, because the effect is per-state and unevenly spread:
#            A  per-state cohesion change, diverging bars (the honest summary)
#            B  the gamma each state now uses, against the old rule
#            C  the four biggest movers, before and after, on real geometry
#
# Author:  Claude (Opus 5) for Jorge Roa
# ISS:     -
# Version: gamma-effect v1
# Date:    27-August-2026
#
# Usage:   Rscript plot_gamma_effect.R [old_tag] [new_tag]
# Output:  state_pages/GAMMA_EFFECT_<old>_vs_<new>.png
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(mxmaps)
})
sf_use_s2(FALSE)
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
source("hexlayout_common.R")

MX_CRS <- 6372
# Diverging pair + neutral grey midpoint. Never a hue at the midpoint, and the
# bars carry a direct label so identity is not colour-alone.
COL_UP <- "#1F6FB2"
COL_DN <- "#B2432B"
COL_MID <- "#9AA0A6"
INK <- "grey15"
INK2 <- "grey40"

args <- commandArgs(trailingOnly = TRUE)
OLD <- if (length(args) >= 1) args[1] else "v11_2"
NEW <- if (length(args) >= 2) args[2] else "v11_3"

f_old <- sprintf("data/mxmunicipio_hex_sf_%s.rds", OLD)
f_new <- sprintf("data/mxmunicipio_hex_sf_%s.rds", NEW)
for (f in c(f_old, f_new)) if (!file.exists(f)) stop("missing artifact: ", f)

pm <- readRDS("data/true_pairs_meta.rds")
mx <- build_mxmun_sf()
borders <- readRDS("data/mxstate_borders.rds")
lk <- df_mxmunicipio_2020 |>
  select(state_code, state_abbr) |>
  distinct()

coh_by_state <- function(f) {
  hex <- readRDS(f)
  hm <- st_transform(hex, MX_CRS)
  inter <- st_is_within_distance(hm, hm, dist = 1, sparse = TRUE)
  h_of <- match(mx$region, hex$region)
  pa <- h_of[pm$i]
  pb <- h_of[pm$j]
  ok <- !is.na(pa) & !is.na(pb)
  sv <- rep(NA, nrow(pm))
  sv[ok] <- mapply(function(x, y) y %in% inter[[x]], pa[ok], pb[ok])
  d <- pm
  d$surv <- sv
  d |>
    filter(!is.na(surv)) |>
    group_by(state_code) |>
    summarise(pairs = n(), coh = 100 * mean(surv), .groups = "drop")
}

a <- coh_by_state(f_old) |> rename(coh_old = coh)
b <- coh_by_state(f_new) |> rename(coh_new = coh)
cmp <- a |>
  left_join(b |> select(state_code, coh_new), by = "state_code") |>
  left_join(lk, by = "state_code") |>
  mutate(delta = coh_new - coh_old)

gt <- read.csv("data/gamma_by_state.csv", colClasses = c(state_code = "character"))
szt <- state_size_table(mx)
old_rule <- szt$per_state |> select(state_code, gamma_old = gamma_auto)
cmp <- cmp |>
  left_join(gt |> rename(gamma_new = gamma), by = "state_code") |>
  left_join(old_rule, by = "state_code")

nat_old <- sum(cmp$coh_old / 100 * cmp$pairs) / sum(cmp$pairs) * 100
nat_new <- sum(cmp$coh_new / 100 * cmp$pairs) / sum(cmp$pairs) * 100
cat(sprintf("national cohesion %s %.2f%% -> %s %.2f%%  (%+.2f points)\n",
            OLD, nat_old, NEW, nat_new, nat_new - nat_old))

# ---- panel A: per-state cohesion change ----------------------------------
pa_df <- cmp |>
  mutate(state_abbr = reorder(state_abbr, delta),
         dir = ifelse(abs(delta) < 0.05, "flat",
                      ifelse(delta > 0, "up", "down")))
p_a <- ggplot(pa_df, aes(x = delta, y = state_abbr, fill = dir)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, color = INK, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.1f", delta),
                hjust = ifelse(delta >= 0, -0.15, 1.15)),
            size = 2.9, color = INK) +
  scale_fill_manual(values = c(up = COL_UP, down = COL_DN, flat = COL_MID),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.14, 0.14))) +
  labs(title = "A — cohesion change by state",
       subtitle = sprintf("%s to %s, percentage points of intra-state adjacency kept",
                          OLD, NEW),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey92"),
        axis.text = element_text(color = INK),
        plot.title = element_text(face = "bold", color = INK),
        plot.subtitle = element_text(color = INK2, size = 8.5))

# ---- panel B: gamma, old rule vs chosen ----------------------------------
p_b <- ggplot(cmp, aes(x = gamma_old, y = gamma_new)) +
  geom_abline(slope = 1, intercept = 0, color = "grey75", linetype = "22") +
  geom_point(aes(size = pairs, fill = delta), shape = 21,
             color = "white", stroke = 0.5, alpha = 0.95) +
  ggrepel::geom_text_repel(aes(label = state_abbr), size = 2.6,
                           color = INK, max.overlaps = 40, seed = 1) +
  scale_fill_gradient2(low = COL_DN, mid = COL_MID, high = COL_UP,
                       midpoint = 0, name = "cohesion\nchange (pp)") +
  scale_size_area(max_size = 9, name = "pairs") +
  labs(title = "B — the gamma each state now uses",
       subtitle = paste0("old rule keyed on CELL SIZE (x) against the gamma ",
                         "chosen on measured cohesion (y).\nPoints above the ",
                         "dashed line are states the old rule held back."),
       x = "gamma, shipped rule", y = "gamma, chosen") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", color = INK),
        plot.subtitle = element_text(color = INK2, size = 8.5))

# ---- panel C: the four biggest movers, before and after ------------------
movers <- cmp |>
  arrange(desc(delta)) |>
  head(4)
lay_old <- readRDS(f_old)
lay_new <- readRDS(f_new)
state_sf <- mx |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_make_valid()

one <- function(lay, s, ttl) {
  d <- lay[lay$state_code == s, ]
  bb <- st_bbox(state_sf$geometry[state_sf$state_code == s])
  px <- 0.06 * (bb["xmax"] - bb["xmin"])
  py <- 0.06 * (bb["ymax"] - bb["ymin"])
  ggplot() +
    geom_sf(data = d, fill = "grey93", color = "grey35", linewidth = 0.13) +
    geom_sf(data = borders, fill = NA, color = INK, linewidth = 0.4) +
    coord_sf(xlim = c(bb["xmin"] - px, bb["xmax"] + px),
             ylim = c(bb["ymin"] - py, bb["ymax"] + py), expand = FALSE) +
    labs(title = ttl) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(size = 8.5, color = INK, hjust = 0.5))
}
cells <- list()
for (k in seq_len(nrow(movers))) {
  s <- movers$state_code[k]
  ab <- movers$state_abbr[k]
  cells[[length(cells) + 1]] <- one(
    lay_old, s, sprintf("%s  %s  %.0f%%", ab, OLD, movers$coh_old[k]))
  cells[[length(cells) + 1]] <- one(
    lay_new, s, sprintf("%s  %s  %.0f%%  (%+.1f)", ab, NEW,
                        movers$coh_new[k], movers$delta[k]))
}
p_c <- wrap_plots(cells, nrow = 2, byrow = FALSE) +
  plot_annotation(title = "C — the four biggest movers, before and after")

p <- (p_a | p_b) / p_c +
  plot_layout(heights = c(1, 0.95)) +
  plot_annotation(
    title = sprintf("Re-keying the per-state gamma rule: %s to %s", OLD, NEW),
    subtitle = sprintf(paste0("National cohesion %.2f%% to %.2f%% (%+.2f points). ",
                              "Legibility held as a CONSTRAINT, not a term: every ",
                              "state keeps 0 cells below a quarter of its mean."),
                       nat_old, nat_new, nat_new - nat_old),
    theme = theme(plot.title = element_text(face = "bold", size = 15,
                                            color = INK),
                  plot.subtitle = element_text(color = INK2, size = 9.5)))

fn <- sprintf("state_pages/GAMMA_EFFECT_%s_vs_%s.png", OLD, NEW)
ggsave(fn, p, width = 15, height = 12, dpi = 300, bg = "white")
cat("Wrote", fn, "\n")

write.csv(cmp, "data/gamma_effect_by_state.csv", row.names = FALSE)
cat("Wrote data/gamma_effect_by_state.csv\n")
