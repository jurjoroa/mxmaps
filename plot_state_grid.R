# ---------------------------------------------------------------------------
# Purpose: every state of one layout on a single sheet, each panel labelled
#          with the numbers that decide whether that state is good.
#
#          WHY A CONTACT SHEET. The per-state pages (cls_<tag>_<code>_<abbr>.png)
#          show one state against true territory in detail, but 32 files cannot
#          be compared at a glance, and the thing a reader actually wants to know
#          -- "is this layout uniformly decent or is it carried by a few states"
#          -- is only visible side by side.
#
#          WHY pct_ceil IS ON EVERY PANEL. Raw cohesion is partly a statement
#          about a state's GEOGRAPHY, not its layout: a cell has ~6 usable faces,
#          so a municipio with 22 true neighbours cannot keep more than 6 of them.
#          Oaxaca reads 77% raw and 91% of its own k=6 ceiling. Ranking states on
#          the raw number sends tuning effort at whichever state is densest.
#
# Author:  Claude (Opus 5) for Jorge Roa
# ISS:     -
# Version: state-grid v1
# Date:    28-August-2026
#
# Usage:   Rscript plot_state_grid.R [tag]
# Output:  state_pages/STATE_GRID_<tag>.png
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
INK <- "grey15"
INK2 <- "grey45"

args <- commandArgs(trailingOnly = TRUE)
TAG <- if (length(args) >= 1) args[1] else "v11_3"
f <- sprintf("data/mxmunicipio_hex_sf_%s.rds", TAG)
if (!file.exists(f)) stop("missing artifact: ", f)

lay <- readRDS(f)
pm <- readRDS("data/true_pairs_meta.rds")
mm <- readRDS("data/true_mun_meta.rds")
mx <- build_mxmun_sf()
gam <- read.csv("data/gamma_by_state.csv", colClasses = c(state_code = "character"))
lk <- df_mxmunicipio_2020 |>
  select(state_code, state_abbr, state_name) |>
  distinct()

# per-state cohesion and ceiling, from the one cached pair set
lm_ <- st_transform(lay, MX_CRS)
inter <- st_is_within_distance(lm_, lm_, dist = 1, sparse = TRUE)
h_of <- match(mx$region, lay$region)
pa <- h_of[pm$i]
pb <- h_of[pm$j]
ok <- !is.na(pa) & !is.na(pb)
pm$surv <- NA
pm$surv[ok] <- mapply(function(x, y) y %in% inter[[x]], pa[ok], pb[ok])

stat <- pm |>
  filter(!is.na(surv)) |>
  group_by(state_code) |>
  summarise(pairs = n(), kept = sum(surv), .groups = "drop") |>
  mutate(coh = 100 * kept / pairs)
cap <- mm |>
  group_by(state_code) |>
  summarise(cap6 = 0.5 * sum(cap6), .groups = "drop")
stat <- stat |>
  left_join(cap, by = "state_code") |>
  mutate(cap6 = pmin(cap6, pairs), pct_ceil = 100 * kept / cap6) |>
  left_join(lk, by = "state_code") |>
  left_join(gam, by = "state_code")

# population per 100k WITHIN the state, so a colour means the same thing in
# every panel relative to that state's own distribution
cs <- df_mxmunicipio_2020 |>
  select(region, state_code, pop) |>
  group_by(state_code) |>
  mutate(p_state = pop / sum(pop) * 1e5) |>
  ungroup() |>
  select(region, p_state)
lay <- lay |> left_join(cs, by = "region")

state_sf <- mx |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_make_valid()

panel <- function(s) {
  d <- lay[lay$state_code == s, ]
  st <- stat[stat$state_code == s, ]
  bb <- st_bbox(state_sf$geometry[state_sf$state_code == s])
  px <- 0.04 * (bb["xmax"] - bb["xmin"])
  py <- 0.04 * (bb["ymax"] - bb["ymin"])
  lw <- if (nrow(d) > 200) 0.06 else if (nrow(d) > 60) 0.12 else 0.22
  ggplot() +
    geom_sf(data = d, aes(fill = p_state), color = "grey30", linewidth = lw) +
    geom_sf(data = state_sf[state_sf$state_code == s, ], fill = NA,
            color = INK, linewidth = 0.35) +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, trans = "log10",
                         na.value = "grey90", guide = "none") +
    coord_sf(xlim = c(bb["xmin"] - px, bb["xmax"] + px),
             ylim = c(bb["ymin"] - py, bb["ymax"] + py), expand = FALSE) +
    labs(title = sprintf("%s  n=%d", st$state_abbr, nrow(d)),
         subtitle = sprintf("coh %.0f%%  ceil %.0f%%  γ %.2f",
                            st$coh, st$pct_ceil, st$gamma)) +
    theme_void(base_size = 8) +
    theme(plot.title = element_text(size = 9, face = "bold", colour = INK,
                                    hjust = 0.5, margin = margin(b = 1)),
          plot.subtitle = element_text(size = 7, colour = INK2, hjust = 0.5,
                                       margin = margin(b = 2)),
          plot.margin = margin(3, 3, 3, 3))
}

# worst pct_ceil first: the panels a reader should look at hardest come first
ord <- stat$state_code[order(stat$pct_ceil)]
plots <- lapply(ord, panel)

nat_coh <- 100 * sum(stat$kept) / sum(stat$pairs)
nat_ceil <- 100 * sum(stat$kept) / sum(stat$cap6)
p <- wrap_plots(plots, ncol = 6) +
  plot_annotation(
    title = sprintf("%s — every state, ordered worst-first by %% of its own k=6 ceiling",
                    TAG),
    subtitle = sprintf(paste0("National cohesion %.1f%%, or %.1f%% of the k=6 ceiling. ",
                              "coh = %% of true intra-state adjacencies kept. ",
                              "ceil = %% of the most any ~6-face tessellation could keep, ",
                              "which is what makes states comparable.\n",
                              "γ = that state's area-variation exponent. ",
                              "Colour is inhabitants per 100k within the state (log scale); ",
                              "cell SIZE carries no data and differs 223× between states."),
                      nat_coh, nat_ceil),
    theme = theme(plot.title = element_text(face = "bold", size = 15, colour = INK),
                  plot.subtitle = element_text(size = 9, colour = INK2)))

fn <- sprintf("state_pages/STATE_GRID_%s.png", TAG)
ggsave(fn, p, width = 16, height = 15, dpi = 300, bg = "white", limitsize = FALSE)
cat("Wrote", fn, "\n")
cat(sprintf("national %.2f%% cohesion | %.2f%% of ceiling\n", nat_coh, nat_ceil))
