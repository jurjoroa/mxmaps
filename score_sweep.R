# ---------------------------------------------------------------------------
# Purpose: score every artifact of a gamma sweep against the SAME cached true
#          pair set, and pick the per-state gamma that maximises cohesion
#          subject to a hard legibility floor.
#
#          WHY A FLOOR AND NOT A WEIGHTED SCORE. v8.8 bought +5 cohesion points
#          by rendering 38 cells under a quarter of their state mean and 8
#          fabricated placeholder discs -- including Ecatepec, population 1.65M,
#          drawn at 1.40 km2 against its true 157.8 km2. A weighted objective
#          will always sell legibility for cohesion at some exchange rate, so
#          legibility is a CONSTRAINT here, not a term: n_tiny must stay 0 and
#          the smallest cell must stay at or above MIN_EQ of its state mean.
#          Within that feasible set, and only there, take the best cohesion.
#
# Author:  Claude (Opus 5) for Jorge Roa
# ISS:     -
# Version: sweep-scorer v1
# Date:    27-August-2026
#
# Usage:   Rscript score_sweep.R <dir_with_gN.rds> [out_csv]
# Output:  <out_csv> (default data/gamma_frontier_full.csv)
#          data/gamma_by_state.csv   the chosen per-state gamma
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(mxmaps)
})
sf_use_s2(FALSE)
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
source("hexlayout_common.R")

MX_CRS <- 6372
TINY <- 0.25    # illegible: cell area below this fraction of its state mean
MIN_EQ <- 0.45  # the smallest cell must stay at or above this. v11.2 ships
                # 0.449, so this holds the shipped floor rather than inventing
                # a new one.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript score_sweep.R <dir> [out_csv]")
dir_in <- args[1]
out_csv <- if (length(args) >= 2) args[2] else "data/gamma_frontier_full.csv"

pm <- readRDS("data/true_pairs_meta.rds")
mm <- readRDS("data/true_mun_meta.rds")
mx <- build_mxmun_sf()
tcn <- st_coordinates(st_centroid(st_transform(mx, MX_CRS),
                                  of_largest_polygon = TRUE))

files <- list.files(dir_in, pattern = "^g[0-9.]+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no gN.rds files in ", dir_in)
cat(sprintf("scoring %d artifact(s) from %s\n", length(files), dir_in))

score_one <- function(f) {
  g <- as.numeric(sub("^g", "", sub("\\.rds$", "", basename(f))))
  hex <- readRDS(f)
  hex_m <- st_transform(hex, MX_CRS)
  areas <- as.numeric(st_area(hex_m))
  pos <- st_coordinates(st_centroid(hex_m))
  inter <- st_is_within_distance(hex_m, hex_m, dist = 1, sparse = TRUE)
  pos_nat <- match(hex$region, mx$region)

  h_of <- match(mx$region, hex$region)
  pa <- h_of[pm$i]
  pb <- h_of[pm$j]
  ok <- !is.na(pa) & !is.na(pb)
  surv <- rep(NA, nrow(pm))
  surv[ok] <- mapply(function(x, y) y %in% inter[[x]], pa[ok], pb[ok])

  do.call(rbind, lapply(sort(unique(hex$state_code)), function(s) {
    idx <- which(hex$state_code == s)
    r_eq <- areas[idx] / (sum(areas[idx]) / length(idx))
    sel <- pm$state_code == s & !is.na(surv)
    tot <- sum(sel)
    pres <- sum(surv[sel])
    cap6 <- min(tot, 0.5 * sum(mm$cap6[mm$state_code == s]))
    dx <- pos[idx, 1] - tcn[pos_nat[idx], 1]
    dy <- pos[idx, 2] - tcn[pos_nat[idx], 2]
    data.frame(gamma = g, state_code = s, n = length(idx),
               pairs = tot, pres = pres,
               cohesion = 100 * pres / max(tot, 1),
               cap6 = cap6,
               pct_ceil = 100 * pres / max(cap6, 1),
               eq_min = min(r_eq),
               n_tiny = sum(r_eq < TINY),
               disp_km = mean(sqrt(dx^2 + dy^2)) / 1000)
  }))
}

res <- do.call(rbind, lapply(files, function(f) {
  cat("  ", basename(f), "\n")
  score_one(f)
}))
res <- res[order(res$state_code, res$gamma), ]
write.csv(res, out_csv, row.names = FALSE)
cat("Wrote", out_csv, "\n")

# ---- pick the per-state gamma ---------------------------------------------
pick <- res |>
  group_by(state_code) |>
  group_modify(function(d, k) {
    feas <- d[d$n_tiny == 0 & d$eq_min >= MIN_EQ, ]
    # Every state has gamma = 0 available: equal areas cannot create a tiny
    # cell except by clipping, so the feasible set is never empty in practice.
    # Fall back loudly rather than silently picking an illegible layout.
    if (!nrow(feas)) {
      feas <- d[which.max(d$eq_min), ]
      feas$note <- "INFEASIBLE-fallback"
    } else {
      feas$note <- ""
    }
    best <- feas[which.max(feas$cohesion), ]
    base <- d[d$gamma == 0, ]
    data.frame(gamma = best$gamma, n = best$n, pairs = best$pairs,
               coh_g0 = base$cohesion[1], coh_best = best$cohesion,
               gain = best$cohesion - base$cohesion[1],
               eq_min = best$eq_min, n_tiny = best$n_tiny,
               pct_ceil = best$pct_ceil, note = best$note)
  }) |>
  ungroup()

write.csv(pick[, c("state_code", "gamma")], "data/gamma_by_state.csv",
          row.names = FALSE)
cat("\n-- chosen gamma per state (feasible set: n_tiny == 0 and eq_min >=",
    MIN_EQ, ") --\n")
print(as.data.frame(pick |>
  mutate(across(c(coh_g0, coh_best, gain, eq_min, pct_ceil),
                ~round(.x, 2)))), row.names = FALSE)

tot_pairs <- sum(pick$pairs)
cat(sprintf(
  "\nNational cohesion: gamma=0 %.2f%% -> chosen %.2f%%  (+%.2f points)\n",
  100 * sum(pick$coh_g0 / 100 * pick$pairs) / tot_pairs,
  100 * sum(pick$coh_best / 100 * pick$pairs) / tot_pairs,
  100 * sum((pick$coh_best - pick$coh_g0) / 100 * pick$pairs) / tot_pairs))
cat("Wrote data/gamma_by_state.csv\n")
