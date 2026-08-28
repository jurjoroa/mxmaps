# ---------------------------------------------------------------------------
# Purpose: score a municipio layout artifact per state and nationally, on the
#          FOUR things that matter -- displacement, area evenness INCLUDING the
#          tail, cohesion, and shape -- plus, new in v2, cohesion stratified by
#          the two properties of the TERRITORY that dominate it.
#
#          WHY THE STRATIFICATION. Raw per-state cohesion partly measures which
#          pairs a state happens to have, not how good the layout is. Measured
#          over v11.2's 6288 pairs: break rate runs 3.5% at max(deg) <= 4 to
#          57.6% at max(deg) >= 11, and 58.1% -> 92.4% across quintiles of
#          shared border length. Restricting to max(deg) <= 6 pairs drops the
#          between-state sd from 7.52 to 5.71 pp and moves Oaxaca from 71.3% to
#          91.4% -- it stops being the worst state. So a version that "improves
#          cohesion" may only have been handed easier pairs.
#
#          pct_ceil is the honest per-state number: a single tessellation cell
#          gets ~6 usable faces, so a municipio with 22 true neighbours cannot
#          keep more than 6 of them. The k = 6 pair ceiling is
#          0.5 * sum(min(deg, 6)) = 90.5% nationally. pct_ceil reports how much
#          of a state's OWN reachable ceiling the layout achieved, which is
#          comparable across states in a way raw cohesion is not.
#
#          Usage: Rscript evaluate_states.R <rds_path>
#
# Author:  Claude (Opus 5) for Jorge Roa
# ISS:     -
# Version: evaluator v2
# Date:    27-August-2026
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
TINY <- 0.25  # "illegible" threshold as a fraction of the state mean cell area

args <- commandArgs(trailingOnly = TRUE)
rds <- if (length(args) >= 1) args[1] else "data/mxmunicipio_hex_sf_v11_2.rds"
hex <- readRDS(rds)
cat("Evaluating:", rds, "\n\n")

# The polygon builder used to be transcribed here as well; it lived in 52 root
# scripts and every bug in this project survived through that duplication.
mxmun_sf <- build_mxmun_sf()

meta_f <- "data/true_pairs_meta.rds"
if (!file.exists(meta_f))
  stop("missing ", meta_f, " -- run: Rscript pair_metadata.R")
pm <- readRDS(meta_f)
mm <- readRDS("data/true_mun_meta.rds")

true_km2 <- mxmun_sf |>
  st_transform(MX_CRS) |>
  mutate(true_km2 = as.numeric(st_area(geometry)) / 1e6) |>
  st_drop_geometry() |>
  select(region, true_km2)

hex_m <- st_transform(hex, MX_CRS)
areas <- as.numeric(st_area(hex_m))
pos <- st_coordinates(st_centroid(hex_m))
tkm <- true_km2$true_km2[match(hex$region, true_km2$region)]
inter <- st_is_within_distance(hex_m, hex_m, dist = 1, sparse = TRUE)
tcn <- st_coordinates(st_centroid(st_transform(mxmun_sf, MX_CRS),
                                  of_largest_polygon = TRUE))
pos_nat <- match(hex$region, mxmun_sf$region)

# Pair survival, vectorised against the cached pair list. The previous version
# rebuilt the pair set with a nested loop over cells and re-derived it per
# state; scoring every artifact against ONE cached pair set is what makes two
# artifacts comparable pair-by-pair.
h_of <- match(mxmun_sf$region, hex$region)
pa <- h_of[pm$i]
pb <- h_of[pm$j]
has_both <- !is.na(pa) & !is.na(pb)
pm$surv <- NA
pm$surv[has_both] <- mapply(function(x, y) y %in% inter[[x]],
                            pa[has_both], pb[has_both])

gam_tab <- state_size_table(mxmun_sf)$per_state

codes <- sort(unique(hex$state_code))
out <- lapply(codes, function(s) {
  idx <- which(hex$state_code == s)
  n <- length(idx)
  share_s <- sum(areas[idx]) / n
  r_eq <- areas[idx] / share_s
  # The area-vs-target reference used to hardcode gamma = 0.35, but the
  # generators set gamma PER STATE, so for 31 of 32 states this column was
  # scoring against a gamma the layout never used. Derive the same rule here.
  # This still cannot see the 12 hand overrides -- the real fix is for
  # generators to ship their targets alongside the artifact.
  gam_s <- gam_tab$gamma_auto[match(s, gam_tab$state_code)]
  if (is.na(gam_s)) gam_s <- 0.35
  r_tg <- areas[idx] / (sum(areas[idx]) * tkm[idx]^gam_s / sum(tkm[idx]^gam_s))
  dx <- pos[idx, 1] - tcn[pos_nat[idx], 1]
  dy <- pos[idx, 2] - tcn[pos_nat[idx], 2]
  disp <- sqrt(dx^2 + dy^2) / 1000

  ps <- pm[pm$state_code == s & !is.na(pm$surv), ]
  tot <- nrow(ps)
  pres <- sum(ps$surv)
  # k = 6 ceiling for this state, capped by the pairs that actually exist.
  # Each municipio can keep at most min(deg, 6) of its neighbours and every
  # kept pair is counted by both endpoints, hence the halving.
  cap6 <- min(tot, 0.5 * sum(mm$cap6[mm$state_code == s]))

  iq <- vapply(idx, function(i) {
    co <- st_coordinates(hex_m$geometry[i])
    co <- co[!is.na(co[, 1]), ]
    if (nrow(co) < 3) return(NA_real_)
    x <- co[, 1]
    y <- co[, 2]
    p2 <- sum(sqrt((x - c(x[-1], x[1]))^2 + (y - c(y[-1], y[1]))^2))
    4 * pi * areas[i] / (p2^2)
  }, 0)

  lo <- ps[ps$deg_max <= 6, ]
  hi <- ps[ps$deg_max > 6, ]
  thin <- ps[ps$border_q == 1, ]
  data.frame(
    state_code = s, n = n,
    disp_km = round(mean(disp)),
    eq_p5 = round(quantile(r_eq, .05), 2),
    eq_p95 = round(quantile(r_eq, .95), 2),
    tg_p5 = round(quantile(r_tg, .05), 2),
    tg_p95 = round(quantile(r_tg, .95), 2),
    eq_min = round(min(r_eq), 3),
    n_tiny = sum(r_eq < TINY),
    pairs = tot,
    cohesion = round(100 * pres / max(tot, 1)),
    pres_n = pres,
    cap6 = round(cap6),
    pct_ceil = round(100 * pres / max(cap6, 1)),
    coh_lo = if (nrow(lo)) round(100 * mean(lo$surv)) else NA_integer_,
    coh_hi = if (nrow(hi)) round(100 * mean(hi$surv)) else NA_integer_,
    n_hi = nrow(hi),
    coh_thin = if (nrow(thin)) round(100 * mean(thin$surv)) else NA_integer_,
    iq = round(mean(iq, na.rm = TRUE), 2))
})
out_df <- do.call(rbind, out)

geom_cols <- c("state_code", "n", "disp_km", "eq_p5", "eq_p95",
               "tg_p5", "tg_p95", "eq_min", "n_tiny", "iq")
coh_cols <- c("state_code", "n", "pairs", "cohesion", "cap6", "pct_ceil",
              "coh_lo", "coh_hi", "n_hi", "coh_thin")
cat("-- geometry --\n")
print(out_df[, geom_cols], row.names = FALSE)
cat("\n-- cohesion, stratified --\n")
print(out_df[, coh_cols], row.names = FALSE)

# National cohesion, two ways. The pair-weighted figure is the real one: it is
# the fraction of ALL true intra-state municipio adjacencies that survive as
# touching cells. The cell-weighted figure is retained only for comparability
# with the v8.3-v8.5 commit messages, which quoted it while labelling it
# "weighted by pairs" -- dense states carry more pairs per cell, so it reads
# optimistically.
pv <- pm[!is.na(pm$surv), ]
cat(sprintf("\nNational cohesion (pair-weighted, exact): %.1f %% (%d/%d)\n",
            100 * mean(pv$surv), sum(pv$surv), nrow(pv)))
cat(sprintf("National cohesion (cell-weighted, legacy): %.1f %%\n",
            sum(out_df$cohesion * out_df$n) / sum(out_df$n)))
cat(sprintf("National k=6 ceiling: %.1f %% | achieved %.1f %% of it\n",
            100 * sum(out_df$cap6) / nrow(pv),
            100 * sum(out_df$pres_n) / sum(out_df$cap6)))

cat("\n-- national cohesion by endpoint degree --\n")
print(pv |>
        mutate(bucket = cut(deg_max, c(0, 4, 6, 8, 10, Inf),
                            labels = c("<=4", "5-6", "7-8", "9-10", "11+"))) |>
        group_by(bucket) |>
        summarise(pairs = n(), cohesion = round(100 * mean(surv), 1),
                  .groups = "drop"))
cat("\n-- national cohesion by shared-border quintile --\n")
print(pv |>
        group_by(border_q) |>
        summarise(pairs = n(), med_km = round(median(border_km), 1),
                  cohesion = round(100 * mean(surv), 1), .groups = "drop"))

# LEGIBILITY TAIL. The p5/p95 area-ratio columns cannot see a handful of
# near-invisible cells, and a version can buy cohesion by creating them --
# which is exactly what v8.8 did (38 cells under 0.25x their state mean against
# 2 for v8.6). A municipio rendered at 1/10 of its state's mean cell area is
# not legible, and legibility is the entire point of the layout, so this is a
# first-class number and not a footnote.
cat(sprintf("\nLegibility tail: %d cells < %.2fx their state mean, %d < 0.10x, smallest %.4fx\n",
            sum(out_df$n_tiny), TINY,
            sum(vapply(codes, function(s) {
              a <- areas[hex$state_code == s]
              sum(a / mean(a) < 0.10)
            }, 0L)),
            min(out_df$eq_min)))
