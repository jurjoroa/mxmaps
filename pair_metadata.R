# ---------------------------------------------------------------------------
# Purpose: cache the per-PAIR properties of the TRUE municipio geometry that
#          every layout must be scored against, so no evaluator recomputes
#          them and no comparison is confounded by them.
#
#          WHY THIS EXISTS. Measured 2026-08-27 over 6288 true intra-state
#          adjacency pairs of v11.2: whether a pair survives is dominated by
#          two properties of the TERRITORY, not of the layout.
#            - endpoint DEGREE: pairs break 3.5% when max(deg) <= 4 and 57.6%
#              when max(deg) >= 11, a 16.4x gradient. Pairs with max(deg) > 6
#              are 58.6% of all pairs but 81.5% of all broken ones.
#            - SHARED BORDER LENGTH: bottom quintile (median 2.4 km) survives
#              58.1%, top quintile (median 43.2 km) survives 92.4%. It is
#              INDEPENDENT of degree (r = -0.15 against log border length) and
#              carries MORE unique deviance than degree (17.4% vs 10.0%).
#
#          Consequence: raw per-state cohesion partly measures which pairs a
#          state happens to have. The between-state sd falls 7.52 -> 5.71 pp
#          when restricted to max(deg) <= 6 pairs, and Oaxaca goes 71.3% ->
#          91.4% and stops being the worst state. Any claim that layout A beats
#          layout B has to be stratified on these or it is partly a claim about
#          geography.
#
#          These are properties of the true polygons alone, so they are the
#          same for every artifact and are computed exactly once.
#
# Author:  Claude (Opus 5) for Jorge Roa
# ISS:     -
# Version: pair-meta v1
# Date:    27-August-2026
#
# Usage:   Rscript pair_metadata.R
# Output:  data/true_pairs_meta.rds   one row per intra-state adjacency pair
#          data/true_pairs_cross.rds  one row per CROSS-state adjacency pair
#          data/true_mun_meta.rds     one row per municipio
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(mxmaps)
})
sf_use_s2(FALSE)
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
source("hexlayout_common.R")

MX_CRS <- 6372  # Mexico LCC. st_area/st_length THROW on lon/lat here (no lwgeom).

mx <- build_mxmun_sf()
mx_m <- st_transform(mx, MX_CRS)
n <- nrow(mx)
cat(sprintf("municipios: %d\n", n))

# Adjacency from st_intersects, then keep same-state ordered pairs. This is the
# same construction evaluate_states.R uses, so the pair set matches exactly.
adj <- st_intersects(mx, mx)
pairs_df <- do.call(rbind, lapply(seq_len(n), function(i) {
  js <- adj[[i]]
  js <- js[js > i & mx$state_code[js] == mx$state_code[i]]
  if (length(js) == 0) return(NULL)
  data.frame(i = i, j = js)
}))
cat(sprintf("intra-state adjacency pairs: %d\n", nrow(pairs_df)))

# Intra-state degree. The k = 6 face ceiling is a statement about same-state
# neighbours, so cross-state adjacency must NOT be counted here.
deg <- tabulate(c(pairs_df$i, pairs_df$j), nbins = n)

# Shared border length. st_intersection of two touching polygons returns the
# shared boundary; a pair that meets only at a corner yields a POINT and 0 m,
# which is the honest answer -- those 96 pairs survive at 46.9% and are worth
# seeing as their own class rather than being dropped as degenerate.
shared_m <- vapply(seq_len(nrow(pairs_df)), function(k) {
  g <- suppressWarnings(
    st_intersection(mx_m$geometry[pairs_df$i[k]], mx_m$geometry[pairs_df$j[k]]))
  if (length(g) == 0) return(0)
  ln <- suppressWarnings(tryCatch(
    st_collection_extract(g, "LINESTRING"), error = function(e) NULL))
  if (is.null(ln) || length(ln) == 0) return(0)
  sum(as.numeric(st_length(ln)))
}, 0)

pairs_df <- pairs_df |>
  mutate(
    region_i = mx$region[i],
    region_j = mx$region[j],
    state_code = mx$state_code[i],
    deg_i = deg[i],
    deg_j = deg[j],
    deg_max = pmax(deg_i, deg_j),
    deg_min = pmin(deg_i, deg_j),
    border_km = shared_m / 1000)

# Quintiles of border length, computed NATIONALLY so the bucket means the same
# thing in every state (a within-state quintile would make each state its own
# yardstick, which is the mistake the quantile-classed choropleths made).
pairs_df$border_q <- as.integer(cut(
  pairs_df$border_km,
  breaks = quantile(pairs_df$border_km, probs = seq(0, 1, 0.2)),
  include.lowest = TRUE, labels = FALSE))

# ---- CROSS-STATE pairs -----------------------------------------------------
# The seam is a measurement hole, not a modelling choice. Every optimiser in
# this repo works one state at a time, and this cache was intra-state by
# construction, so nothing could see what happens where two states meet.
# Measured on v11.3: 49.4% of true cross-state adjacencies survive against
# 82.9% of intra-state ones -- a 33.5 pp gap -- plus fabricated seam contacts
# where cells touch and municipios do not. A reader at a state line is shown
# neighbours that are not neighbours, and no shipped metric reports it.
cross_df <- do.call(rbind, lapply(seq_len(n), function(i) {
  js <- adj[[i]]
  js <- js[js > i & mx$state_code[js] != mx$state_code[i]]
  if (length(js) == 0) return(NULL)
  data.frame(i = i, j = js)
}))
cat(sprintf("cross-state adjacency pairs: %d\n", nrow(cross_df)))

cross_shared_m <- vapply(seq_len(nrow(cross_df)), function(k) {
  g <- suppressWarnings(
    st_intersection(mx_m$geometry[cross_df$i[k]], mx_m$geometry[cross_df$j[k]]))
  if (length(g) == 0) return(0)
  ln <- suppressWarnings(tryCatch(
    st_collection_extract(g, "LINESTRING"), error = function(e) NULL))
  if (is.null(ln) || length(ln) == 0) return(0)
  sum(as.numeric(st_length(ln)))
}, 0)

# Degree here is TOTAL (both sides of the seam), because a seam pair competes
# for faces against a municipio's intra-state neighbours too.
deg_total <- tabulate(c(pairs_df$i, pairs_df$j, cross_df$i, cross_df$j),
                      nbins = n)
cross_df <- cross_df |>
  mutate(
    region_i = mx$region[i],
    region_j = mx$region[j],
    state_i = mx$state_code[i],
    state_j = mx$state_code[j],
    deg_total_i = deg_total[i],
    deg_total_j = deg_total[j],
    deg_total_max = pmax(deg_total_i, deg_total_j),
    border_km = cross_shared_m / 1000)

saveRDS(cross_df, "data/true_pairs_cross.rds")

mun_df <- data.frame(
  region = mx$region,
  state_code = mx$state_code,
  deg = deg,
  deg_total = deg_total,
  # per-municipio contribution to the k = 6 pair ceiling. NOTE this is a
  # BUDGET, not a proven bound -- see ceiling_check; cells do achieve more than
  # 6 same-state contacts in places.
  cap6 = pmin(deg, 6L))

saveRDS(pairs_df, "data/true_pairs_meta.rds")
saveRDS(mun_df, "data/true_mun_meta.rds")

cat("\n-- degree --\n")
print(summary(deg))
cat(sprintf("municipios with deg > 6: %d (%.1f%%)\n",
            sum(deg > 6), 100 * mean(deg > 6)))
cat(sprintf("national k=6 pair ceiling: %.2f%% (%d of %d pairs)\n",
            100 * 0.5 * sum(pmin(deg, 6)) / nrow(pairs_df),
            round(0.5 * sum(pmin(deg, 6))), nrow(pairs_df)))
cat("\n-- shared border km by quintile --\n")
print(pairs_df |>
        group_by(border_q) |>
        summarise(n = n(),
                  min_km = round(min(border_km), 2),
                  med_km = round(median(border_km), 2),
                  max_km = round(max(border_km), 1),
                  .groups = "drop"))
cat(sprintf("point-contact pairs (0 m shared border): %d\n",
            sum(pairs_df$border_km == 0)))
cat(sprintf("cor(deg_max, log1p(border_km)) = %.3f\n",
            cor(pairs_df$deg_max, log1p(pairs_df$border_km))))
cat("\n-- cross-state --\n")
cat(sprintf("pairs: %d | point-contact: %d | median shared border %.2f km\n",
            nrow(cross_df), sum(cross_df$border_km == 0),
            median(cross_df$border_km)))
cat("\nWrote data/true_pairs_meta.rds, data/true_mun_meta.rds, data/true_pairs_cross.rds\n")
