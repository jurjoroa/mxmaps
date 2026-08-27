# ------------------------------------------------------------------------------
# Purpose: Quantify the national-scale ink bias of the municipio hex layout.
#          Because real state outlines are fixed, a state's total rendered area
#          cannot change; only the split of that area among its municipios can.
#          This script measures how far each state's ink share is from its
#          municipio-count share, for the v8.8 layout and for the true
#          territory map, and decomposes the inequality into between-state and
#          within-state parts.
# Author : Jorge Roa (with Claude Code)
# Version: 1.0
# Date   : 2026-08-26
# Usage  : Rscript quantify_national_bias.R [rds_path]
# Output : console tables + state_pages/national_bias_table.csv
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(mxmaps)
})
sf_use_s2(FALSE)

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
args <- commandArgs(trailingOnly = TRUE)
rds <- if (length(args) >= 1) args[1] else "data/mxmunicipio_hex_sf_v8_8.rds"

hex <- readRDS(rds)
cat("Layout:", rds, " cells:", nrow(hex), "\n\n")

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

mx_m <- st_transform(mx, 6372)
true_km2 <- as.numeric(st_area(mx_m)) / 1e6
hex_m <- st_transform(hex, 6372)
hex_km2 <- as.numeric(st_area(hex_m)) / 1e6

data("df_mxmunicipio_2020")
cen <- df_mxmunicipio_2020

per_mun <- data.frame(region = mx$region, state_code = mx$state_code,
                      true_km2 = true_km2) |>
  left_join(data.frame(region = hex$region, hex_km2 = hex_km2),
            by = "region") |>
  left_join(cen |> select(region, pop, state_abbr), by = "region")
stopifnot(!any(is.na(per_mun$hex_km2)), !any(is.na(per_mun$pop)))

n_mun <- nrow(per_mun)
tot_true <- sum(per_mun$true_km2)
tot_hex <- sum(per_mun$hex_km2)
tot_pop <- sum(per_mun$pop)
cat(sprintf("Total area: true %.0f km2   hex %.0f km2   (hex/true %.4f)\n",
            tot_true, tot_hex, tot_hex / tot_true))

# --- per-state shares ----------------------------------------------------------
st <- per_mun |>
  group_by(state_code) |>
  summarise(abbr = first(state_abbr), n = n(),
            true_km2 = sum(true_km2), hex_km2 = sum(hex_km2),
            pop = sum(pop), .groups = "drop") |>
  mutate(
    share_n = n / n_mun,
    share_true = true_km2 / tot_true,
    share_hex = hex_km2 / tot_hex,
    share_pop = pop / tot_pop,
    km2_per_cell_true = true_km2 / n,
    km2_per_cell_hex = hex_km2 / n,
    ink_ratio_hex = share_hex / share_n,
    ink_ratio_true = share_true / share_n
  ) |>
  arrange(desc(ink_ratio_hex))

cat("\n=== PER-STATE INK SHARE vs MUNICIPIO-COUNT SHARE (v8.8) ===\n")
cat(sprintf("%-5s %-6s %5s %8s %8s %8s %9s %8s %8s\n",
            "code", "abbr", "n", "shr_n%", "shrInk%", "shrPop%",
            "km2/cell", "ink/n", "true/n"))
for (i in seq_len(nrow(st))) {
  r <- st[i, ]
  cat(sprintf("%-5s %-6s %5d %8.2f %8.2f %8.2f %9.0f %8.2f %8.2f\n",
              r$state_code, r$abbr, r$n, 100 * r$share_n, 100 * r$share_hex,
              100 * r$share_pop, r$km2_per_cell_hex, r$ink_ratio_hex,
              r$ink_ratio_true))
}

f <- function(x, y) c(pearson = cor(x, y), spearman = cor(x, y, method = "spearman"))
cat("\n=== CORRELATIONS ACROSS THE 32 STATES ===\n")
cat(sprintf("ink share (hex)  ~ municipio-count share : r=%.3f  rho=%.3f\n",
            f(st$share_hex, st$share_n)[1], f(st$share_hex, st$share_n)[2]))
cat(sprintf("ink share (true) ~ municipio-count share : r=%.3f  rho=%.3f\n",
            f(st$share_true, st$share_n)[1], f(st$share_true, st$share_n)[2]))
cat(sprintf("ink share (hex)  ~ population share      : r=%.3f  rho=%.3f\n",
            f(st$share_hex, st$share_pop)[1], f(st$share_hex, st$share_pop)[2]))
cat(sprintf("ink share (true) ~ population share      : r=%.3f  rho=%.3f\n",
            f(st$share_true, st$share_pop)[1], f(st$share_true, st$share_pop)[2]))

tvd <- function(a, b) sum(abs(a - b)) / 2
cat("\n=== HOW FAR IS INK FROM ONE-CELL-ONE-VOICE? (TVD, 0=perfect) ===\n")
cat(sprintf("state level, hex  vs count share : %.4f  (%.1f%% of ink misallocated)\n",
            tvd(st$share_hex, st$share_n), 100 * tvd(st$share_hex, st$share_n)))
cat(sprintf("state level, true vs count share : %.4f  (%.1f%%)\n",
            tvd(st$share_true, st$share_n), 100 * tvd(st$share_true, st$share_n)))
cat(sprintf("state level, hex  vs pop share   : %.4f\n",
            tvd(st$share_hex, st$share_pop)))
cat(sprintf("state level, true vs pop share   : %.4f\n",
            tvd(st$share_true, st$share_pop)))
eq <- rep(1 / n_mun, n_mun)
cat(sprintf("\nMUNICIPIO level, hex  vs equal    : %.4f  (%.1f%%)\n",
            tvd(per_mun$hex_km2 / tot_hex, eq),
            100 * tvd(per_mun$hex_km2 / tot_hex, eq)))
cat(sprintf("MUNICIPIO level, true vs equal    : %.4f  (%.1f%%)\n",
            tvd(per_mun$true_km2 / tot_true, eq),
            100 * tvd(per_mun$true_km2 / tot_true, eq)))

gini <- function(x) {
  x <- sort(x); n <- length(x)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}
cat(sprintf("\nGini of ink per municipio: true %.3f -> hex %.3f\n",
            gini(per_mun$true_km2), gini(per_mun$hex_km2)))

# --- variance decomposition of log(ink per cell) -------------------------------
dec <- function(a, grp) {
  l <- log(a)
  gm <- tapply(l, grp, mean)
  nn <- tapply(l, grp, length)
  between <- sum(nn * (gm[as.character(grp)][!duplicated(grp)] -
                         mean(l))^2) # placeholder, recomputed below
  bw <- sum(nn * (gm - mean(l))^2) / length(l)
  wi <- sum((l - gm[as.character(grp)])^2) / length(l)
  c(total = bw + wi, between = bw, within = wi)
}
dt <- dec(per_mun$true_km2, per_mun$state_code)
dh <- dec(per_mun$hex_km2, per_mun$state_code)
cat("\n=== VARIANCE OF log(ink per municipio), decomposed ===\n")
cat(sprintf("%-6s %8s %8s %8s %10s\n", "map", "total", "between", "within",
            "within%"))
cat(sprintf("%-6s %8.3f %8.3f %8.3f %9.1f%%\n", "true", dt[1], dt[2], dt[3],
            100 * dt[3] / dt[1]))
cat(sprintf("%-6s %8.3f %8.3f %8.3f %9.1f%%\n", "hex", dh[1], dh[2], dh[3],
            100 * dh[3] / dh[1]))
cat(sprintf("between-state variance change: %+.1f%%   within-state: %+.1f%%\n",
            100 * (dh[2] / dt[2] - 1), 100 * (dh[3] / dt[3] - 1)))

# --- the 223x claim, checked ---------------------------------------------------
cat("\n=== CELL SIZE RANGE BETWEEN STATES (hex) ===\n")
o <- st |> arrange(desc(km2_per_cell_hex))
cat(sprintf("largest : %s %.0f km2/cell (n=%d)\n", o$abbr[1],
            o$km2_per_cell_hex[1], o$n[1]))
cat(sprintf("smallest: %s %.0f km2/cell (n=%d)\n", o$abbr[nrow(o)],
            o$km2_per_cell_hex[nrow(o)], o$n[nrow(o)]))
cat(sprintf("ratio   : %.0fx\n",
            o$km2_per_cell_hex[1] / o$km2_per_cell_hex[nrow(o)]))
cat("\n=== WORST OFFENDERS, ink per data point vs national average ===\n")
avg <- tot_hex / n_mun
o2 <- st |> mutate(x = km2_per_cell_hex / avg) |> arrange(desc(x))
for (i in c(1:5, (nrow(o2) - 4):nrow(o2))) {
  cat(sprintf("  %-6s %5d cells  %8.0f km2/cell  %6.2fx national avg\n",
              o2$abbr[i], o2$n[i], o2$km2_per_cell_hex[i], o2$x[i]))
}

# --- within-state legibility gain ---------------------------------------------
w <- per_mun |>
  group_by(state_code) |>
  summarise(abbr = first(state_abbr),
            rng_true = max(true_km2) / min(true_km2),
            rng_hex = max(hex_km2) / min(hex_km2),
            cv_true = sd(true_km2) / mean(true_km2),
            cv_hex = sd(hex_km2) / mean(hex_km2), .groups = "drop")
cat("\n=== WITHIN-STATE max/min cell area (median over 32 states) ===\n")
cat(sprintf("true %.1fx  ->  hex %.1fx\n", median(w$rng_true),
            median(w$rng_hex)))
cat(sprintf("CV: true %.2f -> hex %.2f (median)\n", median(w$cv_true),
            median(w$cv_hex)))

write.csv(st, "state_pages/national_bias_table.csv", row.names = FALSE)
cat("\nWrote state_pages/national_bias_table.csv\n")
