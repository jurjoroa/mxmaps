# -----------------------------------------------------------------------------
# Purpose: Data-based validation of municipio layout artifacts. Instead of the
#          purely geometric "cohesion" proxy, measure whether the SPATIAL
#          STRUCTURE of real census variables survives the layout: global
#          Moran's I (national + 6 largest states), Spearman rank correlation
#          of local Moran (LISA) values, and the share of municipios whose
#          LISA quadrant is unchanged between the true map and the layout.
#          Finally, test whether cohesion predicts any of that.
# Author:  Claude (exploratory scratch, repo root)
# Version: 1.0
# Date:    2026-08-26
# Usage:   Rscript validate_spatial_structure.R [rds ...]
#          default: data/mxmunicipio_hex_sf_{v8_5,v8_6,v8_8,v11,v12}.rds
# Output:  printed tables only (no files written)
#
# Notes
#  - Moran's I implemented directly (spdep is not installed):
#      I = (n / S0) * (z' W z) / (z' z),  W row-standardised, S0 = sum(W).
#  - Contiguity (queen) is taken from the geometry itself, identically for the
#    true polygons and for every layout: st_is_within_distance(dist = 1 m) in
#    EPSG:6372, self removed. The 1 m tolerance is what evaluate_states.R uses
#    for hex cells; applying it to both sides keeps the comparison symmetric.
#  - Isolated units (no neighbour) get an undefined lag and are dropped from
#    LISA statistics; their count is reported.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(mxmaps)
})
sf_use_s2(FALSE)
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
options(warn = -1)
options(width = 200)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  args <- file.path("data", paste0("mxmunicipio_hex_sf_",
                                   c("v8_5", "v8_6", "v8_8", "v11", "v12"),
                                   ".rds"))
}

# ---- true municipio polygons ------------------------------------------------
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
               lapply(rings[sapply(rings, "[[", "hole")], "[[", "coords")))
})
mx <- st_make_valid(st_sf(region = regions,
                          geometry = st_sfc(polys, crs = 4326)))
mx$state_code <- substr(mx$region, 1, 2)
mx_m <- st_transform(mx, 6372)
mx$true_km2 <- as.numeric(st_area(mx_m)) / 1e6

# ---- real variables ---------------------------------------------------------
data("df_mxmunicipio_2020")
d <- as.data.frame(df_mxmunicipio_2020)
d <- d[match(mx$region, d$region), ]
stopifnot(identical(d$region, mx$region))

vars <- data.frame(
  # 1. how many people per km2 -- the classic density gradient
  log_pop_density = log10(d$pop / mx$true_km2),
  # 2. share speaking an indigenous language (strongly regionalised, south)
  indigenous_pct = 100 * d$indigenous_language / d$pop,
  # 3. share identifying as Afro-Mexican (Costa Chica / Veracruz clusters)
  afromexican_pct = 100 * d$afromexican / d$pop,
  # 4. sex ratio: males per 100 females (migration signature, north/centre)
  sex_ratio = 100 * d$pop_male / d$pop_female,
  # 5. absolute municipio population (urban hierarchy, not a rate)
  log_pop = log10(d$pop)
)
var_names <- names(vars)

# ---- helpers ----------------------------------------------------------------
nb_list <- function(g_sf, tol = 1) {
  gm <- st_transform(g_sf, 6372)
  lst <- st_is_within_distance(gm, gm, dist = tol, sparse = TRUE)
  lapply(seq_along(lst), function(i) lst[[i]][lst[[i]] != i])
}

# row-standardised spatial lag of z; NA where a unit has no neighbour
row_std_lag <- function(nb, z) {
  vapply(seq_along(nb), function(i) {
    j <- nb[[i]]
    if (length(j) == 0) return(NA_real_)
    mean(z[j])
  }, 0)
}

# I = (n / S0) * (z'Wz)/(z'z); with row-standardisation S0 = #non-isolated rows
moran_i <- function(nb, x) {
  keep <- is.finite(x)
  if (sum(keep) < 3) return(NA_real_)
  z <- rep(NA_real_, length(x))
  z[keep] <- x[keep] - mean(x[keep])
  lag <- row_std_lag(lapply(nb, function(j) j[keep[j]]), z)
  ok <- keep & !is.na(lag)
  s0 <- sum(ok)
  if (s0 == 0) return(NA_real_)
  (sum(keep) / s0) * sum(z[ok] * lag[ok]) / sum(z[keep]^2)
}

# local Moran on standardised z; returns list(lisa, lag, z)
lisa <- function(nb, x) {
  keep <- is.finite(x)
  z <- rep(NA_real_, length(x))
  z[keep] <- (x[keep] - mean(x[keep])) / sd(x[keep])
  lag <- row_std_lag(lapply(nb, function(j) j[keep[j]]), z)
  list(lisa = z * lag, lag = lag, z = z)
}

quad <- function(z, lag) {
  q <- rep(NA_character_, length(z))
  ok <- !is.na(z) & !is.na(lag)
  q[ok] <- ifelse(z[ok] > 0,
                  ifelse(lag[ok] > 0, "HH", "HL"),
                  ifelse(lag[ok] > 0, "LH", "LL"))
  q
}

# subset a neighbour list to a set of indices, reindexed 1..length(idx)
nb_subset <- function(nb, idx) {
  map <- rep(NA_integer_, length(nb))
  map[idx] <- seq_along(idx)
  lapply(idx, function(i) {
    j <- map[nb[[i]]]
    j[!is.na(j)]
  })
}

# ---- true-map reference -----------------------------------------------------
cat("Building neighbour graphs (queen contiguity, 1 m tolerance)...\n")
nb_true <- nb_list(mx)
cat(sprintf("  true map: %d units, mean %.2f neighbours, %d isolated\n",
            nrow(mx), mean(lengths(nb_true)), sum(lengths(nb_true) == 0)))

big6 <- names(sort(table(mx$state_code), decreasing = TRUE))[1:6]
st_lab <- setNames(d$state_abbr[match(big6, d$state_code)], big6)
st_idx <- lapply(big6, function(s) which(mx$state_code == s))
names(st_idx) <- big6

true_ref <- lapply(var_names, function(v) {
  x <- vars[[v]]
  l <- lisa(nb_true, x)
  list(i_nat = moran_i(nb_true, x),
       i_state = vapply(big6, function(s) {
         moran_i(nb_subset(nb_true, st_idx[[s]]), x[st_idx[[s]]])
       }, 0),
       lisa = l$lisa,
       q = quad(l$z, l$lag))
})
names(true_ref) <- var_names

# ---- true intra-state adjacency pairs (for the cohesion recomputation) ------
adj_int <- st_intersects(mx, mx)
true_pairs <- do.call(rbind, lapply(seq_len(nrow(mx)), function(i) {
  js <- adj_int[[i]]
  js <- js[js > i & mx$state_code[js] == mx$state_code[i]]
  if (length(js) == 0) return(NULL)
  data.frame(i = i, j = js)
}))
cat(sprintf("  true intra-state adjacency pairs: %d\n", nrow(true_pairs)))

cohesion_of <- function(nb, ord) {
  # ord[k] = row of the layout holding true-map unit k
  key <- new.env(hash = TRUE, size = 4e5)
  for (a in seq_along(nb)) {
    for (b in nb[[a]]) {
      if (b > a) assign(paste0(a, "_", b), TRUE, envir = key)
    }
  }
  ai <- ord[true_pairs$i]
  bj <- ord[true_pairs$j]
  lo <- pmin(ai, bj)
  hi <- pmax(ai, bj)
  hit <- vapply(seq_along(lo), function(k) {
    exists(paste0(lo[k], "_", hi[k]), envir = key, inherits = FALSE)
  }, TRUE)
  100 * mean(hit)
}

# ---- per-layout evaluation --------------------------------------------------
# A control layout is appended: the best artifact's cells with the region
# labels randomly permuted WITHIN each state. It has the same outlines, the
# same cell shapes and the same per-state counts as a real layout, so it fixes
# the floor for "a layout that is valid but carries no spatial information".
ctrl_src <- if (any(grepl("v8_8", args))) {
  grep("v8_8", args, value = TRUE)[1]
} else {
  args[length(args)]
}
lay_names <- c(sub("^mxmunicipio_hex_sf_", "", tools::file_path_sans_ext(
  basename(args))), "SHUFFLED")
args <- c(args, ctrl_src)
n_real <- length(args) - 1L
res <- list()
coh <- setNames(rep(NA_real_, length(args)), lay_names)

for (k in seq_along(args)) {
  hx <- readRDS(args[k])
  if (lay_names[k] == "SHUFFLED") {
    set.seed(20260826)
    for (s in unique(hx$state_code)) {
      i <- which(hx$state_code == s)
      hx$region[i] <- sample(hx$region[i])
    }
  }
  ord <- match(mx$region, hx$region)
  stopifnot(!any(is.na(ord)))
  nb_h_raw <- nb_list(hx)
  # reindex the layout graph into true-map order
  nb_h <- nb_subset(nb_h_raw, ord)
  cat(sprintf("  %-6s: mean %.2f neighbours, %d isolated\n", lay_names[k],
              mean(lengths(nb_h)), sum(lengths(nb_h) == 0)))
  coh[k] <- cohesion_of(nb_h_raw, ord)

  res[[lay_names[k]]] <- lapply(var_names, function(v) {
    x <- vars[[v]]
    l <- lisa(nb_h, x)
    q <- quad(l$z, l$lag)
    tr <- true_ref[[v]]
    ok <- !is.na(tr$lisa) & !is.na(l$lisa)
    list(
      i_nat = moran_i(nb_h, x),
      i_state = vapply(big6, function(s) {
        moran_i(nb_subset(nb_h, st_idx[[s]]), x[st_idx[[s]]])
      }, 0),
      rho = cor(tr$lisa[ok], l$lisa[ok], method = "spearman"),
      q_agree = 100 * mean(tr$q[ok] == q[ok]),
      n_ok = sum(ok)
    )
  })
  names(res[[lay_names[k]]]) <- var_names
}

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d, width = 6)

cat("\n============================================================\n")
cat("TABLE 1. Global Moran's I, national (2469 municipios)\n")
cat("============================================================\n")
t1 <- data.frame(variable = var_names,
                 true = fmt(vapply(true_ref, function(z) z$i_nat, 0)))
for (ln in lay_names) {
  t1[[ln]] <- fmt(vapply(res[[ln]], function(z) z$i_nat, 0))
}
print(t1, row.names = FALSE)
cat("\n% of true Moran's I retained by the layout:\n")
t1r <- data.frame(variable = var_names)
for (ln in lay_names) {
  t1r[[ln]] <- round(100 * vapply(res[[ln]], function(z) z$i_nat, 0) /
                       vapply(true_ref, function(z) z$i_nat, 0))
}
print(t1r, row.names = FALSE)

cat("\n============================================================\n")
cat("TABLE 2. Global Moran's I within the 6 largest states\n")
cat("============================================================\n")
for (s in big6) {
  cat(sprintf("\n-- %s (%s), n = %d municipios\n", st_lab[[s]], s,
              length(st_idx[[s]])))
  t2 <- data.frame(variable = var_names,
                   true = fmt(vapply(true_ref, function(z) z$i_state[[s]], 0)))
  for (ln in lay_names) {
    t2[[ln]] <- fmt(vapply(res[[ln]], function(z) z$i_state[[s]], 0))
  }
  print(t2, row.names = FALSE)
}

cat("\n============================================================\n")
cat("TABLE 3. Local Moran (LISA): Spearman rho, true vs layout\n")
cat("============================================================\n")
t3 <- data.frame(variable = var_names)
for (ln in lay_names) {
  t3[[ln]] <- fmt(vapply(res[[ln]], function(z) z$rho, 0))
}
print(t3, row.names = FALSE)

cat("\n============================================================\n")
cat("TABLE 4. HEADLINE: % of municipios whose LISA quadrant is\n")
cat("         UNCHANGED between the true map and the layout\n")
cat("         (chance level with 4 quadrants of unequal size is\n")
cat("          reported as 'null' -- see notes below)\n")
cat("============================================================\n")
t4 <- data.frame(variable = var_names)
for (ln in lay_names) {
  t4[[ln]] <- round(vapply(res[[ln]], function(z) z$q_agree, 0), 1)
}
# expected agreement if the layout's lag sign were independent of the true one,
# holding each side's marginal quadrant distribution fixed
t4$null <- round(vapply(var_names, function(v) {
  tq <- true_ref[[v]]$q
  tq <- tq[!is.na(tq)]
  p <- table(tq) / length(tq)
  100 * sum(p^2)
}, 0), 1)
print(t4, row.names = FALSE)
cat("\nmean over variables:\n")
mq <- vapply(lay_names, function(ln)
  mean(vapply(res[[ln]], function(z) z$q_agree, 0)), 0)
print(round(mq, 2))

cat("\n============================================================\n")
cat("TABLE 5. Does 'cohesion' predict spatial-structure survival?\n")
cat("============================================================\n")
mrho <- vapply(lay_names, function(ln)
  mean(vapply(res[[ln]], function(z) z$rho, 0)), 0)
mret <- vapply(lay_names, function(ln)
  mean(100 * vapply(res[[ln]], function(z) z$i_nat, 0) /
         vapply(true_ref, function(z) z$i_nat, 0)), 0)
t5 <- data.frame(layout = lay_names,
                 cohesion_pct = round(coh, 1),
                 lisa_quad_agree_pct = round(mq, 2),
                 lisa_spearman_rho = round(mrho, 3),
                 moran_retained_pct = round(mret, 1))
print(t5, row.names = FALSE)

pr <- function(a, b, m) {
  if (length(unique(a)) < 3) return(NA_real_)
  cor(a, b, method = m)
}
r <- seq_len(n_real)  # the SHUFFLED control is excluded from correlations
cat(sprintf("\nn real layouts = %d (SHUFFLED control excluded)\n", n_real))
cat(sprintf("cohesion span %.1f-%.1f pts; quad agreement %.2f-%.2f\n",
            min(coh[r]), max(coh[r]), min(mq[r]), max(mq[r])))
cat(sprintf("cor(cohesion, LISA quadrant agreement)  Pearson  = %+.3f\n",
            pr(coh[r], mq[r], "pearson")))
cat(sprintf("cor(cohesion, LISA quadrant agreement)  Spearman = %+.3f\n",
            pr(coh[r], mq[r], "spearman")))
cat(sprintf("cor(cohesion, LISA Spearman rho)        Pearson  = %+.3f\n",
            pr(coh[r], mrho[r], "pearson")))
cat(sprintf("cor(cohesion, Moran's I retained)       Pearson  = %+.3f\n",
            pr(coh[r], mret[r], "pearson")))
sl <- coef(lm(mq[r] ~ coh[r]))[2]
cat(sprintf("OLS slope: +1 cohesion point buys %+.3f points of quadrant\n",
            sl))
cat(sprintf("           agreement. Full observed cohesion span (%.1f pts)\n",
            max(coh[r]) - min(coh[r])))
span <- max(coh[r]) - min(coh[r])
cat(sprintf("           buys %.2f points.\n", sl * span))
cat(sprintf("Control: SHUFFLED cohesion %.1f%%, quadrant agreement %.2f%%\n",
            coh[["SHUFFLED"]], mq[["SHUFFLED"]]))

# per-variable, so a single variable cannot drive the verdict
cat("\nPer-variable cor(cohesion, quadrant agreement) across layouts:\n")
pv <- data.frame(
  variable = var_names,
  pearson = round(vapply(var_names, function(v)
    pr(coh[r], vapply(lay_names[r],
                      function(ln) res[[ln]][[v]]$q_agree, 0),
       "pearson"), 0), 3),
  spearman = round(vapply(var_names, function(v)
    pr(coh[r], vapply(lay_names[r], function(ln) res[[ln]][[v]]$q_agree, 0),
       "spearman"), 0), 3))
print(pv, row.names = FALSE)

cat("\nNOTES\n")
cat("- z_i is the same on both maps (same data, same municipio), so the\n")
cat("  quadrant can only change through the sign of the neighbour lag.\n")
cat("  Quadrant agreement is therefore exactly: the share of municipios\n")
cat("  whose neighbourhood stays on the same side of the national mean.\n")
cat("- 'null' in Table 4 is sum(p_q^2) over the true quadrant marginals:\n")
cat("  the agreement expected from an independent layout with the same\n")
cat("  quadrant mix. Scores at or below it carry no local information.\n")
cat("- Spearman rho in Table 3 is over LISA values, which are dominated by\n")
cat("  |z_i| (identical across maps); it is therefore an optimistic,\n")
cat("  high-floor measure. Table 4 is the discriminating one.\n")
