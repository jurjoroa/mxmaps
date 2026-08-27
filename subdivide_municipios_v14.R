# ---------------------------------------------------------------------------
# Purpose: v14 -- RECURSIVE SUBDIVISION (mosaic / rectangular-cartogram in the
#          Van Kreveld & Speckmann sense, KD-tree flavour) of each real state
#          polygon into one cell per municipio.
#
#          Why this family. Every previous attempt in this repo (v3-v8 hex grid
#          + clip, Hungarian, CVT, adjacency-aware CVT, v9 mesh deformation,
#          v10 uniform lattice, v11-v13 power diagram by damped Newton on the
#          OT dual) is an ITERATIVE optimiser: it can stall, it can leave a
#          cell empty, it needs a coverage-repair pass. Recursive subdivision
#          gets the same hard constraints by CONSTRUCTION, with no solver:
#
#            exact target area  -- each cut is solved by bisection on the
#                                  clipped area, so the ratio is met to
#                                  tolerance and the leaves partition exactly
#            exact containment  -- leaves are intersections of the REAL state
#                                  polygon with half-planes, so the state and
#                                  national outlines are reproduced vertex for
#                                  vertex
#            no gaps / overlaps -- a binary space partition is a partition
#            contiguity         -- each leaf is a convex-cut piece of its
#                                  parent; disconnected only where the state
#                                  polygon itself is (islands, isthmuses)
#
#          The price is shape: cells are orthogonal-ish slabs, not hexagons,
#          so isoperimetric quotient falls. That is the thing this script is
#          written to MEASURE.
#
#          The algorithm, per state:
#            1. domain  = st_union of the state's true municipio polygons
#                         (full multipolygon, islands kept), in EPSG:6372
#            2. targets = true_km2^gamma, normalised to the domain area, with
#                         the same 3-round [0.45, 2.2] x share clamp as v8/v11
#                         and the same per-state gamma (auto rule + the v8.8
#                         override table) -- so the area brief is identical to
#                         v8.8 and v11 and the comparison is apples to apples
#            3. recurse(P, S):
#                 |S| == 1            -> P is that municipio's cell
#                 else                -> axis = longer side of bbox(P)
#                                        sort S by centroid coord on that axis
#                                        k = split index nearest half the
#                                            target mass of S
#                                        r = target mass of S[1..k] / total
#                                        c = bisection root of
#                                            area(P & {axis <= c}) = r*area(P)
#                                        recurse both halves
#            4. leaves -> sf(region, state_code, geometry), back to EPSG:4326
#
# Author:  Claude (Opus 5) for Jorge Roa
# ISS:     -
# Version: v14
# Date:    2026-08-26
#
# Usage:   Rscript subdivide_municipios_v14.R [state_code ...]
#          V14_TOL=<rel>   area tolerance per cut (default 1e-4 of the parent)
#          V14_MODE=bbox|spread|pca   cut-direction rule. bbox (default) is
#                          the longer side of the sub-polygon bbox; spread
#                          the wider municipio-centroid range; pca the
#                          centroid cloud's principal axis, i.e. a straight
#                          cut at an arbitrary angle.
#          V14_OUT=<path>  output path override
# Output:  data/mxmunicipio_hex_sf_v14.rds  (only the states actually run)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(mxmaps)
})
sf_use_s2(FALSE)
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

args <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
ONLY <- if (length(args) >= 1 && !any(is.na(args))) {
  sprintf("%02d", args)
} else NULL
.env <- function(nm, d) {
  v <- Sys.getenv(nm); if (nzchar(v)) as.numeric(v) else d
}
TOL <- .env("V14_TOL", 1e-4)
MAXIT <- as.integer(.env("V14_MAXIT", 40))
MODE <- { v <- Sys.getenv("V14_MODE"); if (nzchar(v)) v else "bbox" }
OUT <- { v <- Sys.getenv("V14_OUT")
         if (nzchar(v)) v else "data/mxmunicipio_hex_sf_v14.rds" }

# ---- true municipio polygons (same construction as the v3-v13 line) -------
data("mxmunicipio.map")
df <- mxmunicipio.map
regions <- unique(df$region)
polys <- lapply(regions, function(r) {
  sub <- df[df$region == r, ]
  groups <- unique(sub$group)
  rings <- lapply(groups, function(g) {
    pts <- sub[sub$group == g, ]
    pts <- pts[order(pts$order), ]
    coords <- cbind(pts$long, pts$lat)
    if (nrow(coords) < 4) return(NULL)
    if (!all(coords[1, ] == coords[nrow(coords), ]))
      coords <- rbind(coords, coords[1, ])
    list(coords = coords, hole = any(pts$hole))
  })
  rings <- Filter(Negate(is.null), rings)
  if (length(rings) == 0) return(st_polygon())
  exterior <- lapply(rings[!sapply(rings, `[[`, "hole")], `[[`, "coords")
  holes    <- lapply(rings[sapply(rings,  `[[`, "hole")], `[[`, "coords")
  st_polygon(c(exterior, holes))
})
mxmun_sf <- st_make_valid(st_sf(region = regions,
                                geometry = st_sfc(polys, crs = 4326)))
mxmun_sf$state_code <- substr(mxmun_sf$region, 1, 2)

mxmun_m <- st_transform(mxmun_sf, 6372)
mxmun_m$true_km2 <- as.numeric(st_area(mxmun_m)) / 1e6
cn <- st_coordinates(st_centroid(mxmun_m, of_largest_polygon = TRUE))
mxmun_m$cx <- cn[, 1]
mxmun_m$cy <- cn[, 2]

# ---- gamma: identical rule to v8.8 / v11 ---------------------------------
state_areas <- mxmun_m |>
  st_drop_geometry() |>
  group_by(state_code) |>
  summarise(state_km2 = sum(true_km2), n = dplyr::n(), .groups = "drop") |>
  mutate(cs_km = sqrt(state_km2 / n / 0.866))
cs_national <- sqrt(sum(state_areas$state_km2) / sum(state_areas$n) / 0.866)
state_areas$gamma <- pmax(pmin(0.18 * (state_areas$cs_km / cs_national - 0.75),
                               0.45), 0)
gamma_override <- c("02" = 0.35, "03" = 0.45, "05" = 0.30, "08" = 0.30,
                    "26" = 0.25, "10" = 0.25, "07" = 0.12, "04" = 0.15,
                    "24" = 0.12, "01" = 0.10, "09" = 0.10, "06" = 0.12)
state_areas$gamma[match(names(gamma_override), state_areas$state_code)] <-
  gamma_override

# ---- geometry helpers ----------------------------------------------------
poly_only <- function(g) {
  # keep only the polygonal part of an intersection result; an empty result
  # is a legitimate outcome when the trial cut lands off the sub-polygon
  if (length(g) == 0) return(st_sfc(st_polygon(), crs = st_crs(g)))
  ty <- as.character(st_geometry_type(g))
  if (ty %in% c("GEOMETRYCOLLECTION")) {
    g <- tryCatch(st_collection_extract(g, "POLYGON"), error = function(e) g)
  }
  ty <- as.character(st_geometry_type(g))
  if (!ty %in% c("POLYGON", "MULTIPOLYGON")) return(st_sfc(st_polygon(),
                                                           crs = st_crs(g)))
  g
}
a_of <- function(g) {
  if (length(g) == 0 || length(g[[1]]) == 0) return(0)
  sum(as.numeric(st_area(g)))
}

# Generalised half-plane {x : <x,u> <= c} clipped to a box big enough to
# cover P. u = (1,0) or (0,1) recovers the axis-aligned KD-tree cut; MODE
# "pca" lets u be the principal axis of the municipio centroid cloud, which
# is still a straight cut, so every guarantee is untouched.
half_box <- function(u, c_val, side, R, ctr) {
  # ctr anchors the rectangle on the domain: EPSG:6372 northings are ~2.5e6,
  # so a box built around the coordinate origin would miss the state.
  v <- c(-u[2], u[1])
  s <- if (side == "lo") -1 else 1
  base <- c_val * u + sum(ctr * v) * v
  a <- base - R * v
  b <- base + R * v
  a2 <- a + s * 2 * R * u
  b2 <- b + s * 2 * R * u
  st_polygon(list(rbind(a, b, b2, a2, a)))
}

cut_poly <- function(P, u, c_val, side, R, crs, ctr) {
  hb <- st_sfc(half_box(u, c_val, side, R, ctr), crs = crs)
  poly_only(st_intersection(P, hb))
}

# Bisection with a false-position accelerator, safeguarded so the bracket
# always contracts. area(P & {<x,u> <= c}) is monotone in c, so the root is
# unique up to flat stretches (which only occur off the polygon).
solve_cut <- function(P, u, ratio, R, crs, ctr, lo, hi) {
  A <- a_of(P)
  target <- ratio * A
  flo <- -target; fhi <- A - target
  it <- 0L
  best <- NULL; best_err <- Inf
  while (it < MAXIT) {
    it <- it + 1L
    span <- hi - lo
    if (span < 1e-6) break
    cand <- lo + span * flo / (flo - fhi)
    if (!is.finite(cand)) cand <- lo + 0.5 * span
    cand <- min(max(cand, lo + 0.05 * span), hi - 0.05 * span)
    glo <- cut_poly(P, u, cand, "lo", R, crs, ctr)
    fc <- a_of(glo) - target
    if (abs(fc) < best_err) {
      best_err <- abs(fc); best <- list(c = cand, glo = glo)
    }
    if (abs(fc) <= TOL * A) break
    if (fc < 0) { lo <- cand; flo <- fc } else { hi <- cand; fhi <- fc }
  }
  ghi <- cut_poly(P, u, best$c, "hi", R, crs, ctr)
  list(c = best$c, lo = best$glo, hi = ghi, err = best_err / A, it = it)
}

# ---- the recursion -------------------------------------------------------
subdivide_state <- function(s) {
  mm <- mxmun_m[mxmun_m$state_code == s, ]
  n <- nrow(mm)
  crs_m <- st_crs(mm)
  dom <- st_make_valid(st_union(st_geometry(mm)))
  A_dom <- a_of(dom)
  bbd <- st_bbox(dom)
  R_dom <- 1.5 * sqrt(diff(bbd[c(1, 3)])^2 + diff(bbd[c(2, 4)])^2)
  ctr_dom <- c(mean(bbd[c(1, 3)]), mean(bbd[c(2, 4)]))

  gam <- state_areas$gamma[state_areas$state_code == s]
  tgt <- pmax(mm$true_km2, 1e-12)^gam
  tgt <- tgt / sum(tgt) * A_dom
  share <- A_dom / n
  for (r in 1:3) {
    tgt <- pmax(pmin(tgt, 2.2 * share), 0.45 * share)
    tgt <- tgt / sum(tgt) * A_dom
  }

  cells <- vector("list", n)
  worst <- 0; ncut <- 0L; iters <- 0L
  stack <- list(list(P = dom, S = seq_len(n)))
  while (length(stack) > 0) {
    node <- stack[[length(stack)]]; stack[[length(stack)]] <- NULL
    S <- node$S; P <- node$P
    if (length(S) == 1L) { cells[[S]] <- P; next }
    bb <- st_bbox(P)
    cxy <- cbind(mm$cx[S], mm$cy[S])
    u <- if (MODE == "pca" && length(S) > 2L) {
      ev <- eigen(stats::cov(cxy), symmetric = TRUE)$vectors[, 1]
      if (any(!is.finite(ev)) || sum(ev^2) < 0.5) c(1, 0) else ev
    } else if (MODE == "spread") {
      if (diff(range(cxy[, 1])) >= diff(range(cxy[, 2]))) c(1, 0) else c(0, 1)
    } else {
      if (diff(bb[c(1, 3)]) >= diff(bb[c(2, 4)])) c(1, 0) else c(0, 1)
    }
    key <- as.numeric(cxy %*% u)
    ord <- order(key); S <- S[ord]
    cm <- cumsum(tgt[S]); tot <- cm[length(cm)]
    k <- which.min(abs(cm[-length(cm)] - tot / 2))
    ratio <- cm[k] / tot
    pv <- st_coordinates(P)[, 1:2, drop = FALSE]
    pr <- as.numeric(pv %*% u)
    cut <- solve_cut(P, u, ratio, R_dom, crs_m, ctr_dom,
                     min(pr), max(pr))
    worst <- max(worst, cut$err); ncut <- ncut + 1L; iters <- iters + cut$it
    stack[[length(stack) + 1L]] <- list(P = cut$lo, S = S[seq_len(k)])
    stack[[length(stack) + 1L]] <- list(P = cut$hi,
                                        S = S[seq.int(k + 1L, length(S))])
  }

  geom <- do.call(c, lapply(cells, function(g) st_geometry(g)))
  out <- st_sf(region = mm$region, state_code = mm$state_code,
               geometry = geom, crs = crs_m)
  a <- as.numeric(st_area(out))
  list(sf = out, n = n, gamma = gam,
       worst_cut = worst, ncut = ncut, mean_it = iters / max(ncut, 1L),
       area_gap = abs(sum(a) - A_dom) / A_dom,
       tgt_err = max(abs(a / tgt - 1)), empty = sum(a <= 0))
}

# ---- run -----------------------------------------------------------------
state_codes <- sort(unique(mxmun_sf$state_code))
if (!is.null(ONLY)) state_codes <- intersect(state_codes, ONLY)

res <- list()
for (s in state_codes) {
  t0 <- proc.time()[["elapsed"]]
  r <- withCallingHandlers(subdivide_state(s),
        warning = function(w) invokeRestart("muffleWarning"))
  el <- proc.time()[["elapsed"]] - t0
  cat(sprintf(paste0("[", MODE, "] state %s n=%3d gamma=%.2f | ",
                     "cuts=%3d mean_it=%4.1f | ",
                     "worst_cut=%.2e area_gap=%.2e tgt_err=%5.2f%% empty=%d",
                     " | %6.1fs\n"),
              s, r$n, r$gamma, r$ncut, r$mean_it, r$worst_cut, r$area_gap,
              100 * r$tgt_err, r$empty, el))
  res[[s]] <- r$sf
}

all_sf <- do.call(rbind, res)
all_sf <- st_transform(all_sf, 4326)
all_sf <- all_sf[, c("region", "state_code", "geometry")]
saveRDS(all_sf, OUT)
cat("\nwrote", OUT, "with", nrow(all_sf), "cells across",
    length(unique(all_sf$state_code)), "states\n")
