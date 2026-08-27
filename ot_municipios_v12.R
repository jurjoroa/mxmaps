# ---------------------------------------------------------------------------
# Purpose: v12 — v11's optimal-transport solver plus CENTROID TARGETING, so
#          each cell's CENTROID is driven onto its municipio's true centroid.
#
#          Why v11 was not enough: v11 pins the generating SITES on the true
#          centroids, but the metric measures the CELL centroid, and pinning a
#          site does not pin its cell's centroid -- a large or irregular cell's
#          centroid drifts well away from its site. v11's mean displacement is
#          41 km, identical to v8.6, despite the sites being exactly right.
#
#          The system solved here is exactly determined: n area equations and
#          2n centroid equations in n weights and 2n site coordinates.
#            area(C_i)     = T_i                  <- the WEIGHTS' job (Newton)
#            centroid(C_i) = true_centroid_i      <- the SITES' job
#          Outer loop alternates: Newton the weights to hit the areas, then
#          correct each site by the residual centroid error, warm-starting the
#          weights each time. Unlike v8's fixed-strength `anchor`, this is a
#          correction proportional to the MEASURED error, so it self-limits.
#
#          Why this and not v8's CVT. A power diagram's cell adjacency is the
#          dual of a regular (weighted Delaunay) triangulation of its sites.
#          Delaunay of the TRUE municipio centroids already reproduces 83.1% of
#          real municipio adjacency, versus v8.6's measured 73.2% -- so ~10
#          points are sitting in the site positions, for free. v8 then MOVES the
#          sites (CVT) to equalise areas and breaks those edges.
#
#          It never had to. Aurenhammer-Hoffmann-Aronov: for ANY distinct sites
#          and any positive target areas summing to the domain, weights EXIST
#          making the power-diagram cell areas exactly the targets. So area is
#          the weights' job and the sites can stay on the centroids.
#
#          The dual (Merigot; Kitagawa-Merigot-Thibert):
#            Phi(w) = sum_i [ int_{C_i(w)} (|x-P_i|^2 - w_i) dx + w_i T_i ]
#          is CONCAVE, with
#            d Phi / d w_i = T_i - area(C_i(w))
#            d2Phi/dw_i dw_j = |dC_i ^ dC_j| / (2 |P_i - P_j|)   (i != j)
#            d2Phi/dw_i^2    = -sum_{j!=i} (same)
#          i.e. the Hessian is MINUS a graph Laplacian L. Newton step:
#            L dw = T - area(w)
#          L is singular (Phi is invariant under w -> w + c), so solve on the
#          complement of the constant vector. Damping keeps every cell's area
#          above a floor, which is what makes it converge globally rather than
#          oscillate the way a fixed-rate gradient step does.
#
#          Constraints preserved: real state polygons, real national outline,
#          per-state cell sizes, one cell per municipio. Assignment is the
#          IDENTITY -- no Hungarian, no swap search, no restarts.
#
#          NOTE the inherent trade: displacement cannot reach zero while cell
#          areas are near-equal, because equalisation forces dense clusters to
#          inflate. GAMMA is that dial (0 = equal cells, 1 = true areas), so the
#          useful output is the displacement-vs-legibility curve, not a point.
#
# Author:  Claude (Opus 5) for Jorge Roa
# ISS:     -
# Version: v12
# Date:    2026-08-25
#
# Usage:   Rscript ot_municipios_v12.R [state_code ...]
#          V11_GAMMA=<g>   override the per-state gamma with one global value
#          V12_OUTER=<k>   centroid-targeting outer iterations (default 25)
#          V12_BETA=<b>    centroid correction step (default 0.7)
#          V12_TAG=<name>  output suffix, for gamma sweeps
# Output:  data/mxmunicipio_hex_sf_v12<TAG>.rds or /tmp/v12_states.rds
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(mxmaps)
})
sf_use_s2(FALSE)
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

args <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
ONLY <- if (length(args) >= 1 && !any(is.na(args))) sprintf("%02d", args) else NULL
.env <- function(nm, d) { v <- Sys.getenv(nm); if (nzchar(v)) as.numeric(v) else d }
GAMMA_FIX <- { v <- Sys.getenv("V11_GAMMA"); if (nzchar(v)) as.numeric(v) else NA }
KNN     <- as.integer(.env("V11_KNN", 40))   # sites each cell is clipped against
MAXIT   <- as.integer(.env("V11_MAXIT", 250))
DEBUG   <- nzchar(Sys.getenv("V11_DEBUG"))
TOL     <- .env("V11_TOL", 0.01)             # stop at max |A/T - 1| < TOL
OUTER   <- as.integer(.env("V12_OUTER", 25))  # centroid-targeting iterations
BETA    <- .env("V12_BETA", 0.15)              # centroid correction step.
# 0.10-0.30 all measured equivalent once the best iterate is kept; 0.4
# diverges by outer iteration 4 and 0.7 degrades
# cohesion below v11. The centroid error itself is nearly irreducible (2%
# reduction on Aguascalientes) because equalising cell areas forces dense
# clusters to inflate -- so most of the gain here is the gentle site motion
# finding a better ADJACENCY configuration, not a displacement fix.
TAG     <- Sys.getenv("V12_TAG")

# ---- source polygons (same construction as the v3-v8 line) ----------------
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
  st_polygon(c(lapply(rings[!sapply(rings, `[[`, "hole")], `[[`, "coords"),
               lapply(rings[ sapply(rings, `[[`, "hole")], `[[`, "coords")))
})
mxmun_sf <- st_make_valid(st_sf(region = regions,
                                geometry = st_sfc(polys, crs = 4326)))
mxmun_sf$state_code <- substr(mxmun_sf$region, 1, 2)

mxmun_km2 <- mxmun_sf |>
  st_transform(6372) |>
  mutate(true_km2 = as.numeric(st_area(geometry)) / 1e6) |>
  st_drop_geometry() |> select(region, state_code, true_km2)
state_areas <- mxmun_km2 |>
  group_by(state_code) |>
  summarise(state_km2 = sum(true_km2), n = dplyr::n(), .groups = "drop") |>
  mutate(cs_km = sqrt(state_km2 / n / 0.866))
cs_national <- sqrt(sum(state_areas$state_km2) / sum(state_areas$n) / 0.866)
# same gamma rule as v8, so the comparison is apples to apples
state_areas$gamma <- pmax(pmin(0.18 * (state_areas$cs_km / cs_national - 0.75),
                               0.45), 0)
gamma_override <- c("02"=0.35,"03"=0.45,"05"=0.30,"08"=0.30,"26"=0.25,
                    "10"=0.25,"07"=0.12,"04"=0.15,"24"=0.12,"01"=0.10,
                    "09"=0.10,"06"=0.12)
state_areas$gamma[match(names(gamma_override), state_areas$state_code)] <-
  gamma_override

state_sf <- mxmun_sf |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |> st_make_valid()

# ---- geometry helpers -----------------------------------------------------
poly_area <- function(g) {
  if (inherits(g, "sfc")) g <- g[[1]]
  if (length(g) == 0) return(0)
  if (inherits(g, "MULTIPOLYGON"))
    return(sum(vapply(g, function(p) poly_area(p[[1]]), numeric(1))))
  if (inherits(g, "POLYGON")) g <- g[[1]]
  if (is.null(dim(g))) return(0)
  x <- g[, 1]; y <- g[, 2]
  abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
}
largest_piece <- function(g) {
  if (!inherits(g[[1]], "MULTIPOLYGON")) return(g[[1]])
  pcs <- st_cast(g, "POLYGON")
  pcs[[which.max(vapply(seq_along(pcs), function(j) poly_area(pcs[j]), numeric(1)))]]
}

# Sutherland-Hodgman clip of a CONVEX polygon by a half-plane a.x <= b, keeping
# a label per edge so shared-face lengths (the Hessian entries) fall out.
# lab[i] labels the edge V[i] -> V[i+1]; `newlab` labels the edge cut along the
# clip line.
clip_lab <- function(V, lab, a, b, newlab) {
  n <- nrow(V)
  oV <- matrix(NA_real_, n + 1, 2); oL <- integer(n + 1); m <- 0L
  d <- V %*% a
  for (i in seq_len(n)) {
    j <- if (i == n) 1L else i + 1L
    ci <- d[i]; cj <- d[j]
    if (ci <= b) {
      m <- m + 1L; oV[m, ] <- V[i, ]; oL[m] <- lab[i]
      if (cj > b) {
        t <- (b - ci) / (cj - ci)
        m <- m + 1L; oV[m, ] <- V[i, ] + t * (V[j, ] - V[i, ]); oL[m] <- newlab
      }
    } else if (cj <= b) {
      t <- (b - ci) / (cj - ci)
      m <- m + 1L; oV[m, ] <- V[i, ] + t * (V[j, ] - V[i, ]); oL[m] <- lab[i]
    }
  }
  if (m < 3L) return(NULL)
  list(V = oV[seq_len(m), , drop = FALSE], lab = oL[seq_len(m)])
}

# ---- one power-diagram evaluation ----------------------------------------
# Returns exact clipped areas (the GRADIENT, which must be right) and the
# convex-cell face lengths (the HESSIAN, which only preconditions the step --
# taking them before the state-polygon clip is a deliberate approximation).
build_pd <- function(P, w, knn, sq, state_clip) {
  n <- nrow(P)
  cells <- vector("list", n); ar <- numeric(n)
  faces <- vector("list", n)
  for (i in seq_len(n)) {
    V <- sq; lab <- rep(0L, nrow(sq))
    ok <- TRUE
    for (j in knn[, i]) {
      a <- 2 * (P[j, ] - P[i, ])
      b <- sum(P[j, ]^2) - sum(P[i, ]^2) + w[i] - w[j]
      r <- clip_lab(V, lab, a, b, j)
      if (is.null(r)) { ok <- FALSE; break }
      V <- r$V; lab <- r$lab
    }
    if (!ok) { ar[i] <- 0; faces[[i]] <- numeric(0); next }
    cl <- suppressWarnings(st_intersection(
      st_sfc(st_polygon(list(rbind(V, V[1, ]))), crs = 4326), state_clip))
    if (length(cl) == 0) { ar[i] <- 0; faces[[i]] <- numeric(0); next }
    ar[i] <- poly_area(cl)
    cells[[i]] <- cl[[1]]

    # Face lengths must be taken on the CLIPPED cell -- the Hessian entry is
    # |dC_i ^ dC_j| for the clipped cells. Reading them off the convex cell
    # before the state clip (an earlier version) makes the Hessian meaningless
    # for states whose power cells stick far outside the polygon, and Newton
    # then stalls with tau collapsing to ~1e-8 (observed on Zacatecas).
    #
    # Label each clipped edge without any further geometry calls: an edge lying
    # on the i|j bisector has |m-P_i|^2 - w_i == |m-P_j|^2 - w_j at its
    # midpoint m. Anything matching no j is state-boundary.
    co <- st_coordinates(cl)
    fl <- numeric(0)
    if (nrow(co) >= 3) {
      rid <- if ("L2" %in% colnames(co)) paste(co[, "L1"], co[, "L2"]) else co[, "L1"]
      for (rg in unique(rid)) {
        M <- co[rid == rg, c("X", "Y"), drop = FALSE]
        k <- nrow(M)
        if (k < 3) next
        if (all(M[1, ] == M[k, ])) { M <- M[-k, , drop = FALSE]; k <- k - 1L }
        if (k < 3) next
        nx <- c(2:k, 1L)
        mid <- (M + M[nx, , drop = FALSE]) / 2
        len <- sqrt(rowSums((M[nx, , drop = FALSE] - M)^2))
        pdi <- rowSums((mid - matrix(P[i, ], k, 2, byrow = TRUE))^2) - w[i]
        js <- knn[, i]
        D <- vapply(js, function(j)
          abs(rowSums((mid - matrix(P[j, ], k, 2, byrow = TRUE))^2) - w[j] - pdi),
          numeric(k))
        if (k == 1L) D <- matrix(D, nrow = 1L)
        bestj <- js[apply(D, 1, which.min)]
        bestd <- apply(D, 1, min)
        tolm <- 1e-6 * max(pmax(abs(pdi), 1e-12))
        keep <- bestd <= tolm
        if (any(keep)) {
          agg <- tapply(len[keep], bestj[keep], sum)
          for (nmj in names(agg))
            fl[nmj] <- (if (is.na(fl[nmj])) 0 else fl[nmj]) + as.numeric(agg[nmj])
        }
      }
    }
    faces[[i]] <- fl[!is.na(fl)]
  }
  list(cells = cells, ar = ar, faces = faces)
}

state_codes <- sort(unique(mxmun_sf$state_code))
if (!is.null(ONLY)) state_codes <- intersect(state_codes, ONLY)
results <- vector("list", length(state_codes)); names(results) <- state_codes

for (s in state_codes) {
  muns <- mxmun_sf[mxmun_sf$state_code == s, ]
  n <- nrow(muns)
  sp <- state_sf$geometry[state_sf$state_code == s]
  pcs <- st_cast(sp, "POLYGON")
  pa <- vapply(seq_along(pcs), function(j) poly_area(pcs[j]), numeric(1))
  state_clip <- st_sfc(pcs[[which.max(pa)]], crs = 4326)
  A_dom <- poly_area(state_clip)

  # SITES PINNED AT THE TRUE CENTROIDS -- never moved.
  # Take the centroid in the projected CRS: of_largest_polygon = TRUE calls
  # st_area internally, which needs lwgeom on lon/lat multipart geometries
  # (fine on Aguascalientes, fails on Baja California's islands).
  P <- st_coordinates(st_transform(
         st_centroid(st_transform(muns, 6372), of_largest_polygon = TRUE),
         4326))[, 1:2, drop = FALSE]
  # nudge any duplicate sites apart (power diagram needs distinct sites)
  dup <- duplicated(round(P, 10))
  if (any(dup)) P[dup, ] <- P[dup, ] + 1e-7

  gam <- if (is.na(GAMMA_FIX)) state_areas$gamma[state_areas$state_code == s] else GAMMA_FIX
  ta <- vapply(seq_len(n), function(i) poly_area(muns$geometry[i]), numeric(1))
  wt <- pmax(ta, 1e-12) ^ gam
  tgt <- wt / sum(wt) * A_dom
  share <- A_dom / n
  for (r in 1:3) {                       # same clamping as v8
    tgt <- pmax(pmin(tgt, 2.2 * share), 0.45 * share)
    tgt <- tgt / sum(tgt) * A_dom
  }

  bb <- st_bbox(state_clip)
  pad <- 0.6 * max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"])
  sq <- matrix(c(bb["xmin"] - pad, bb["ymin"] - pad,
                 bb["xmax"] + pad, bb["ymin"] - pad,
                 bb["xmax"] + pad, bb["ymax"] + pad,
                 bb["xmin"] - pad, bb["ymax"] + pad), ncol = 2, byrow = TRUE)
  dm <- as.matrix(dist(P))
  kk <- min(KNN, n - 1L)
  knn <- if (n > 1) apply(dm, 1, function(r) order(r)[2:(kk + 1)]) else matrix(integer(0), 0, 1)
  if (is.null(dim(knn))) knn <- matrix(knn, nrow = 1)

  # cell centroid of each power cell, for the outer correction
  cell_ctr <- function(cells) {
    out <- matrix(NA_real_, n, 2)
    for (i in seq_len(n)) {
      if (is.null(cells[[i]])) next
      cc <- suppressWarnings(st_coordinates(
        st_centroid(st_sfc(cells[[i]], crs = 4326))))
      out[i, ] <- cc[1, 1:2]
    }
    out
  }

  w <- rep(0, n)
  P_tgt <- P                       # the true centroids: what cells must land on
  b <- build_pd(P, w, knn, sq, state_clip)
  disp0 <- NA_real_
  # keep the best iterate by mean centroid error, so the result is the best
  # configuration found rather than wherever the loop happened to stop
  best <- list(P = P, w = w, b = b, dsp = Inf)
  # cost per outer iteration is a full Newton solve, so cap the count for large
  # states. The best iterate is almost always found early, so this costs little.
  n_outer <- max(8L, min(OUTER, as.integer(ceiling(4000 / n))))
  for (outer in seq_len(n_outer)) {
    hist <- c()
    for (it in seq_len(MAXIT)) {
      ar <- b$ar
      res <- tgt - ar
      live <- ar > 0
      err <- max(abs(ar[live] / tgt[live] - 1))
      # Merit for the line search is the L2 residual ||A - T||, NOT max|A/T-1|.
      # With the max-norm a single pathological cell vetoes any step that helps
      # all the others, which deadlocked the iteration at ~256% on Zacatecas.
      rn <- sqrt(sum((ar - tgt)^2))
      hist <- c(hist, err)
      if (all(live) && err < TOL) break

      # ---- Hessian = -L, L a graph Laplacian with c_ij = |face| / (2|Pi-Pj|)
      L <- matrix(0, n, n)
      for (i in seq_len(n)) {
        fl <- b$faces[[i]]
        if (!length(fl)) next
        js <- as.integer(names(fl))
        cij <- fl / (2 * pmax(dm[i, js], 1e-12))
        L[i, js] <- L[i, js] - cij
        L[i, i] <- L[i, i] + sum(cij)
      }
      L <- (L + t(L)) / 2                       # symmetrise
      # empty cells have no faces: give them a self term so the system is solvable
      dead <- which(!live | diag(L) <= 0)
      if (length(dead)) for (i in dead) L[i, i] <- max(L[i, i], mean(diag(L)[diag(L) > 0], na.rm = TRUE))
      # L is singular on the constant vector -> pin the mean to zero
      reg <- mean(abs(diag(L))) * 1e-8
      dw <- tryCatch(solve(L + matrix(1 / n, n, n) + diag(reg, n), res),
                     error = function(e) NULL)
      if (is.null(dw) || any(!is.finite(dw))) {
        dw <- res / pmax(diag(L), 1e-12)        # Jacobi fallback
      }
      dw <- dw - mean(dw)

      # ---- damped Newton: keep every cell above an area floor (KMT)
      floor_a <- 0.5 * min(min(ar[live]), min(tgt))
      tau <- 1
      accepted <- FALSE
      for (bt in 1:25) {
        wn <- w + tau * dw
        wn <- wn - mean(wn)
        bn <- build_pd(P, wn, knn, sq, state_clip)
        if (all(bn$ar > floor_a)) {
          # KMT damping: keep every cell above an area floor AND require a
          # genuine decrease of the L2 residual.
          rnn <- sqrt(sum((bn$ar - tgt)^2))
          if (rnn < rn) { w <- wn; b <- bn; accepted <- TRUE; break }
        }
        tau <- tau / 2
      }
      if (DEBUG && (it %% 20 == 0 || it <= 5 || !accepted))
        cat(sprintf("    it%-4d maxerr=%7.2f%%  L2=%10.4g  tau=%9.3g  %s\n",
                    it, 100 * err, rn, tau, if (accepted) "ok" else "STALL"))
      if (!accepted) break
    }

    # ---- centroid targeting: move each site by its cell's centroid error
    cc <- cell_ctr(b$cells)
    ok2 <- !is.na(cc[, 1])
    err_v <- P_tgt[ok2, , drop = FALSE] - cc[ok2, , drop = FALSE]
    dsp <- mean(sqrt(rowSums(err_v^2)))
    if (is.na(disp0)) disp0 <- dsp
    if (dsp < best$dsp) best <- list(P = P, w = w, b = b, dsp = dsp)
    if (DEBUG) cat(sprintf("    outer%-3d mean centroid err = %.5f deg\n", outer, dsp))
    if (outer == n_outer || dsp < 1e-5) break
    # Backtrack the outer step. Taking the full correction with STALE weights
    # can empty a cell, and an earlier version simply gave up when that
    # happened -- which broke the loop at outer=1 for 13 of 32 states, i.e.
    # they got no centroid targeting at all.
    okmove <- FALSE; step <- 1
    for (bt in 1:7) {
      Pn <- P
      Pn[ok2, ] <- P[ok2, , drop = FALSE] + step * BETA * err_v
      # a site must stay inside the state, else its cell can vanish
      inside <- lengths(st_intersects(
        st_sfc(lapply(seq_len(n), function(i) st_point(Pn[i, ])), crs = 4326),
        state_clip)) > 0
      Pn[!inside, ] <- P[!inside, ]
      bn <- build_pd(Pn, w, knn, sq, state_clip)
      if (all(bn$ar > 0)) { okmove <- TRUE; break }
      step <- step / 2
    }
    if (!okmove) break
    P <- Pn; b <- bn
    dm <- as.matrix(dist(P))
    knn <- apply(dm, 1, function(r) order(r)[2:(kk + 1)])
    if (is.null(dim(knn))) knn <- matrix(knn, nrow = 1)
  }

  P <- best$P; w <- best$w; b <- best$b
  ar <- b$ar
  rr <- ar[ar > 0] / tgt[ar > 0]

  # ---- cohesion: assignment is the IDENTITY (cell i is municipio i)
  keep <- which(!vapply(b$cells, is.null, TRUE))
  gc_ <- st_transform(st_sfc(b$cells[keep], crs = 4326), 6372)
  cadj <- st_is_within_distance(gc_, gc_, dist = 1, sparse = TRUE)
  ints <- st_intersects(muns, muns)
  pres <- 0L; tot <- 0L
  for (i in seq_len(n)) {
    js <- ints[[i]]; js <- js[js > i]
    for (j in js) {
      tot <- tot + 1L
      a <- match(i, keep); bq <- match(j, keep)
      if (!is.na(a) && !is.na(bq) && bq %in% cadj[[a]]) pres <- pres + 1L
    }
  }
  cat(sprintf("state %s n=%3d gamma=%.2f | outer=%2d | area p5=%.2f p95=%.2f maxerr=%5.1f%% | cells %d/%d | COHESION %5.1f%%\n",
              s, n, gam, outer, quantile(rr, .05), quantile(rr, .95),
              100 * max(abs(rr - 1)), length(keep), n, 100 * pres / max(tot, 1)))

  geoms <- lapply(seq_len(n), function(i) {
    if (!is.null(b$cells[[i]])) return(largest_piece(st_sfc(b$cells[[i]], crs = 4326)))
    st_buffer(st_sfc(st_point(P[i, ]), crs = 4326), 0.05 * sqrt(tgt[i]))[[1]]
  })
  results[[s]] <- st_sf(region = muns$region, state_code = muns$state_code,
                        geometry = st_sfc(geoms, crs = 4326))
}

out <- do.call(rbind, results[!vapply(results, is.null, TRUE)])
out <- out[order(out$region), ]
if (is.null(ONLY)) {
  saveRDS(out, sprintf("data/mxmunicipio_hex_sf_v12%s.rds", TAG))
  cat(sprintf("Saved data/mxmunicipio_hex_sf_v12%s.rds (%d cells)\n", TAG, nrow(out)))
} else {
  saveRDS(out, "/tmp/v12_states.rds")
  cat(sprintf("Saved /tmp/v12_states.rds (%d cells)\n", nrow(out)))
}
