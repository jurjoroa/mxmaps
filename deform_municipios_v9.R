# ---------------------------------------------------------------------------
# Purpose: v9 municipio layout — equalize municipio areas by DEFORMING the
#          true polygons on a shared-vertex mesh, instead of re-tessellating
#          and assigning (v3-v8.5). Adjacency is preserved exactly because a
#          shared vertex is shared by construction: moving it moves every
#          incident municipio together, so no true neighbour can come apart.
#          State outlines are preserved exactly by pinning boundary vertices.
#
#          Objective (settled 2026-08-25): every municipio equally legible.
#          Area evenness is the deliverable; shape regularity is the price.
#
#          ONE explicit objective, minimised by gradient descent:
#
#            J = sum_i (A_i / T_i - 1)^2  +  MU * sum_e (|e| / |e0| - 1)^2
#                \_______ area term _____/     \___ edge-spring term ____/
#
#          The area term is the deliverable; it is scale-free (a ratio) so
#          every municipio counts equally regardless of size. The spring term
#          is the only regulariser -- it keeps edges near their original
#          length so cells deform rather than collapse. MU is ONE global
#          number. There is deliberately no per-state parameter table: the
#          v8.5 line ended up hardcoding 12 of 32 states plus Edomex-specific
#          constants, which is curve-fitting, not an algorithm.
#
# Author:  Claude (Opus 5) for Jorge Roa
# ISS:     -
# Version: v9
# Date:    2026-08-25
#
# Usage:   Rscript deform_municipios_v9.R [state_code ...]
# Output:  data/mxmunicipio_equalarea_v9.rds (all states)
#          /tmp/v9_states.rds                (subset run)
#
# Notes:   Deformation runs in EPSG:6372 (metres) so "equal area" means equal
#          km2 -- the same thing evaluate_states.R measures. Output is 4326.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(mxmaps)
})
sf_use_s2(FALSE)

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

args <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
ONLY <- if (length(args) >= 1 && !any(is.na(args))) sprintf("%02d", args) else NULL

# ---- tuning: ONE global set, no per-state table ---------------------------
GAMMA   <- 0.00   # 0 = fully equal area; 1 = true area
DENSIFY <- 3L     # extra shared vertices per mesh edge. The source polygons
                  # carry only ~12 vertices per municipio in dense states,
                  # which is too few degrees of freedom to reshape into equal
                  # areas. Subdividing edges adds DOF without moving any
                  # geometry, and midpoints of a shared edge stay shared.
MU      <- 0.00   # edge-spring weight in J. MUST be 0: any weight makes
                  # the spring a barrier (it starts at 0 and can only grow),
                  # so the line search rejects every step by iteration ~4.
                  # J is therefore the pure objective; shape is regularised in
                  # the search direction (MU_DIR) and by the flip guard.
MU_DIR  <- 0.35   # edge-spring weight in the search direction
N_ITER  <- 600
STEP0   <- 0.35   # initial gradient step (relative, adapted per iteration)
CAP_FRAC <- 0.25  # per-vertex step cap as a fraction of shortest incident
                  # edge -- keeps a vertex from crossing a neighbouring edge
AMIN    <- 0.05   # a ring may not shrink below this fraction of its
                  # original area -- prevents the degeneracy that deadlocks
                  # the guard
TOL     <- 0.01   # stop when max |A/T - 1| < TOL

# env overrides (sweeping/diagnostics only; defaults above are the real config)
.envnum <- function(nm, d) {
  v <- Sys.getenv(nm); if (nzchar(v)) as.numeric(v) else d
}
GAMMA    <- .envnum("V9_GAMMA",   GAMMA)
DENSIFY  <- as.integer(.envnum("V9_DENSIFY", DENSIFY))
MU       <- .envnum("V9_MU",      MU)
MU_DIR   <- .envnum("V9_MU_DIR",  MU_DIR)
N_ITER   <- as.integer(.envnum("V9_ITER", N_ITER))
STEP0    <- .envnum("V9_STEP",    STEP0)
CAP_FRAC <- .envnum("V9_CAP",     CAP_FRAC)
AMIN     <- .envnum("V9_AMIN",    AMIN)
DEBUG    <- nzchar(Sys.getenv("V9_DEBUG"))

# ---- source polygons (same construction as the v3-v8.5 line) --------------

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
  holes    <- lapply(rings[ sapply(rings, `[[`, "hole")], `[[`, "coords")
  st_polygon(c(exterior, holes))
})
mxmun_sf <- st_make_valid(st_sf(region = regions,
                                geometry = st_sfc(polys, crs = 4326)))
mxmun_sf$state_code <- substr(mxmun_sf$region, 1, 2)

# ---- mesh -----------------------------------------------------------------

# Build a shared-vertex mesh. Neighbouring municipios in mxmunicipio.map share
# bit-exact vertices (measured: 98% of adjacent Oaxaca pairs share >=2, median
# 4), so deduping coordinates recovers a planar mesh and adjacency becomes a
# property of the data structure rather than something to be optimised for.
build_mesh <- function(mun_sf) {
  g <- st_cast(st_make_valid(mun_sf$geometry), "MULTIPOLYGON")
  co <- st_coordinates(g)                   # X Y L1(ring) L2(poly) L3(feature)
  rkey <- paste(co[, "L3"], co[, "L2"], co[, "L1"])
  ring_id <- match(rkey, unique(rkey))

  # st_coordinates closes every ring (last == first); drop the duplicate
  is_last <- c(ring_id[-1] != ring_id[-length(ring_id)], TRUE)
  co <- co[!is_last, , drop = FALSE]
  ring_id <- ring_id[!is_last]

  key <- paste(co[, "X"], co[, "Y"])
  uk <- unique(key)
  vi <- match(key, uk)
  V <- cbind(co[match(uk, key), "X"], co[match(uk, key), "Y"])

  ord <- order(ring_id)
  ring_id <- ring_id[ord]; vi <- vi[ord]

  rk <- unique(rkey)
  parts <- do.call(rbind, strsplit(rk, " ", fixed = TRUE))
  list(V = V, vi = vi, ring = ring_id,
       ring_mun  = as.integer(parts[, 1]),
       ring_hole = as.integer(parts[, 3]) > 1L,
       n_mun = nrow(mun_sf))
}

# next-row-within-ring index (wrapping)
mk_nxt <- function(ring) {
  rl <- rle(ring)
  ends <- cumsum(rl$lengths)
  starts <- ends - rl$lengths + 1L
  nxt <- seq_along(ring) + 1L
  nxt[ends] <- starts
  nxt
}

# Subdivide every mesh edge into DENSIFY+1 pieces. New vertices land exactly on
# the existing segment, so geometry is unchanged; they are keyed by the
# undirected edge so the two municipios sharing an edge get the SAME new
# vertices and sharing (hence adjacency) survives densification.
densify_mesh <- function(m, k) {
  if (k <= 0) return(m)
  nxt <- mk_nxt(m$ring)
  a <- m$vi; b <- m$vi[nxt]
  lo <- pmin(a, b); hi <- pmax(a, b)
  ekey <- paste(lo, hi)
  ue <- unique(ekey)
  eid <- match(ekey, ue)
  e_lo <- lo[match(ue, ekey)]; e_hi <- hi[match(ue, ekey)]

  nV <- nrow(m$V)
  fr <- seq_len(k) / (k + 1)
  # new vertex ids: block per fraction, so id = nV + (f-1)*n_edges + eid
  newV <- do.call(rbind, lapply(fr, function(f)
    (1 - f) * m$V[e_lo, , drop = FALSE] + f * m$V[e_hi, , drop = FALSE]))
  V2 <- rbind(m$V, newV)
  n_e <- length(ue)
  ids <- vapply(seq_len(k), function(j) nV + (j - 1L) * n_e + eid,
                integer(length(eid)))
  if (length(eid) == 1L) ids <- matrix(ids, nrow = 1L)

  # expand each directed edge a->b into a, then the k interior vertices in the
  # direction of travel (ids run lo->hi, so reverse when a > b)
  fwd <- a < b
  seq_ids <- ids
  seq_ids[!fwd, ] <- ids[!fwd, k:1, drop = FALSE]

  out_vi <- as.vector(t(cbind(a, seq_ids)))
  out_ring <- rep(m$ring, each = k + 1L)
  list(V = V2, vi = out_vi, ring = out_ring,
       ring_mun = m$ring_mun, ring_hole = m$ring_hole, n_mun = m$n_mun)
}

# ---- geometry + objective (no sf calls; these run every iteration) --------

# signed ring areas (shoelace) and municipio areas
ring_signed <- function(m, V, nxt) {
  x1 <- V[m$vi, 1]; y1 <- V[m$vi, 2]
  x2 <- V[m$vi[nxt], 1]; y2 <- V[m$vi[nxt], 2]
  as.numeric(rowsum(x1 * y2 - x2 * y1, m$ring, reorder = FALSE)) / 2
}

mun_area_of <- function(S, coef, ring_mun, n_mun) {
  a <- numeric(n_mun)
  v <- as.numeric(rowsum(coef * S, ring_mun, reorder = TRUE))
  a[sort(unique(ring_mun))] <- v
  a
}

# ---- per-state deformation ------------------------------------------------

state_codes <- sort(unique(mxmun_sf$state_code))
if (!is.null(ONLY)) state_codes <- intersect(state_codes, ONLY)

results <- vector("list", length(state_codes))
names(results) <- state_codes

for (s in state_codes) {
  mun <- mxmun_sf[mxmun_sf$state_code == s, ]
  mun_m <- st_transform(mun, 6372)
  n <- nrow(mun)

  m <- build_mesh(mun_m)
  m <- densify_mesh(m, DENSIFY)
  V <- m$V
  nxt <- mk_nxt(m$ring)
  prv <- integer(length(nxt)); prv[nxt] <- seq_along(nxt)
  nv <- nrow(V)

  # ring_coef folds in hole-vs-exterior and the ring's initial orientation, so
  # municipio area is a plain signed sum. The flip guard keeps these valid.
  S0 <- ring_signed(m, V, nxt)
  ring_coef <- ifelse(m$ring_hole, -1, 1) * sign(S0)
  A0 <- mun_area_of(S0, ring_coef, m$ring_mun, m$n_mun)
  mun_of_row <- m$ring_mun[m$ring]

  # targets: equal share, optionally nudged toward true area by GAMMA
  w_t <- A0 ^ GAMMA
  Tg  <- w_t / sum(w_t) * sum(A0)

  # pin the state boundary so the outline is preserved exactly
  bnd <- st_boundary(st_union(st_geometry(mun_m)))
  vpts <- st_sfc(lapply(seq_len(nv), function(i) st_point(V[i, ])), crs = 6372)
  pinned <- lengths(st_is_within_distance(vpts, bnd, dist = 1)) > 0

  # rest lengths for the spring term, and the per-vertex step cap
  L0 <- sqrt((V[m$vi, 1] - V[m$vi[nxt], 1])^2 +
             (V[m$vi, 2] - V[m$vi[nxt], 2])^2)
  L0 <- pmax(L0, 1e-6)
  cap <- CAP_FRAC * pmin(
    as.numeric(tapply(L0, m$vi, min))[as.character(seq_len(nv))],
    as.numeric(tapply(L0, m$vi[nxt], min))[as.character(seq_len(nv))])
  cap[is.na(cap)] <- CAP_FRAC * stats::median(L0)

  sgn0 <- sign(S0)

  objective <- function(V) {
    S <- ring_signed(m, V, nxt)
    A <- mun_area_of(S, ring_coef, m$ring_mun, m$n_mun)
    L <- sqrt((V[m$vi, 1] - V[m$vi[nxt], 1])^2 +
              (V[m$vi, 2] - V[m$vi[nxt], 2])^2)
    # area-weighted vertex-mean centroid per municipio: cheap, stable, and
    # only used as the centre to scale about
    cx <- numeric(m$n_mun); cy <- numeric(m$n_mun)
    rx <- rowsum(V[m$vi, 1], mun_of_row, reorder = TRUE)
    ry <- rowsum(V[m$vi, 2], mun_of_row, reorder = TRUE)
    ct <- rowsum(rep(1, length(m$vi)), mun_of_row, reorder = TRUE)
    ii <- as.integer(rownames(rx))
    cx[ii] <- rx[, 1] / ct[, 1]; cy[ii] <- ry[, 1] / ct[, 1]
    list(J = sum((A / Tg - 1)^2) + MU * sum((L / L0 - 1)^2),
         A = A, S = S, L = L, cx = cx, cy = cy)
  }

  st <- STEP0
  fails <- 0L
  o <- objective(V)
  cxr <- o$cx; cyr <- o$cy
  cat(sprintf("state %s  n=%3d  vertices=%6d  pinned=%5d (%2.0f%%)  J0=%.3g\n",
              s, n, nv, sum(pinned), 100 * mean(pinned), o$J))

  for (it in seq_len(N_ITER)) {
    A <- o$A; L <- o$L
    if (max(abs(A / Tg - 1)) < TOL) break

    # --- search direction -------------------------------------------------
    # Scaling municipio i about its centroid by sqrt(T_i/A_i) fixes its area
    # EXACTLY, so that displacement is a near-Newton direction for the area
    # term -- far better conditioned than the raw gradient, whose magnitude
    # just tracks edge length. Vertices are shared, so a vertex takes the
    # urgency-weighted mean of its incident municipios' demands. (These
    # demands agree rather than cancel on a shared border: moving a vertex
    # away from a grower's centroid is the same direction as moving it
    # toward a shrinker's.) J below still decides whether the step is taken.
    k <- sqrt(Tg / pmax(A, 1e-9)) - 1
    ki <- k[mun_of_row]
    wi <- abs(ki) + 1e-3
    dx <- wi * ki * (V[m$vi, 1] - cxr[mun_of_row])
    dy <- wi * ki * (V[m$vi, 2] - cyr[mun_of_row])

    acc <- function(val, idx) {
      z <- numeric(nv)
      r <- rowsum(val, idx, reorder = TRUE)
      z[as.integer(rownames(r))] <- r[, 1]
      z
    }
    D <- cbind(acc(dx, m$vi), acc(dy, m$vi))
    wsum <- pmax(acc(wi, m$vi), 1e-9)
    D <- D / wsum

    # edge springs: restore length, keeping cells deforming rather than
    # collapsing into slivers
    ux <- (V[m$vi, 1] - V[m$vi[nxt], 1]) / L
    uy <- (V[m$vi, 2] - V[m$vi[nxt], 2]) / L
    fs <- MU_DIR * (1 - L / L0) * L0   # metres, matching the area term
    D[, 1] <- D[, 1] + acc(fs * ux, m$vi) - acc(fs * ux, m$vi[nxt])
    D[, 2] <- D[, 2] + acc(fs * uy, m$vi) - acc(fs * uy, m$vi[nxt])
    D[pinned, ] <- 0
    D[!is.finite(D)] <- 0

    # --- flip guard, applied PER VERTEX. Clipping each vertex to its own cap
    #     preserves direction and keeps relative magnitudes below the cap.
    #     Scaling the whole field by the worst ratio instead lets one outlier
    #     vertex (a tiny municipio needing ~13x linear growth, sitting among
    #     short edges) throttle every other vertex to zero -- which is exactly
    #     what stalled Oaxaca at ~40 iterations while Aguascalientes converged.
    dm <- sqrt(D[, 1]^2 + D[, 2]^2)
    D <- D * ifelse(dm > cap, cap / pmax(dm, 1e-12), 1)
    if (max(dm) <= 0 || !all(is.finite(D))) break

    # --- accept only if J decreases and no ring flipped orientation
    # --- guard: no ring may flip orientation OR collapse below AMIN of its
    #     original area. Degeneracy is the real failure mode -- once one ring
    #     reaches ~zero area, every step flips it, and an all-or-nothing global
    #     test then lets that single ring freeze all 570 municipios (this is
    #     exactly what stalled Oaxaca at iteration 29, at st=2e-5). So the
    #     response is LOCALISED: damp only the offending rings' vertices and
    #     let the rest of the map keep moving.
    ok <- FALSE; why <- "flip"; nbad <- 0L
    for (bt in 1:10) {
      Dw <- D
      for (tries in 1:4) {
        Vn <- V + st * Dw
        on <- objective(Vn)
        bad <- which(sign(on$S) != sgn0 | abs(on$S) < AMIN * abs(S0))
        if (!length(bad)) break
        Dw[unique(m$vi[m$ring %in% bad]), ] <- 0
      }
      nbad <- length(bad)
      if (nbad == 0L && on$J < o$J) { ok <- TRUE; break }
      why <- if (nbad > 0L) "flip" else "J-up"
      st <- st / 2
    }
    if (DEBUG && (it <= 12 || it %% 50 == 0 || !ok)) {
      Ln <- o$L
      cat(sprintf("    it%-4d J=%-10.5g area=%-10.5g spring=%-10.5g st=%-9.3g bt=%-3d %s maxD=%.1fm medcap=%.1fm\n",
        it, o$J, sum((o$A/Tg - 1)^2), MU * sum((Ln/L0 - 1)^2), st, bt,
        if (ok) "ok" else paste0("STALL:", why, "/", nbad),
        max(sqrt(D[,1]^2 + D[,2]^2)), stats::median(cap)))
    }
    if (!ok) {
      fails <- fails + 1L
      if (fails >= 8L) break
      st <- STEP0
      next
    }
    fails <- 0L
    V <- Vn; o <- on
    cxr <- o$cx; cyr <- o$cy
    st <- min(st * 1.15, 1.0)
  }

  r <- o$A / Tg
  cat(sprintf("  -> iters=%3d  J=%.4g  area ratio p5=%.2f p95=%.2f  max err=%.1f%%\n",
              it, o$J, quantile(r, .05), quantile(r, .95),
              100 * max(abs(r - 1))))

  # ---- rebuild sf from the deformed vertices
  ord <- order(m$ring)
  rl <- rle(m$ring[ord]); vi_o <- m$vi[ord]
  ends <- cumsum(rl$lengths); starts <- ends - rl$lengths + 1L
  by_mun <- split(seq_along(rl$values), m$ring_mun[rl$values])
  ring_coords <- function(r) {
    cc <- V[vi_o[starts[r]:ends[r]], , drop = FALSE]
    rbind(cc, cc[1, ])
  }
  geoms <- lapply(seq_len(m$n_mun), function(i) {
    rs <- by_mun[[as.character(i)]]
    if (is.null(rs)) return(st_multipolygon())
    ext <- rs[!m$ring_hole[rl$values[rs]]]
    hol <- rs[ m$ring_hole[rl$values[rs]]]
    st_multipolygon(lapply(seq_along(ext), function(k)
      if (k == 1L) c(list(ring_coords(ext[k])), lapply(hol, ring_coords))
      else list(ring_coords(ext[k]))))
  })
  out <- st_sf(region = mun$region, state_code = mun$state_code,
               geometry = st_sfc(geoms, crs = 6372))
  results[[s]] <- st_transform(st_make_valid(out), 4326)
}

hex_v9 <- do.call(rbind, results[!vapply(results, is.null, TRUE)])
hex_v9 <- hex_v9[order(hex_v9$region), ]

if (is.null(ONLY)) {
  saveRDS(hex_v9, "data/mxmunicipio_equalarea_v9.rds")
  cat(sprintf("Saved data/mxmunicipio_equalarea_v9.rds (%d cells)\n",
              nrow(hex_v9)))
} else {
  saveRDS(hex_v9, "/tmp/v9_states.rds")
  cat(sprintf("Saved /tmp/v9_states.rds (%d cells)\n", nrow(hex_v9)))
}
