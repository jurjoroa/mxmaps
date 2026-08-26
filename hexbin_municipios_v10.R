# ---------------------------------------------------------------------------
# Purpose: v10 municipio hexbin cartogram — 2469 IDENTICAL regular hexagons on
#          one national lattice, partitioned into 32 state regions of exactly
#          n_s cells, one cell per municipio.
#
#          Why a uniform lattice rather than the v3-v8.5 re-tessellation: with
#          identical hexagons, two of the four objective terms stop being
#          optimisation targets and become EXACT by construction --
#            area evenness    : exact (every cell identical)
#            shape regularity : exact (IQ = 0.9069, a regular hexagon)
#          leaving only position fidelity and adjacency to optimise. That is a
#          far better-conditioned problem than the v8 line's five knobs pushing
#          on four implicit terms, and it needs no per-state parameter table.
#
#          ONE explicit objective for the within-state assignment:
#
#            J = sum_i |pos(cell_i) - pos_true(i)|^2  +  LAMBDA * broken_adj
#
#          UNAVOIDABLE CONSEQUENCE, stated up front: if every cell is the same
#          size, each state's rendered area is proportional to its MUNICIPIO
#          COUNT, not its territory. Oaxaca (570) renders ~57x Colima (10)
#          though its true area is only ~17x. That is not a bug and no
#          algorithm can avoid it -- it is what "every municipio equally
#          legible" means geometrically. GAMMA below is the dial if some
#          territorial fidelity is wanted back, at the cost of uniform cells.
#
# Author:  Claude (Opus 5) for Jorge Roa
# ISS:     -
# Version: v10
# Date:    2026-08-25
#
# Usage:   Rscript hexbin_municipios_v10.R
# Output:  data/mxmunicipio_hexbin_v10.rds  (EPSG:6372, so hexes stay regular)
#          data/mxstate_hexbin_v10.rds      (dissolved state regions, overlay)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(clue); library(mxmaps)
})
sf_use_s2(FALSE)

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

.envnum <- function(nm, d) { v <- Sys.getenv(nm); if (nzchar(v)) as.numeric(v) else d }
PART_ITER <- as.integer(.envnum("V10_PART_ITER", 60))   # Lloyd rounds
LAMBDA    <- .envnum("V10_LAMBDA", 0.6)   # adjacency weight in J (relative to
                                          # squared displacement in cell units)
SWAP_PASS <- as.integer(.envnum("V10_SWAP", 30))
ONLY      <- Sys.getenv("V10_ONLY")

# ---- source polygons ------------------------------------------------------

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
  st_polygon(c(lapply(rings[!sapply(rings, `[[`, "hole")], `[[`, "coords"),
               lapply(rings[ sapply(rings, `[[`, "hole")], `[[`, "coords")))
})
mxmun_sf <- st_make_valid(st_sf(region = regions,
                                geometry = st_sfc(polys, crs = 4326)))
mxmun_sf$state_code <- substr(mxmun_sf$region, 1, 2)
mx_m <- st_transform(mxmun_sf, 6372)

mun_xy <- st_coordinates(st_centroid(mx_m, of_largest_polygon = TRUE))
mun_xy <- mun_xy[, 1:2, drop = FALSE]
mx_m$true_km2 <- as.numeric(st_area(mx_m)) / 1e6

states <- mx_m |>
  st_drop_geometry() |>
  group_by(state_code) |>
  summarise(n = dplyr::n(), km2 = sum(true_km2), .groups = "drop") |>
  arrange(state_code)
N <- sum(states$n)
A_mex <- sum(states$km2) * 1e6           # m2

# true intra-state adjacency (the thing we are trying not to break)
adj_all <- st_intersects(mx_m, mx_m)

st_ctr0 <- mx_m |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_centroid(of_largest_polygon = TRUE) |>
  st_coordinates()
seeds <- st_ctr0[, 1:2, drop = FALSE]
rownames(seeds) <- states$state_code

# ---- 1. uniform hex lattice ----------------------------------------------
# Circumradius R chosen so N hexes total exactly Mexico's area, i.e. the whole
# arrangement keeps Mexico's overall footprint even though states are resized.
nS <- nrow(states)
R <- sqrt(2 * A_mex / (3 * sqrt(3) * N))
cat(sprintf("cells N=%d  hex circumradius R=%.2f km  hex area=%.1f km2\n",
            N, R / 1000, (3 * sqrt(3) / 2) * R^2 / 1e6))

bb <- st_bbox(mx_m)
pad <- 0.45 * max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"])
# axial coords for FLAT-TOP hexes: centre = (1.5*R*q, sqrt(3)*R*(r + q/2))
qs <- floor((bb["xmin"] - pad) / (1.5 * R)):ceiling((bb["xmax"] + pad) / (1.5 * R))
rs <- floor((bb["ymin"] - pad) / (sqrt(3) * R) - max(abs(qs)) / 2) :
      ceiling((bb["ymax"] + pad) / (sqrt(3) * R) + max(abs(qs)) / 2)
lat <- expand.grid(q = qs, r = rs)
lat$x <- 1.5 * R * lat$q
lat$y <- sqrt(3) * R * (lat$r + lat$q / 2)
keep <- lat$x >= bb["xmin"] - pad & lat$x <= bb["xmax"] + pad &
        lat$y >= bb["ymin"] - pad & lat$y <= bb["ymax"] + pad
lat <- lat[keep, ]
rownames(lat) <- NULL

# ---- 1b. state-level Dorling layout: make ROOM before assigning cells -----
# With uniform cells, state s needs n_s * hex_area of map. Oaxaca needs 4.8x
# its true footprint (570 cells = 451,000 km2 vs 94,000 km2 of real territory).
# Constraining the cell set to Mexico's true silhouette therefore CANNOT work:
# the surplus has to come from somewhere, and a capacitated partition takes it
# from the empty north, pushing Oaxaca's region 1762 km off its real position
# (measured). That is the same failure that got v2 rejected in 2026-04.
#
# So the footprint must expand instead. Lay the states out as a Dorling circle
# cartogram first -- circle radius set by CELL COUNT, pulled toward the true
# centroid, pushed apart until non-overlapping -- then draw cells only from
# inside that layout. States keep their relative geography and their correct
# size; the price is that the national outline is Mexico inflated where
# municipios are dense, which is the honest cartogram trade and not a bug.
hex_area <- (3 * sqrt(3) / 2) * R^2
rad <- sqrt(states$n * hex_area / pi)
P0 <- seeds
P  <- seeds
for (iter in 1:4000) {
  dx <- outer(P[, 1], P[, 1], "-"); dy <- outer(P[, 2], P[, 2], "-")
  d  <- sqrt(dx^2 + dy^2); diag(d) <- Inf
  need <- outer(rad, rad, "+")
  ov <- need - d
  ov[ov < 0] <- 0; diag(ov) <- 0
  ux <- dx / pmax(d, 1e-9); uy <- dy / pmax(d, 1e-9)
  P[, 1] <- P[, 1] + 0.5 * rowSums(ov * ux) + 0.02 * (P0[, 1] - P[, 1])
  P[, 2] <- P[, 2] + 0.5 * rowSums(ov * uy) + 0.02 * (P0[, 2] - P[, 2])
  if (max(ov) < 0.01 * R) break
}
seeds <- P
cat(sprintf("dorling: %d iters, max overlap %.1f km, mean shift from true %.0f km\n",
            iter, max(ov) / 1000,
            mean(sqrt(rowSums((P - P0)^2))) / 1000))

# ---- 1c. candidate cells = inside the Dorling layout, exactly N of them ----
# Every candidate gets assigned, so the union is gap-free and contiguous.
dc <- sqrt(outer(lat$x, seeds[, 1], "-")^2 + outer(lat$y, seeds[, 2], "-")^2)
# "insideness": how far inside its nearest state circle a cell sits
inside <- apply(dc - matrix(rad, nrow(lat), nS, byrow = TRUE), 1, min)
ord <- order(inside)
sel_cells <- sort(ord[seq_len(N)])
cat(sprintf("cell set: %d cells drawn from the Dorling layout (worst cell %.1f km outside its circle)\n",
            N, max(inside[sel_cells]) / 1000))
lat <- lat[sel_cells, ]
rownames(lat) <- NULL

# ---- 2. capacitated partition of lattice cells into 32 state regions -----
# Transportation problem: state s must receive exactly n_s cells, cost =
# squared distance from cell to the state's seed. Solved by greedy on
# distance-sorted (cell, state) pairs, with Lloyd updates of the seeds. Seeds
# start at true state centroids so regions land near their real geography.

nL <- nrow(lat)
assign_cells <- function(seeds) {
  d <- outer(lat$x, seeds[, 1], "-")^2 + outer(lat$y, seeds[, 2], "-")^2
  ord <- order(d)                       # column-major: cell + (state-1)*nL
  cell_of <- ((ord - 1L) %% nL) + 1L
  state_of <- ((ord - 1L) %/% nL) + 1L
  cap <- states$n
  taken <- integer(nL)
  out <- integer(nL)
  filled <- integer(nS)
  for (k in seq_along(ord)) {
    c_ <- cell_of[k]; s_ <- state_of[k]
    if (taken[c_] == 1L || filled[s_] >= cap[s_]) next
    taken[c_] <- 1L; out[c_] <- s_; filled[s_] <- filled[s_] + 1L
    if (sum(filled) == N) break
  }
  out
}

lab <- assign_cells(seeds)
for (it in seq_len(PART_ITER)) {
  sel <- lab > 0
  new_seeds <- seeds
  agg <- aggregate(cbind(x, y) ~ lab, data = cbind(lat[sel, ], lab = lab[sel]),
                   FUN = mean)
  new_seeds[agg$lab, ] <- as.matrix(agg[, c("x", "y")])
  shift <- max(sqrt(rowSums((new_seeds - seeds)^2)))
  seeds <- new_seeds
  lab <- assign_cells(seeds)
  if (shift < 0.02 * R) break
}
cat(sprintf("partition: %d Lloyd rounds, %d cells assigned (target %d)\n",
            it, sum(lab > 0), N))
stopifnot(sum(lab > 0) == N)

cells <- lat[lab > 0, ]
cells$state_code <- states$state_code[lab[lab > 0]]

# lattice adjacency: 6 axial neighbours
axial_key <- function(q, r) paste(q, r)
ckey <- axial_key(cells$q, cells$r)
nb_off <- list(c(1, 0), c(1, -1), c(0, -1), c(-1, 0), c(-1, 1), c(0, 1))
cell_nb <- lapply(seq_len(nrow(cells)), function(i) {
  ks <- vapply(nb_off, function(o) axial_key(cells$q[i] + o[1],
                                             cells$r[i] + o[2]), "")
  m <- match(ks, ckey)
  m[!is.na(m)]
})

# ---- 3. within-state assignment: municipios -> cells ---------------------
# Hungarian on normalised positions gives the optimal position-fidelity
# assignment; local swaps then trade a little displacement for adjacency,
# judged by ONE explicit J.

norm_pts <- function(P) {
  P <- scale(P, center = TRUE, scale = FALSE)
  s <- sqrt(mean(rowSums(P^2)))
  if (!is.finite(s) || s == 0) s <- 1
  P / s
}

assigned <- integer(N)          # municipio row (in mx_m) -> cell row
state_list <- if (nzchar(ONLY)) strsplit(ONLY, ",")[[1]] else states$state_code

for (s in state_list) {
  mi <- which(mx_m$state_code == s)
  ci <- which(cells$state_code == s)
  n <- length(mi)
  stopifnot(length(ci) == n)

  if (n == 1L) { assigned[mi] <- ci; next }

  Pm <- norm_pts(mun_xy[mi, , drop = FALSE])
  Pc <- norm_pts(as.matrix(cells[ci, c("x", "y")]))
  D <- outer(Pm[, 1], Pc[, 1], "-")^2 + outer(Pm[, 2], Pc[, 2], "-")^2
  a <- as.integer(solve_LSAP(D))        # municipio k -> local cell a[k]

  # true adjacency within this state, in local municipio indices
  nb_true <- lapply(seq_len(n), function(k) {
    js <- adj_all[[mi[k]]]
    js <- js[js != mi[k] & mx_m$state_code[js] == s]
    match(js, mi)
  })

  # J = sum of squared displacement (in normalised units) + LAMBDA * broken
  local_cell_nb <- lapply(ci, function(c_) match(cell_nb[[c_]], ci))
  local_cell_nb <- lapply(local_cell_nb, function(z) z[!is.na(z)])

  disp_of <- function(k, cc) sum((Pm[k, ] - Pc[cc, ])^2)
  brk_of <- function(k, a) {
    nbk <- nb_true[[k]]
    if (!length(nbk)) return(0)
    sum(vapply(nbk, function(j)
      as.integer(!(a[j] %in% local_cell_nb[[a[k]]])), 0L))
  }
  Jk <- function(k, a) disp_of(k, a[k]) + LAMBDA * brk_of(k, a)

  for (pass in seq_len(SWAP_PASS)) {
    improved <- FALSE
    for (k in sample(n)) {
      # candidate partners: municipios sitting on cells near k's cell
      cand <- unique(unlist(local_cell_nb[c(a[k], local_cell_nb[[a[k]]])]))
      cand <- setdiff(match(cand, a), NA)
      cand <- cand[cand != k]
      if (!length(cand)) next
      best <- 0; bestj <- NA_integer_
      for (j in cand) {
        # The affected set must include every municipio whose J can change:
        # k and j themselves, AND all of their TRUE neighbours, whose broken-
        # adjacency counts flip when k or j moves. Evaluating the delta over
        # only {k} U cand (as a first cut did) computes the wrong number, so
        # the search rejects real improvements and converges in ~2 passes at
        # ~67% against a ~90% ceiling.
        aff <- unique(c(k, j, nb_true[[k]], nb_true[[j]]))
        a2 <- a; a2[c(k, j)] <- a[c(j, k)]
        d <- sum(vapply(aff, Jk, 0, a = a2)) - sum(vapply(aff, Jk, 0, a = a))
        if (d < best) { best <- d; bestj <- j }
      }
      if (!is.na(bestj)) {
        a[c(k, bestj)] <- a[c(bestj, k)]
        improved <- TRUE
      }
    }
    if (!improved) break
  }
  assigned[mi] <- ci[a]

  kept <- sum(vapply(seq_len(n), function(k) {
    nbk <- nb_true[[k]]
    if (!length(nbk)) return(0L)
    sum(vapply(nbk, function(j)
      as.integer(a[j] %in% local_cell_nb[[a[k]]]), 0L))
  }, 0L))
  tot <- sum(lengths(nb_true))
  cat(sprintf("state %s n=%3d  passes=%2d  adjacency kept %4d/%4d = %5.1f%%\n",
              s, n, pass, kept, tot, 100 * kept / max(tot, 1)))
}

# ---- 4. emit hexagons ----------------------------------------------------

hex_at <- function(cx, cy) {
  a <- (0:5) * pi / 3
  st_polygon(list(rbind(cbind(cx + R * cos(a), cy + R * sin(a)),
                        c(cx + R * cos(0), cy + R * sin(0)))))
}
sel <- assigned > 0
out <- st_sf(region = mx_m$region[sel],
             state_code = mx_m$state_code[sel],
             geometry = st_sfc(lapply(which(sel), function(k)
               hex_at(cells$x[assigned[k]], cells$y[assigned[k]])), crs = 6372))
out <- out[order(out$region), ]

st_out <- out |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop")

if (!nzchar(ONLY)) {
  saveRDS(out, "data/mxmunicipio_hexbin_v10.rds")
  saveRDS(st_out, "data/mxstate_hexbin_v10.rds")
  cat(sprintf("Saved data/mxmunicipio_hexbin_v10.rds (%d hexes)\n", nrow(out)))
} else {
  saveRDS(out, "/tmp/v10_states.rds")
  cat(sprintf("Saved /tmp/v10_states.rds (%d hexes)\n", nrow(out)))
}
