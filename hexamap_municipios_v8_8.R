suppressPackageStartupMessages({
  library(sf); library(dplyr); library(clue); library(mxmaps)
})
sf_use_s2(FALSE)

# v8.8: position anchor floored so it applies to DENSE states too (79% of
# municipios previously had anchor == 0 exactly).
# v8.6: local_search now optimises the MEASURED objective -- true cell
# adjacency (cadj) instead of a site-distance proxy -- with a correct swap
# delta and the cell/municipio index confusion fixed. Everything else is
# identical to v8.5 so the comparison is clean.
# v8.5: Edomex individual tuning (stronger pull, more rounds).
# v8.4: alternating search <-> equalization rounds to close the dense-state
# cohesion gap (Oaxaca/Veracruz/Puebla/Michoacan ~54-62% in v8.3).
#
# Rounds: [short CVT touch-up (geometry adapts)] -> [owner-swap local search
# (cohesion recovers)] x3, finishing on a search pass. Anchors, gamma
# weighting, coverage repair, sliver cleanup all carried from v8.3.
#
# Usage: Rscript hexamap_municipios_v8_4.R [state_code ...]
# Output: data/mxmunicipio_hex_sf_v8_8.rds (full) or /tmp/v88_states.rds.

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
ANCHOR_FLOOR <- { v <- Sys.getenv("V88_ANCHOR"); if (nzchar(v)) as.numeric(v) else 0.30 }

args <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
ONLY <- if (length(args) >= 1 && !any(is.na(args))) sprintf("%02d", args) else NULL

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
    if (!all(coords[1,] == coords[nrow(coords),])) coords <- rbind(coords, coords[1,])
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

mxmun_km2 <- mxmun_sf |>
  st_transform(6372) |>
  mutate(true_km2 = as.numeric(st_area(geometry)) / 1e6) |>
  st_drop_geometry() |>
  select(region, state_code, true_km2)

state_areas <- mxmun_km2 |>
  group_by(state_code) |>
  summarise(state_km2 = sum(true_km2), n = dplyr::n(), .groups = "drop") |>
  mutate(cs_km = sqrt(state_km2 / n / 0.866))
cs_national <- sqrt(sum(state_areas$state_km2) / sum(state_areas$n) / 0.866)
state_areas$gamma_auto <- pmax(pmin(0.18 * (state_areas$cs_km / cs_national - 0.75),
                                    0.45), 0)
gamma_override <- c("02" = 0.35, "03" = 0.45, "05" = 0.30, "08" = 0.30,
                    "26" = 0.25, "10" = 0.25, "07" = 0.12, "04" = 0.15,
                    "24" = 0.12, "01" = 0.10, "09" = 0.10, "06" = 0.12)
state_areas$gamma <- state_areas$gamma_auto
state_areas$gamma[match(names(gamma_override), state_areas$state_code)] <-
  gamma_override

state_sf <- mxmun_sf |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_make_valid()

adj_all <- st_intersects(mxmun_sf, mxmun_sf)
adj_pairs <- list()
for (i in seq_along(regions)) {
  js <- adj_all[[i]]
  js <- js[js != i & mxmun_sf$state_code[js] == mxmun_sf$state_code[i]]
  adj_pairs[[i]] <- js
}

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
  pieces <- st_cast(g, "POLYGON")
  ar <- vapply(seq_along(pieces), function(j) poly_area(pieces[j]), numeric(1))
  pieces[[which.max(ar)]]
}

clip_hp <- function(pts, a, b) {
  n <- nrow(pts); out <- vector("list", n + 1); m <- 0
  for (i in seq_len(n)) {
    cur <- pts[i, ]; nxt <- pts[(i %% n) + 1, ]
    dc <- sum(a * cur); dn <- sum(a * nxt)
    if (dc <= b) {
      m <- m + 1; out[[m]] <- cur
      if (dn > b) { t <- (b - dc) / (dn - dc); m <- m + 1
                    out[[m]] <- cur + t * (nxt - cur) }
    } else if (dn <= b) {
      t <- (b - dc) / (dn - dc); m <- m + 1
      out[[m]] <- cur + t * (nxt - cur)
    }
  }
  if (m < 3) return(NULL)
  do.call(rbind, out[seq_len(m)])
}

local_search <- function(assigned, pos, st_pairs, mun_xy, lambda,
                         cadj_m = NULL, touch_d, max_pass = 25) {
  # v8.6 rewrite. Three defects in the v8.4/v8.5 version, all of them the same
  # thing -- the search was not optimising the quantity that gets measured:
  #
  #  1. PROXY OBJECTIVE. viol_of() tested whether two SITES were further apart
  #     than touch_d. The metric tests whether two CELLS share a boundary. In a
  #     power diagram neither implies the other: two sites can be close with a
  #     third cell wedged between them (no shared face), and two elongated
  #     cells can share a long face with distant centroids. cadj_m -- the true
  #     cell adjacency -- was already computed by the caller and left unused.
  #
  #  2. WRONG DELTA. Swapping municipios a and c changes the broken-adjacency
  #     count of every TRUE NEIGHBOUR of a and of c, not just a and c. The old
  #     dJ summed only over {a, c}, so it accepted and rejected on a number
  #     that was not the change in J.
  #
  #  3. INDEX CONFUSION. cand was built from which(dd < ...), i.e. CELL
  #     indices, then used as MUNICIPIO indices in disp_of(c) and in the swap.
  #     That is only correct while `assigned` is still the identity permutation
  #     (which it is on entry), so it silently swapped the wrong municipios
  #     from the first accepted move onward.
  n <- length(assigned)
  inv <- integer(n); inv[assigned] <- seq_len(n)      # cell -> municipio

  disp_of <- function(k, a) sqrt(sum((mun_xy[k, ] - pos[a[k], ])^2))
  # cadj_m = NULL is the pre-tessellation call site, where no cells exist yet
  # and only seed positions are available -- the distance proxy is all there is
  # at that stage, and that use of it is legitimate.
  viol_of <- if (is.null(cadj_m)) function(k, a) {
    js <- st_pairs[[k]]
    if (length(js) == 0) return(0L)
    sum(vapply(js, function(j)
      as.integer(sqrt(sum((pos[a[k], ] - pos[a[j], ])^2)) > touch_d), 0L))
  } else function(k, a) {
    js <- st_pairs[[k]]
    if (length(js) == 0) return(0L)
    sum(!cadj_m[a[k], a[js]])
  }
  Jk <- function(k, a) disp_of(k, a) + lambda * viol_of(k, a)

  improved <- TRUE; pass <- 0
  while (improved && pass < max_pass) {
    improved <- FALSE; pass <- pass + 1
    for (k in sample(seq_len(n))) {
      if (length(st_pairs[[k]]) == 0) next
      # candidates: municipios currently holding a cell near k's neighbours'
      # cells. touch_d is only a proximity heuristic for generating candidates
      # here -- it no longer enters the objective.
      near_cells <- unique(unlist(lapply(st_pairs[[k]], function(b) {
        nb <- assigned[b]
        dd <- sqrt((pos[, 1] - pos[nb, 1])^2 + (pos[, 2] - pos[nb, 2])^2)
        which(dd < 6 * touch_d)
      })))
      cand <- inv[near_cells[near_cells >= 1 & near_cells <= n]]
      cand <- unique(cand[cand > 0 & cand != k])
      if (length(cand) == 0) next

      best_dJ <- -1e-9; best_c <- NA_integer_
      for (c in cand) {
        aff <- unique(c(k, c, st_pairs[[k]], st_pairs[[c]]))
        a2 <- assigned; a2[c(k, c)] <- assigned[c(c, k)]
        dJ <- sum(vapply(aff, Jk, 0, a = a2)) -
              sum(vapply(aff, Jk, 0, a = assigned))
        if (!is.na(dJ) && dJ < best_dJ) { best_dJ <- dJ; best_c <- c }
      }
      if (!is.na(best_c)) {
        assigned[c(k, best_c)] <- assigned[c(best_c, k)]
        inv[assigned[c(k, best_c)]] <- c(k, best_c)
        improved <- TRUE
      }
    }
  }
  list(assigned = assigned,
       viol = sum(vapply(seq_len(n), viol_of, 0L, a = assigned)) %/% 2)
}

# CVT with explicit returns (no closure state)
cvt_run <- function(P, w, tgt, alpha, n_iter, damp, tol,
                    state_clip, sq, mun_xy, anchor, true_nb = NULL) {
  n <- nrow(P)
  cells <- vector("list", n); areas <- numeric(n)
  cents <- matrix(NA_real_, n, 2)
  for (it in seq_len(n_iter)) {
    d_mat <- as.matrix(dist(P))
    knn <- apply(d_mat, 1, function(r) order(r)[2:min(ncol(d_mat), 17)])
    cells <- vector("list", n); areas <- numeric(n)
    cents <- matrix(NA_real_, n, 2)
    for (i in seq_len(n)) {
      pts <- sq
      for (j in knn[, i]) {
        a <- 2 * (P[j, ] - P[i, ])
        b2 <- sum(P[j, ]^2) - sum(P[i, ]^2) + w[i] - w[j]
        pts <- clip_hp(pts, a, b2)
        if (is.null(pts)) break
      }
      if (is.null(pts)) { cells[[i]] <- NULL; areas[i] <- 0; next }
      poly <- st_polygon(list(rbind(pts, pts[1, ])))
      clipped <- suppressWarnings(st_intersection(st_sfc(poly, crs = 4326),
                                                   state_clip))
      if (length(clipped) == 0) { cells[[i]] <- NULL; areas[i] <- 0; next }
      areas[i] <- poly_area(clipped)
      cc <- st_coordinates(suppressWarnings(st_centroid(clipped)))
      cents[i, ] <- c(cc[1, 1], cc[1, 2])
      cells[[i]] <- clipped[[1]]
    }
    ok <- areas > 0
    w[ok] <- pmax(pmin(w[ok] + alpha * (tgt[ok] - areas[ok]), 4 * max(tgt)),
                  -4 * max(tgt))
    w <- w - mean(w)
    P[ok, ] <- P[ok, ] + damp * (cents[ok, ] - P[ok, ]) +
      anchor * (mun_xy[ok, ] - P[ok, ])
    if (!is.null(true_nb)) {
      for (i in which(ok)) {
        nb <- true_nb[[i]]
        if (length(nb) > 0) {
          pull <- colMeans(P[nb, , drop = FALSE])
          P[i, ] <- P[i, ] + pull_mu * (pull - P[i, ])
        }
      }
    }
    empty <- which(!ok)
    if (length(empty) > 0) {
      for (i in empty) {
        np <- st_nearest_points(st_sfc(st_point(P[i, ]), crs = 4326),
                                state_clip)
        cc2 <- st_coordinates(np)[2, ]
        P[i, ] <- P[i, ] + 0.5 * (c(cc2[1], cc2[2]) - P[i, ])
        w[i] <- w[i] + 0.3 * tgt[i]
      }
      w <- w - mean(w)
    }
    rel <- areas[ok] / tgt[ok]
    if (length(empty) == 0 && all(abs(rel - 1) < tol)) break
  }
  list(P = P, w = w, cells = cells, areas = areas, cents = cents)
}

coverage_repair <- function(cells, state_clip, st_area_deg2) {
  nonnull <- which(!vapply(cells, is.null, TRUE))
  cov <- st_union(st_sfc(cells[nonnull], crs = 4326))
  uncovered <- suppressWarnings(st_difference(state_clip, cov))
  if (inherits(uncovered, "sfc") && length(uncovered) > 0) {
    uncovered <- uncovered[[1]]
  }
  if (length(uncovered) > 0 && poly_area(uncovered) > 1e-10) {
    upieces <- st_cast(st_sfc(uncovered, crs = 4326), "POLYGON")
    cell_geoms <- st_sfc(cells[nonnull], crs = 4326)
    cell_cxy <- st_coordinates(st_centroid(cell_geoms))
    for (up in seq_len(length(upieces))) {
      pc <- upieces[up]
      pxy <- st_coordinates(st_centroid(pc))[1, ]
      dd <- sqrt((cell_cxy[, 1] - pxy[1])^2 + (cell_cxy[, 2] - pxy[2])^2)
      tgt_cell <- nonnull[which.min(dd)]
      cells[[tgt_cell]] <- st_union(st_sfc(cells[[tgt_cell]], crs = 4326),
                                    pc)[[1]]
    }
  }
  for (i in seq_along(cells)) {
    if (is.null(cells[[i]])) next
    g <- st_sfc(cells[[i]], crs = 4326)
    if (inherits(g[[1]], "MULTIPOLYGON")) cells[[i]] <- largest_piece(g)
  }
  cells
}

largest_piece <- function(g) {
  if (!inherits(g[[1]], "MULTIPOLYGON")) return(g[[1]])
  pieces <- st_cast(g, "POLYGON")
  ar <- vapply(seq_along(pieces), function(j) poly_area(pieces[j]), numeric(1))
  pieces[[which.max(ar)]]
}

viol_base <- function(assigned, seeds, st_pairs, cs) {
  sum(vapply(seq_len(length(assigned)), function(k) {
    if (length(st_pairs[[k]]) == 0) return(0L)
    sum(vapply(st_pairs[[k]], function(j)
      as.integer(sqrt(sum((seeds[assigned[k], ] - seeds[assigned[j], ])^2)) >
                   1.2 * cs), 0L))
  }, 0L)) %/% 2
}

state_codes <- sort(unique(mxmun_sf$state_code))
if (!is.null(ONLY)) state_codes <- intersect(state_codes, ONLY)

results <- vector("list", length(state_codes))
names(results) <- state_codes

for (s in state_codes) {
  state_poly <- state_sf$geometry[state_sf$state_code == s]
  state_muns <- mxmun_sf[mxmun_sf$state_code == s, ]
  n_target   <- nrow(state_muns)
  mun_idx    <- match(state_muns$region, mxmun_sf$region)
  st_area_deg2 <- poly_area(state_poly)
  cellsize <- sqrt(st_area_deg2 / n_target / 0.866)
  gamma_s <- state_areas$gamma[state_areas$state_code == s]
  cs_ratio <- state_areas$cs_km[state_areas$state_code == s] / cs_national
  # v8.8: give the position anchor a FLOOR. The v8.2-v8.6 form
  #   anchor = max(min(0.35*(cs_ratio - 1), 0.7), 0)
  # keys off cs_ratio - 1, so it is exactly 0 for every state whose cells are
  # smaller than the national average -- i.e. every DENSE state. That is 17 of
  # 32 states covering 1954 of 2469 municipios (79%): Oaxaca 0.46, Puebla 0.45,
  # Edomex 0.47, Hidalgo 0.56, Veracruz 0.65, Yucatan 0.68, Michoacan 0.81,
  # Chiapas 0.87, Jalisco 0.89, Guerrero 1.00 -- all unanchored. The anchor
  # arrived in v8.2 as a SPARSE-state fix (BC/BCS) and was never generalised.
  # Sparse states are unaffected (they sit above the floor).
  anchor <- pmax(pmin(0.35 * (cs_ratio - 1), 0.7), ANCHOR_FLOOR)
  # Per-state cohesion tuning: Edomex (hole + arms) needs a stronger
  # neighbor pull and more alternating rounds to converge
  pull_mu <- if (s == "15") 0.35 else 0.25
  n_rounds <- if (s == "15") 5 else 3

  ta <- mxmun_km2$true_km2[mun_idx]
  w_t <- ta^gamma_s
  share <- st_area_deg2 / n_target
  tgt <- w_t / sum(w_t) * st_area_deg2
  for (r in 1:3) {
    tgt <- pmax(pmin(tgt, 2.2 * share), 0.45 * share)
    tgt <- tgt / sum(tgt) * st_area_deg2
  }

  if (length(state_poly[[1]]) > 1 || inherits(state_poly[[1]], "MULTIPOLYGON")) {
    pieces <- st_cast(state_poly, "POLYGON")
    p_areas <- vapply(seq_along(pieces), function(j) poly_area(pieces[j]),
                      numeric(1))
    state_clip <- st_sfc(pieces[[which.max(p_areas)]], crs = 4326)
  } else {
    state_clip <- st_sfc(state_poly[[1]], crs = 4326)
  }

  seeds <- NULL; cs <- cellsize
  for (factor in c(1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60)) {
    cs <- cellsize * factor
    grid <- st_make_grid(state_poly, cellsize = cs, square = FALSE)
    cents0 <- st_centroid(grid, warn = FALSE)
    suppressWarnings({
      d_poly <- as.numeric(st_distance(st_transform(cents0, 6372),
                                       st_transform(state_poly, 6372)))
    })
    inside <- lengths(st_within(cents0, state_poly)) > 0
    ord <- order(-as.integer(inside), d_poly)
    if (length(grid) >= n_target) {
      pick <- ord[seq_len(n_target)]
      seeds <- st_coordinates(cents0)[pick, , drop = FALSE]
      break
    }
  }
  if (is.null(seeds)) stop(sprintf("State %s: no lattice fit", s))

  suppressWarnings(mun_xy <- st_coordinates(st_centroid(state_muns)))
  set.seed(sum(as.integer(charToRaw(s))) * 7)  # deterministic per state
  cost <- sqrt(outer(mun_xy[, 1], seeds[, 1], "-")^2 +
               outer(mun_xy[, 2], seeds[, 2], "-")^2)
  assigned <- as.integer(solve_LSAP(cost))

  st_pairs <- list()
  for (k in seq_len(n_target)) {
    gi <- mun_idx[k]
    js_reg <- mxmun_sf$region[adj_pairs[[gi]]]
    js_st <- match(js_reg, state_muns$region)
    st_pairs[[k]] <- js_st[!is.na(js_st)]
  }
  lambda <- 0.5 * cs
  if (n_target <= 9) {
    all_perms <- function(v) {
      if (length(v) == 1) return(matrix(v, nrow = 1))
      do.call(rbind, lapply(seq_along(v), function(i) cbind(v[i], all_perms(v[-i]))))
    }
    pm <- all_perms(seq_len(n_target))
    disp_sum <- rowSums(matrix(cost[cbind(rep(seq_len(n_target), each = nrow(pm)),
                                          as.vector(t(pm)))],
                               nrow = nrow(pm)))
    sd_mat <- as.matrix(dist(seeds))
    vp <- rep(0, nrow(pm))
    for (k in seq_len(n_target)) {
      for (j in st_pairs[[k]]) {
        if (j <= k) next
        vp <- vp + as.numeric(sd_mat[cbind(pm[, k], pm[, j])] > 1.2 * cs)
      }
    }
    assigned <- pm[which.min(disp_sum + lambda * vp), ]
  } else {
    ls1 <- local_search(assigned, seeds, st_pairs, mun_xy, lambda,
                      touch_d = 1.2 * cs)
    assigned <- ls1$assigned
  }
  P <- seeds[assigned, , drop = FALSE]
  rownames(P) <- NULL

  bb <- st_bbox(state_poly)
  pad <- 0.5 * max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"]) + cs
  sq <- matrix(c(bb["xmin"] - pad, bb["ymin"] - pad,
                 bb["xmax"] + pad, bb["ymin"] - pad,
                 bb["xmax"] + pad, bb["ymax"] + pad,
                 bb["xmin"] - pad, bb["ymax"] + pad), ncol = 2, byrow = TRUE)
  w <- rep(0, n_target)

  r1 <- cvt_run(P, w, tgt, 0.6, 60, 1.0, 0.05, state_clip, sq, mun_xy, anchor, st_pairs)
  r2 <- cvt_run(r1$P, r1$w, tgt, 0.05, 10, 0.4, 1.0, state_clip, sq, mun_xy, anchor, st_pairs)
  r3 <- cvt_run(r2$P, r2$w, tgt, 0.3, 20, 0.6, 0.05, state_clip, sq, mun_xy, anchor, st_pairs)
  P <- r3$P; w <- r3$w; cells <- r3$cells; cents <- r3$cents
  cells <- coverage_repair(cells, state_clip, st_area_deg2)

  g_cells <- st_transform(st_sfc(cells, crs = 4326), 6372)
  cadj <- st_is_within_distance(g_cells, g_cells, dist = 1, sparse = TRUE)
  cadj_m <- as.matrix(cadj); diag(cadj_m) <- TRUE
  cell_pos <- cents
  na_rows <- which(!complete.cases(cell_pos))
  if (length(na_rows) > 0) cell_pos[na_rows, ] <- P[na_rows, ]
  spacing <- sqrt(mean(tgt))

  assigned2 <- seq_len(n_target)
  n_restart <- if (s == "15") 3 else 1
  best_coh <- -1
  best <- list(assigned2 = assigned2, cells = cells, cell_pos = cell_pos)
  for (rs in seq_len(n_restart)) {
    # strtoi() defaults to base = 0L, which reads a leading zero as OCTAL, so
    # strtoi("08") and strtoi("09") are NA and set.seed(NA) errors. This is why
    # v8.5 crashes at Chihuahua and cannot reproduce its own committed
    # artifact. as.integer is base 10 and safe for the "01".."32" codes.
    set.seed(rs * 1000 + as.integer(s))
    if (rs > 1) {
      ru <- cvt_run(P, w, tgt, 0.25, 12, 0.5, 1.0, state_clip, sq, mun_xy,
                    anchor, st_pairs)
      P <- ru$P; w <- ru$w; cells <- ru$cells; cents <- ru$cents
      cells <- coverage_repair(cells, state_clip, st_area_deg2)
    }
    g_cells <- st_transform(st_sfc(cells, crs = 4326), 6372)
    cadj <- st_is_within_distance(g_cells, g_cells, dist = 1, sparse = TRUE)
    cadj_m <- matrix(FALSE, n_target, n_target)
    for (ci in seq_len(n_target)) cadj_m[ci, cadj[[ci]]] <- TRUE
    diag(cadj_m) <- TRUE
    cell_pos <- cents
    na_rows <- which(!complete.cases(cell_pos))
    if (length(na_rows) > 0) cell_pos[na_rows, ] <- P[na_rows, ]
    ls <- local_search(assigned2, cell_pos, st_pairs, mun_xy,
                       1.0 * spacing, cadj_m, 1.4 * spacing, max_pass = 15)
    assigned2 <- ls$assigned
    pres <- 0L; tot <- 0L
    for (k in seq_len(n_target)) {
      js <- st_pairs[[k]]
      if (length(js) == 0) next
      for (j in js) {
        if (j <= k) next
        tot <- tot + 1L
        ak <- assigned2[k]; aj <- assigned2[j]
        if (length(cadj[[ak]][cadj[[ak]] == aj]) > 0) pres <- pres + 1L
      }
    }
    coh <- pres / max(tot, 1)
    cat(sprintf("  state %s restart %d: cohesion %.1f%%\n", s, rs, 100 * coh))
    if (coh > best_coh) {
      best_coh <- coh
      best <- list(assigned2 = assigned2, cells = cells, cell_pos = cell_pos)
    }
  }
  assigned2 <- best$assigned2; cells <- best$cells; cell_pos <- best$cell_pos

  output_geoms <- st_sfc(lapply(seq_len(n_target), function(i) {
    ci <- assigned2[i]
    if (!is.null(cells[[ci]])) return(largest_piece(st_sfc(cells[[ci]], crs = 4326)))
    st_buffer(st_sfc(st_point(cell_pos[ci, ]), crs = 4326),
              0.05 * sqrt(tgt[ci]))[[1]]
  }), crs = 4326)

  results[[s]] <- st_sf(region = state_muns$region, state_code = s,
                        geometry = output_geoms, crs = 4326)

  fa <- as.numeric(st_area(st_transform(results[[s]], 6372)))
  ratio_tg <- fa / (sum(fa) * tgt / sum(tgt))
  pres <- 0L; tot <- 0L
  g6 <- st_transform(results[[s]], 6372)
  fin <- st_is_within_distance(g6, g6, dist = 1, sparse = TRUE)
  for (k in seq_len(n_target)) {
    js <- st_pairs[[k]]
    if (length(js) == 0) next
    for (j in js) {
      if (j <= k) next
      tot <- tot + 1L
      if (length(fin[[k]][fin[[k]] == j]) > 0) pres <- pres + 1L
    }
  }
  disp <- sqrt((cell_pos[assigned2, 1] - mun_xy[, 1])^2 +
               (cell_pos[assigned2, 2] - mun_xy[, 2])^2)
  cat(sprintf(
    "State %s: area-vs-target p5=%.2f p95=%.2f | adj %.0f%%->%.0f%% | disp=%.3f\n",
    s, quantile(ratio_tg, .05), quantile(ratio_tg, .95),
    100 * viol_base(assigned, seeds, st_pairs, cs) / max(tot, 1),
    100 * pres / max(tot, 1), mean(disp)))
}

hex_v8 <- do.call(rbind, results)
stopifnot(nrow(hex_v8) == sum(mxmun_sf$state_code %in% state_codes))
stopifnot(inherits(hex_v8$geometry, "sfc"))

if (is.null(ONLY)) {
  stopifnot(nrow(hex_v8) == nrow(mxmun_sf))
  saveRDS(hex_v8, "data/mxmunicipio_hex_sf_v8_8.rds")
  cat(sprintf("Saved data/mxmunicipio_hex_sf_v8_8.rds (%d cells)\n", nrow(hex_v8)))
} else {
  saveRDS(hex_v8, "/tmp/v88_states.rds")
  cat(sprintf("Saved /tmp/v88_states.rds (%d cells)\n", nrow(hex_v8)))
}
