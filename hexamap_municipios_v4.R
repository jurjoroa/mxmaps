suppressPackageStartupMessages({
  library(sf); library(dplyr); library(clue); library(mxmaps)
})
sf_use_s2(FALSE)

# v4: pure-hex patch layout — "even space" refinement of v3.1.
#
# New trade-off from Jorge: cells no longer need to match the true state
# border or exact locations. Must hold: pure unclipped hexagons, even cell
# sizes, neighbors stay adjacent, arrangement reads as Mexico.
#
# Per state:
#   1. Largest-fit sizing: search factors 1.30 -> 0.60, take the largest
#      whose grid fits >= n_municipios cells inside the state polygon.
#      Cells are as large as the state can hold (more even than v3.1).
#   2. Patch selection: exactly n cells, ranked inside-first then by
#      distance to the state polygon. No clipping, no absorption, no blobs.
#      Shortage -> patch bulges past the true border; surplus -> patch
#      recedes from the border.
#   3. Cross-state overlap repair: contested cells go to the state whose
#      polygon is nearer; loser swaps in its nearest non-clashing candidate.
#   4. Connectivity repair: isolated patch cells swapped for connectors.
#   5. Assignment: Hungarian (clue::solve_LSAP) on centroid distances.
#
# Output: data/mxmunicipio_hex_sf_v4.rds — same schema as v3/v3.1.

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

# --- 1. Reconstruct municipio sf from mxmunicipio.map -----------------------
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
mxmun_sf <- st_sf(region = regions, geometry = st_sfc(polys, crs = 4326))
mxmun_sf$state_code <- substr(mxmun_sf$region, 1, 2)
mxmun_sf <- st_make_valid(mxmun_sf)

state_sf <- mxmun_sf |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_make_valid()

poly_area <- function(g) {
  if (inherits(g, "sfc")) g <- g[[1]]
  if (length(g) == 0) return(0)
  if (inherits(g, "MULTIPOLYGON")) {
    return(sum(vapply(g, function(p) poly_area(p[[1]]), numeric(1))))
  }
  if (inherits(g, "POLYGON")) g <- g[[1]]
  if (is.null(dim(g))) return(0)
  x <- g[, 1]; y <- g[, 2]
  abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
}

state_codes <- sort(unique(mxmun_sf$state_code))

# --- 2. Per-state patch selection -------------------------------------------
build_patch <- function(s, f_start = 1.30) {
  state_poly <- state_sf$geometry[state_sf$state_code == s]
  n_target   <- sum(mxmun_sf$state_code == s)
  cellsize <- sqrt(poly_area(state_poly) / n_target / 0.866)

  best <- NULL
  for (f in seq(f_start, 0.60, by = -0.05)) {
    cs <- cellsize * f
    grid <- st_make_grid(st_buffer(state_poly, 3 * cs),
                         cellsize = cs, square = FALSE)
    cents <- st_centroid(grid, warn = FALSE)
    inside <- lengths(st_within(cents, state_poly)) > 0
    if (sum(inside) >= n_target) {
      best <- list(cs = cs, grid = grid, inside = inside,
                   xy = st_coordinates(cents))
      break
    }
  }
  if (is.null(best)) {
    stop(sprintf("State %s: no factor fit %d cells inside", s, n_target))
  }

  # Rank: fully-contained hexes first, then centroid-inside, then by distance
  # to the state polygon. Fully-contained cells never clash across borders.
  full <- lengths(st_within(best$grid, state_poly)) > 0
  cents <- st_centroid(best$grid, warn = FALSE)
  suppressWarnings({
    d_poly <- as.numeric(st_distance(st_transform(cents, 6372),
                                     st_transform(state_poly, 6372)))
  })
  ord <- order(-as.integer(full), -as.integer(best$inside), d_poly)
  list(
    pool = list(cs = best$cs, grid = best$grid, xy = best$xy,
                inside = best$inside, full = full,
                d_poly = d_poly, rank = ord),
    sel  = ord[seq_len(n_target)]
  )
}

patches <- setNames(vector("list", length(state_codes)), state_codes)
for (s in state_codes) {
  patches[[s]] <- build_patch(s)
  p <- patches[[s]]
  n_bulge <- sum(!p$pool$inside[p$sel])
  cat(sprintf("State %s: %d cells, cs=%.3f (%d not inside, %d poking)\n",
              s, length(p$sel), p$pool$cs, n_bulge,
              sum(p$pool$inside[p$sel] & !p$pool$full[p$sel])))
}

# --- 3. Cross-state overlap repair ------------------------------------------
# Batch repair: each sweep finds all contested cross-state pairs and resolves
# them (loser = cell farther from its own state polygon; loser swaps in its
# nearest free candidate). Per-state unions cached per sweep, invalidated on
# change. Stuck states regenerate at a smaller cellsize for more slack.
regen <- setNames(integer(length(state_codes)), state_codes)
state_union <- function(s) st_union(patches[[s]]$pool$grid[patches[[s]]$sel])
others_u_of <- function(exclude) {
  st_union(do.call(c, lapply(state_codes[state_codes != exclude],
                             function(z) {
                               if (is.null(u_cache[[z]])) {
                                 u_cache[[z]] <<- state_union(z)
                               }
                               u_cache[[z]]
                             })))
}

for (sweep in seq_len(500)) {
  flat <- do.call(rbind, lapply(state_codes, function(s) {
    data.frame(state = s, grid_idx = patches[[s]]$sel)
  }))
  all_g <- do.call(c, lapply(seq_len(nrow(flat)), function(i) {
    patches[[flat$state[i]]]$pool$grid[flat$grid_idx[i]]
  }))
  inter <- st_intersects(all_g, all_g)
  pairs <- NULL
  for (i in seq_along(inter)) {
    js <- inter[[i]][inter[[i]] > i]
    js <- js[flat$state[js] != flat$state[i]]
    if (length(js) > 0) pairs <- rbind(pairs, cbind(i, js))
  }
  if (is.null(pairs)) {
    cat(sprintf("Cross-state overlaps resolved after %d sweeps\n", sweep - 1))
    break
  }
  if (sweep %% 25 == 0) {
    cat(sprintf("  sweep %d: %d contested pairs\n", sweep, nrow(pairs)))
  }
  u_cache <- setNames(vector("list", length(state_codes)), state_codes)
  changed <- character(0)
  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k, 1]; j <- pairs[k, 2]
    si <- flat$state[i]; sj <- flat$state[j]
    # recheck live: patches may have changed earlier this sweep
    gi <- patches[[si]]$pool$grid[patches[[si]]$sel]
    gj <- patches[[sj]]$pool$grid[patches[[sj]]$sel]
    a <- suppressWarnings(st_intersection(gi, gj))
    if (length(a) == 0 || poly_area(a) < 1e-8) next
    di <- patches[[si]]$pool$d_poly[flat$grid_idx[i]]
    dj <- patches[[sj]]$pool$d_poly[flat$grid_idx[j]]
    loser <- if (di <= dj) sj else si
    p <- patches[[loser]]
    # drop the cell of the loser that still intersects the other cell
    other_g <- if (loser == si) gj else gi
    sel_g <- p$pool$grid[p$sel]
    drop_pos <- which(vapply(seq_along(sel_g), function(m) {
      length(suppressWarnings(st_intersection(sel_g[m], other_g))) > 0
    }, logical(1)))
    if (length(drop_pos) == 0) next
    drop_pos <- drop_pos[1]
    others_u <- others_u_of(loser)
    cand <- p$pool$rank[!p$pool$rank %in% p$sel]
    placed <- FALSE
    for (cand_i in cand) {
      g <- p$pool$grid[cand_i]
      a <- suppressWarnings(st_intersection(g, others_u))
      if (length(a) == 0 || poly_area(a) < 1e-8) {
        p$sel[drop_pos] <- cand_i
        patches[[loser]] <- p
        u_cache[[loser]] <- NULL
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      regen[[loser]] <- regen[[loser]] + 1
      if (regen[[loser]] > 8) {
        stop(sprintf("State %s: overlap repair exhausted regenerations", loser))
      }
      f_new <- 1.30 * 0.90^regen[[loser]]
      patches[[loser]] <- build_patch(loser, f_start = f_new)
      u_cache[[loser]] <- NULL
      cat(sprintf("  regenerated state %s at f=%.2f (was stuck)\n",
                  loser, f_new))
    }
  }
}

# --- 4. Connectivity repair --------------------------------------------------
for (s in state_codes) {
  p <- patches[[s]]
  cs <- p$pool$cs
  # Others' cells can't change during this state's repair: cache the union.
  others_u <- st_union(do.call(c, lapply(
    state_codes[state_codes != s],
    function(z) patches[[z]]$pool$grid[patches[[z]]$sel])))
  n_repair <- 0
  repeat {
    xy <- p$pool$xy[p$sel, , drop = FALSE]
    n <- nrow(xy)
    d_mat <- as.matrix(dist(xy))
    adj_m <- d_mat <= 1.25 * cs
    diag(adj_m) <- FALSE
    seen <- logical(n); seen[1] <- TRUE; queue <- 1
    while (length(queue) > 0) {
      q <- queue[1]; queue <- queue[-1]
      nb <- which(adj_m[q, ] & !seen)
      seen[nb] <- TRUE; queue <- c(queue, nb)
    }
    if (all(seen)) break
    if (n_repair >= 30) {
      cat(sprintf("State %s: connectivity repair hit cap; leaving as is\n", s))
      break
    }
    n_repair <- n_repair + 1
    iso <- which(!seen)
    far <- iso[which.min(apply(d_mat[iso, -iso, drop = FALSE], 1, min))]
    kept <- p$sel[-far]
    kept_xy <- xy[-far, , drop = FALSE]
    cand <- p$pool$rank[!p$pool$rank %in% p$sel]
    placed <- FALSE
    for (cand_i in cand) {
      g <- p$pool$grid[cand_i]
      a <- suppressWarnings(st_intersection(g, others_u))
      if (length(a) > 0 && poly_area(a) > 1e-8) next
      cxy <- p$pool$xy[cand_i, ]
      if (min(sqrt((kept_xy[, 1] - cxy[1])^2 +
                   (kept_xy[, 2] - cxy[2])^2)) >= 1.25 * cs) {
        next
      }
      p$sel[far] <- cand_i
      patches[[s]] <- p
      placed <- TRUE
      cat(sprintf("  state %s: connectivity swap %d\n", s, n_repair))
      break
    }
    if (!placed) {
      cat(sprintf("State %s: connectivity repair found no swap; leaving as is\n", s))
      break
    }
  }
}

# --- 5. Assignment (Hungarian) + output --------------------------------------
results <- vector("list", length(state_codes))
names(results) <- state_codes

for (s in state_codes) {
  p <- patches[[s]]
  state_muns <- mxmun_sf[mxmun_sf$state_code == s, ]
  n_target <- nrow(state_muns)
  cell_geoms <- p$pool$grid[p$sel]
  suppressWarnings({
    mun_xy  <- st_coordinates(st_centroid(state_muns))
    cell_xy <- st_coordinates(st_centroid(cell_geoms))
  })
  cost <- sqrt(outer(mun_xy[, 1], cell_xy[, 1], "-")^2 +
               outer(mun_xy[, 2], cell_xy[, 2], "-")^2)
  assigned <- as.integer(solve_LSAP(cost))
  results[[s]] <- st_sf(
    region     = state_muns$region,
    state_code = s,
    geometry   = st_sfc(cell_geoms[assigned], crs = 4326),
    crs        = 4326
  )
  disp <- sqrt((cell_xy[assigned, 1] - mun_xy[, 1])^2 +
               (cell_xy[assigned, 2] - mun_xy[, 2])^2)
  cat(sprintf("State %s: assigned, disp mean=%.4f max=%.4f deg\n",
              s, mean(disp), max(disp)))
}

hex_v4 <- do.call(rbind, results)
stopifnot(nrow(hex_v4) == nrow(mxmun_sf))
stopifnot(setequal(hex_v4$region, mxmun_sf$region))

# Same-state overlap check (must be zero: one lattice per state)
same_ov <- 0
suppressWarnings({ inter_idx <- st_intersects(hex_v4, hex_v4) })
for (i in seq_along(inter_idx)) {
  js <- inter_idx[[i]]; js <- js[js > i]
  for (j in js) {
    if (hex_v4$state_code[i] == hex_v4$state_code[j]) {
      a <- suppressWarnings(st_intersection(hex_v4$geometry[i],
                                             hex_v4$geometry[j]))
      if (length(a) > 0 && poly_area(a) > 1e-8) same_ov <- same_ov + 1
    }
  }
}
cat(sprintf("Same-state overlaps: %d\n", same_ov))
stopifnot(same_ov == 0)

saveRDS(hex_v4, "data/mxmunicipio_hex_sf_v4.rds")
cat(sprintf("Saved data/mxmunicipio_hex_sf_v4.rds (%d cells)\n", nrow(hex_v4)))
