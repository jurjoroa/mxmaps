suppressPackageStartupMessages({
  library(sf); library(dplyr); library(clue); library(mxmaps)
})
sf_use_s2(FALSE)

# v3.1: refinement of the canonical v3 layout (same philosophy, better solver).
#
# Changes vs v3 (hexamap_municipios_v3.R):
#   1. Optimal assignment: greedy outermost-first nearest-cell loop replaced
#      with the Hungarian algorithm (clue::solve_LSAP) on the centroid
#      distance matrix — minimizes TOTAL municipio-to-cell displacement
#      within each state.
#   2. Balanced surplus absorption: unclaimed cells go to their nearest
#      assigned cell with a load penalty (absorbed_count * cellsize), so
#      surplus territory spreads across neighbors instead of piling onto
#      one unlucky hex (the big irregular blobs in sparse states).
#   3. Sliver/island cleanup: clipped fragments under 10% of a full hex area
#      and secondary pieces of split cells are merged into their nearest kept
#      cell instead of being assigned as real cells or silently dropped.
#      Kills the coastal dot artifacts and the island gaps.
#
# Unchanged: per-state st_make_grid sized to municipio count, shrink loop,
# clip-to-state-polygon, 1:1 municipio-to-cell, full tessellation, output
# schema (region, state_code, geometry).
#
# Output: data/mxmunicipio_hex_sf_v3_1.rds (v3 artifact untouched, for
# side-by-side comparison).

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

# --- 2. State polygons (dissolve municipios) --------------------------------
state_sf <- mxmun_sf |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_make_valid()

# --- 3. Helpers --------------------------------------------------------------
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

centroid_xy <- function(g) {
  cc <- st_coordinates(suppressWarnings(st_centroid(g)))
  c(cc[1, 1], cc[1, 2])
}

# Merge fragment geometries into their nearest kept cell (by centroid dist)
merge_fragments <- function(kept_list, frag_list) {
  if (length(frag_list) == 0) return(kept_list)
  kept_xy <- t(vapply(kept_list, centroid_xy, numeric(2)))
  for (fg in frag_list) {
    f_xy <- centroid_xy(fg)
    d <- sqrt((kept_xy[, 1] - f_xy[1])^2 + (kept_xy[, 2] - f_xy[2])^2)
    tgt <- which.min(d)
    merged <- suppressWarnings(st_union(kept_list[[tgt]], fg))
    if (inherits(merged, "sfc")) merged <- merged[[1]]
    kept_list[[tgt]] <- merged
    kept_xy[tgt, ] <- centroid_xy(kept_list[[tgt]])
  }
  kept_list
}

state_codes <- sort(unique(mxmun_sf$state_code))
results <- vector("list", length(state_codes))
names(results) <- state_codes

for (s in state_codes) {
  state_poly <- state_sf$geometry[state_sf$state_code == s]
  state_muns <- mxmun_sf[mxmun_sf$state_code == s, ]
  n_target   <- nrow(state_muns)

  st_area_deg2 <- poly_area(state_poly)
  cellsize <- sqrt(st_area_deg2 / n_target / 0.866)

  # Shrink loop: clip grid to state; split multipolygons; merge slivers.
  grid_final <- NULL
  cs <- cellsize
  for (factor in c(1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60)) {
    cs <- cellsize * factor
    suppressWarnings({
      grid <- st_make_grid(state_poly, cellsize = cs, square = FALSE)
      grid <- st_sf(geometry = grid, crs = 4326)
      grid <- grid[lengths(st_intersects(grid, state_poly)) > 0, ]
      grid_clip_geom <- st_intersection(grid$geometry, state_poly)
    })

    # Split multipolygons: largest piece stays, other pieces are fragments.
    # Slivers (< 10% of a full hex) also become fragments.
    sliver_min <- 0.10 * 0.866 * cs^2
    kept <- list(); frags <- list()
    for (k in seq_along(grid_clip_geom)) {
      g <- grid_clip_geom[k]
      if (length(g[[1]]) == 0 || st_is_empty(g) || poly_area(g) <= 1e-8) next
      if (inherits(g[[1]], "MULTIPOLYGON")) {
        pieces <- st_cast(g, "POLYGON")
        areas <- vapply(seq_along(pieces),
                        function(j) poly_area(pieces[j]), numeric(1))
        kept[[length(kept) + 1]] <- pieces[[which.max(areas)]]
        for (j in order(-areas)[-1]) frags[[length(frags) + 1]] <- pieces[[j]]
      } else if (poly_area(g) < sliver_min) {
        frags[[length(frags) + 1]] <- g[[1]]
      } else {
        kept[[length(kept) + 1]] <- g[[1]]
      }
    }
    kept <- merge_fragments(kept, frags)
    if (length(kept) >= n_target) {
      grid_final <- st_sfc(kept, crs = 4326)
      break
    }
  }
  if (is.null(grid_final)) {
    stop(sprintf("State %s: could not generate enough cells (%d municipios)",
                 s, n_target))
  }

  # Assignment: Hungarian algorithm on centroid distances (optimal total).
  suppressWarnings({
    mun_xy  <- st_coordinates(st_centroid(state_muns))
    cell_xy <- st_coordinates(st_centroid(grid_final))
  })
  cost <- sqrt(outer(mun_xy[, 1], cell_xy[, 1], "-")^2 +
               outer(mun_xy[, 2], cell_xy[, 2], "-")^2)
  assigned <- as.integer(solve_LSAP(cost))   # cell index per municipio

  # Balanced surplus absorption: nearest assigned cell + mild load penalty
  # (0.15 * cellsize per already-absorbed cell) so surplus spreads instead of
  # piling on one hex, without overriding proximity in sparse states.
  target_for_cell <- integer(length(grid_final))
  target_for_cell[assigned] <- seq_len(n_target)
  unclaimed <- which(target_for_cell == 0)
  absorbed_n <- integer(n_target)
  if (length(unclaimed) > 0) {
    d_nearest <- vapply(unclaimed, function(uc) {
      min(sqrt((cell_xy[assigned, 1] - cell_xy[uc, 1])^2 +
               (cell_xy[assigned, 2] - cell_xy[uc, 2])^2))
    }, numeric(1))
    for (uc in unclaimed[order(d_nearest)]) {
      d <- sqrt((cell_xy[assigned, 1] - cell_xy[uc, 1])^2 +
                (cell_xy[assigned, 2] - cell_xy[uc, 2])^2)
      m <- which.min(d + absorbed_n * 0.15 * cs)
      target_for_cell[uc] <- m
      absorbed_n[m] <- absorbed_n[m] + 1
    }
  }

  output_geoms <- st_sfc(lapply(seq_len(n_target), function(m) {
    cells_for_m <- which(target_for_cell == m)
    if (length(cells_for_m) == 1) return(grid_final[[cells_for_m]])
    g <- suppressWarnings(st_union(grid_final[cells_for_m]))
    if (inherits(g, "sfc")) g <- g[[1]]
    g
  }), crs = 4326)

  results[[s]] <- st_sf(
    region     = state_muns$region,
    state_code = s,
    geometry   = output_geoms,
    crs        = 4326
  )

  disp <- sqrt((cell_xy[assigned, 1] - mun_xy[, 1])^2 +
               (cell_xy[assigned, 2] - mun_xy[, 2])^2)
  cat(sprintf(
    "State %s: %d mun -> %d cells (absorbed %d, cellsize=%.3f, disp mean=%.4f max=%.4f deg)\n",
    s, n_target, length(grid_final), length(unclaimed), cs, mean(disp), max(disp)))
}

# --- 4. Combine + sanity-check ----------------------------------------------
hex_v31 <- do.call(rbind, results)
stopifnot(nrow(hex_v31) == nrow(mxmun_sf))
stopifnot(setequal(hex_v31$region, mxmun_sf$region))

suppressWarnings({ inter_idx <- st_intersects(hex_v31, hex_v31) })
cross <- 0
for (i in seq_along(inter_idx)) {
  js <- inter_idx[[i]]; js <- js[js > i]
  for (j in js) {
    if (hex_v31$state_code[i] != hex_v31$state_code[j]) {
      inter <- suppressWarnings(st_intersection(hex_v31$geometry[i],
                                                 hex_v31$geometry[j]))
      if (length(inter) > 0 && poly_area(inter) > 1e-8) cross <- cross + 1
    }
  }
}
cat(sprintf("Cross-state overlaps: %d\n", cross))

saveRDS(hex_v31, "data/mxmunicipio_hex_sf_v3_1.rds")
cat(sprintf("Saved data/mxmunicipio_hex_sf_v3_1.rds (%d cells)\n", nrow(hex_v31)))
