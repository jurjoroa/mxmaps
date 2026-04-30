suppressPackageStartupMessages({
  library(sf); library(dplyr); library(mxmaps)
})
sf_use_s2(FALSE)

# v3: per-state hex grid, sized to fit municipio count, clipped to state polygon.
#
# Pattern follows the references' methodology (st_make_grid + clip-to-boundary
# + spatial assignment) but applied per-state so each state gets a cell size
# proportional to its municipio density.
#
# Properties:
#   - Each municipio gets exactly one cell (1:1 mapping preserved from v1)
#   - Cells inside state polygon → no overflow, clear state identity
#   - State polygon fully tessellated → no gaps
#   - No overlaps (each state's grid is independent and clipped to its border)
#   - Distortion: cells near state borders get clipped to the polygon edge,
#     becoming hex-with-bite shapes. Dense states have small cells, sparse
#     states have big cells.
#
# Trade-off vs v1: gives up pure hex shape near borders for full coverage and
# clean state identity.
#
# Output: data/mxmunicipio_hex_sf_v3.rds — same schema as v1.

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

# --- 3. Per-state hex grid + assignment -------------------------------------
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
results <- vector("list", length(state_codes))
names(results) <- state_codes

for (s in state_codes) {
  state_poly <- state_sf$geometry[state_sf$state_code == s]
  state_muns <- mxmun_sf[mxmun_sf$state_code == s, ]
  n_target   <- nrow(state_muns)

  st_area_deg2 <- poly_area(state_poly)

  # Initial cellsize: theoretical fit if state were a perfect hex packing.
  # 0.866 = ratio of regular hex area to squared flat-to-flat width.
  # Real states have irregular shape, so we'll shrink iteratively.
  cellsize <- sqrt(st_area_deg2 / n_target / 0.866)

  # Generate grid; shrink cellsize until we have enough cells (with slack)
  grid_clip <- NULL
  for (factor in c(1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60)) {
    cs <- cellsize * factor
    suppressWarnings({
      grid <- st_make_grid(state_poly, cellsize = cs, square = FALSE)
      grid <- st_sf(geometry = grid, crs = 4326)
      grid <- grid[lengths(st_intersects(grid, state_poly)) > 0, ]
      grid_clip_geom <- st_intersection(grid$geometry, state_poly)
    })
    # Filter empty / non-polygon results
    keep <- vapply(seq_along(grid_clip_geom),
                   function(i) {
                     g <- grid_clip_geom[i]
                     if (length(g[[1]]) == 0) return(FALSE)
                     if (st_is_empty(g)) return(FALSE)
                     poly_area(g) > 1e-8
                   }, logical(1))
    grid_clip_geom <- grid_clip_geom[keep]
    if (length(grid_clip_geom) >= n_target) {
      grid_clip <- grid_clip_geom
      break
    }
  }
  if (is.null(grid_clip)) {
    stop(sprintf("State %s: could not generate enough cells (%d municipios)",
                 s, n_target))
  }

  # Recast all to POLYGON (some clips produce MULTIPOLYGON for islands etc.)
  grid_clip <- st_sfc(unlist(lapply(seq_along(grid_clip), function(i) {
    g <- grid_clip[i]
    if (inherits(g[[1]], "MULTIPOLYGON")) {
      # Take the largest piece — small islands get absorbed into nearest cell
      pieces <- st_cast(g, "POLYGON")
      areas <- vapply(seq_along(pieces),
                      function(j) poly_area(pieces[j]), numeric(1))
      list(pieces[which.max(areas)][[1]])
    } else {
      list(g[[1]])
    }
  }), recursive = FALSE), crs = 4326)

  # Assignment: outermost municipios first (they anchor the borders)
  suppressWarnings({
    state_center <- st_coordinates(st_centroid(state_poly))[1, ]
    mun_xy   <- st_coordinates(st_centroid(state_muns))
    cell_xy  <- st_coordinates(st_centroid(grid_clip))
  })
  mun_dist <- sqrt((mun_xy[, 1] - state_center[1])^2 +
                   (mun_xy[, 2] - state_center[2])^2)
  process_order <- order(-mun_dist)  # outermost first

  used     <- rep(FALSE, length(grid_clip))
  assigned <- integer(n_target)

  for (m in process_order) {
    avail   <- which(!used)
    cell_d  <- sqrt((cell_xy[avail, 1] - mun_xy[m, 1])^2 +
                    (cell_xy[avail, 2] - mun_xy[m, 2])^2)
    pick    <- avail[which.min(cell_d)]
    assigned[m] <- pick
    used[pick]  <- TRUE
  }

  # Absorb unclaimed surplus cells: each goes to its nearest assigned cell
  # so the state's territory is fully tessellated by the n_target output cells.
  unclaimed <- which(!used)
  target_for_cell <- integer(length(grid_clip))
  target_for_cell[assigned] <- seq_len(n_target)
  if (length(unclaimed) > 0) {
    assigned_cell_xy <- cell_xy[assigned, , drop = FALSE]
    for (uc in unclaimed) {
      d <- sqrt((assigned_cell_xy[, 1] - cell_xy[uc, 1])^2 +
                (assigned_cell_xy[, 2] - cell_xy[uc, 2])^2)
      target_for_cell[uc] <- which.min(d)
    }
  }

  output_geoms <- st_sfc(lapply(seq_len(n_target), function(m) {
    cells_for_m <- which(target_for_cell == m)
    if (length(cells_for_m) == 1) {
      grid_clip[[cells_for_m]]
    } else {
      g <- suppressWarnings(st_union(grid_clip[cells_for_m]))
      g[[1]]
    }
  }), crs = 4326)

  results[[s]] <- st_sf(
    region     = state_muns$region,
    state_code = s,
    geometry   = output_geoms,
    crs        = 4326
  )

  cat(sprintf("State %s: %d municipios → %d cells (%d absorbed surplus, cellsize=%.3f)\n",
              s, n_target, length(grid_clip), length(unclaimed), cs))
}

# --- 4. Combine + sanity-check ---------------------------------------------
hex_v3 <- do.call(rbind, results)
stopifnot(nrow(hex_v3) == nrow(mxmun_sf))
stopifnot(setequal(hex_v3$region, mxmun_sf$region))

# Verify: zero cross-state overlaps (impossible by construction, but check)
suppressWarnings({ inter_idx <- st_intersects(hex_v3, hex_v3) })
cross <- 0
for (i in seq_along(inter_idx)) {
  js <- inter_idx[[i]]; js <- js[js > i]
  for (j in js) {
    if (hex_v3$state_code[i] != hex_v3$state_code[j]) {
      inter <- suppressWarnings(st_intersection(hex_v3$geometry[i],
                                                 hex_v3$geometry[j]))
      if (length(inter) > 0 && poly_area(inter) > 1e-8) cross <- cross + 1
    }
  }
}
cat(sprintf("Cross-state overlaps: %d\n", cross))

saveRDS(hex_v3, "data/mxmunicipio_hex_sf_v3.rds")
cat(sprintf("Saved data/mxmunicipio_hex_sf_v3.rds (%d cells)\n", nrow(hex_v3)))
