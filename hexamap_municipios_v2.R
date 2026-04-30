library(sf)
library(dplyr)
library(mxmaps)

sf_use_s2(FALSE)

# v2: snap-to-grid layout — eliminates between-state hex overlaps by
# allocating municipios to a single Mexico-wide hex lattice instead of
# per-state independent geogrid runs.

# --- 1. Reconstruct sf from the fortified mxmunicipio.map -------------------
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

# --- 2. State polygons (dissolve municipios into 32 states) -----------------
state_sf <- mxmun_sf |>
  group_by(state_code) |>
  summarise(geometry = st_union(geometry), .groups = "drop")

# --- 3. Build Mexico-wide hex lattice ---------------------------------------
# n_target controls cell density: higher → smaller hexes, more slack per state.
# Mexico's land area is ~30% of its bbox, so we oversample heavily to ensure
# each state has at least ~1.5x as many in-state cells as municipios.
n_target    <- nrow(mxmun_sf) * 5

mexico_bbox <- st_bbox(mxmun_sf)
mexico_area <- (mexico_bbox["xmax"] - mexico_bbox["xmin"]) *
               (mexico_bbox["ymax"] - mexico_bbox["ymin"])
mexico_area <- as.numeric(mexico_area)

# st_make_grid hex cellsize = flat-to-flat distance (short diameter)
# For a hex: area ≈ 0.866 * d^2, so d = sqrt(area / 0.866)
cellsize <- sqrt((mexico_area / n_target) / 0.866)

hex_lattice <- st_make_grid(
  st_as_sfc(mexico_bbox),
  cellsize = cellsize,
  square   = FALSE
) |> st_sf(geometry = _, crs = 4326)
hex_lattice$cell_id <- seq_len(nrow(hex_lattice))

# --- 4. Tag each hex cell with its state (via cell-centroid → state) --------
suppressWarnings({
  cell_centroids <- st_centroid(hex_lattice)
  cell_centroid_xy <- st_coordinates(cell_centroids)
  hex_lattice$state_code <- st_join(cell_centroids, state_sf,
                                    join = st_within)$state_code
})

n_in_state <- sum(!is.na(hex_lattice$state_code))
cat(sprintf("Lattice: %d cells, %d tagged to a state (%.0f%%), cellsize=%.3f deg\n",
  nrow(hex_lattice), n_in_state, 100 * n_in_state / nrow(hex_lattice), cellsize))

# Per-state slack diagnostic
state_capacity <- as.data.frame(table(hex_lattice$state_code, useNA = "no"))
names(state_capacity) <- c("state_code", "cells")
state_capacity$municipios <- as.integer(table(mxmun_sf$state_code)[state_capacity$state_code])
state_capacity$slack      <- state_capacity$cells - state_capacity$municipios
tight <- state_capacity[state_capacity$slack < 1, ]
if (nrow(tight) > 0) {
  cat("Tight states (cells < municipios + 1):\n")
  print(tight)
}

# --- 5. Per-state municipio → cell matching ---------------------------------
# Each municipio claims its nearest available cell. Cells available to a
# state = its in-state cells + (if shortage) nearest unclaimed externals.
# Tightest states processed first so they get the externals they need before
# sparse states would otherwise hoard them.

state_quotas <- as.integer(table(mxmun_sf$state_code))
names(state_quotas) <- names(table(mxmun_sf$state_code))
suppressWarnings({
  state_centroid_xy <- st_coordinates(st_centroid(state_sf))
  rownames(state_centroid_xy) <- state_sf$state_code
  mun_centroid_xy   <- st_coordinates(st_centroid(mxmun_sf))
})

# Tightness = needed externals (positive when in-state cells fall short)
in_state_count <- vapply(names(state_quotas), function(s) {
  sum(hex_lattice$state_code == s, na.rm = TRUE)
}, integer(1))
deficit     <- state_quotas - in_state_count
state_order <- names(deficit)[order(-deficit)]   # tightest first

hex_idx <- integer(nrow(mxmun_sf))
used    <- rep(FALSE, nrow(hex_lattice))

for (s in state_order) {
  state_pos    <- which(mxmun_sf$state_code == s)
  state_xy     <- mun_centroid_xy[state_pos, , drop = FALSE]
  state_center <- state_centroid_xy[s, ]
  n_need       <- length(state_pos)

  in_state_cand <- which(hex_lattice$state_code == s & !used)

  if (length(in_state_cand) >= n_need) {
    cand <- in_state_cand
  } else {
    # Extend with nearest unused external cells (2x slack for matching)
    n_ext  <- n_need - length(in_state_cand)
    ext    <- which(!used &
                    (is.na(hex_lattice$state_code) | hex_lattice$state_code != s))
    ext_xy <- cell_centroid_xy[ext, , drop = FALSE]
    ext_d  <- sqrt((ext_xy[, 1] - state_center[1])^2 +
                   (ext_xy[, 2] - state_center[2])^2)
    take   <- min(length(ext), n_ext * 2)
    cand   <- c(in_state_cand, ext[order(ext_d)[seq_len(take)]])
  }

  # Outermost municipio first — its position is most "load-bearing"
  ord <- order(-sqrt((state_xy[, 1] - state_center[1])^2 +
                     (state_xy[, 2] - state_center[2])^2))

  for (i in ord) {
    free_local  <- which(!used[cand])
    free_global <- cand[free_local]
    free_xy     <- cell_centroid_xy[free_global, , drop = FALSE]
    dists       <- sqrt((free_xy[, 1] - state_xy[i, 1])^2 +
                        (free_xy[, 2] - state_xy[i, 2])^2)
    pick        <- free_global[which.min(dists)]
    hex_idx[state_pos[i]] <- pick
    used[pick] <- TRUE
  }
}

# Diagnostic: how many municipios got an in-state cell vs an overflow cell
in_polygon <- !is.na(hex_lattice$state_code[hex_idx]) &
              hex_lattice$state_code[hex_idx] == mxmun_sf$state_code
cat(sprintf("In-polygon: %d / %d (%.0f%%) — rest are overflow into neighbors\n",
            sum(in_polygon), length(in_polygon),
            100 * mean(in_polygon)))

# --- 6. Build output sf -----------------------------------------------------
mxmunicipio_hex_sf_v2 <- st_sf(
  region     = mxmun_sf$region,
  state_code = mxmun_sf$state_code,
  geometry   = hex_lattice$geometry[hex_idx],
  crs        = 4326
)

# Sanity checks: every municipio gets a unique cell
stopifnot(length(unique(hex_idx)) == length(hex_idx))
stopifnot(nrow(mxmunicipio_hex_sf_v2) == nrow(mxmun_sf))

# --- 7. Save ----------------------------------------------------------------
saveRDS(mxmunicipio_hex_sf_v2, "data/mxmunicipio_hex_sf_v2.rds")
cat("Saved data/mxmunicipio_hex_sf_v2.rds —",
    nrow(mxmunicipio_hex_sf_v2), "hexes, all unique ✓\n")
