suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})
sf_use_s2(FALSE)

# Algorithmic correction of v1 hex overlaps.
#
# Strategy: asymmetric clipping. For every cross-state overlap pair, the cell
# from the state with smaller hexes (i.e. the denser state) is the "winner";
# the larger-cell neighbor has the overlap area subtracted from it. Tie-break
# by alphabetical state code so the operation is deterministic.
#
# Rationale:
#   - v1's purpose is to give each state hex space proportional to its municipio
#     count. Dense states (Oaxaca 570, Puebla 217) get small hexes; sparse
#     states (BC Sur 5, BC 6) get big hexes.
#   - When a small hex overlaps a big hex, clipping the small one would erase
#     it. Clipping the big one keeps every municipio visible at the cost of
#     one bite-mark per overlap.
#   - No movement = no cascading within-state collisions. Single pass, fast.
#
# Output: data/mxmunicipio_hex_sf_corrected.rds — same schema as v1, with
# clipped geometries for affected municipios. Original v1 file untouched.

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
hex <- readRDS("data/mxmunicipio_hex_sf.rds")
hex$state_code <- substr(hex$region, 1, 2)

# Per-cell width — proxy for "winner" priority (smaller wins)
hex$cell_w <- vapply(seq_len(nrow(hex)), function(i) {
  bb <- st_bbox(hex$geometry[i])
  as.numeric(bb["xmax"] - bb["xmin"])
}, numeric(1))

cat(sprintf("Loaded v1: %d hexes across %d states\n",
            nrow(hex), length(unique(hex$state_code))))

# --- Step 1: Detect overlapping pairs ----------------------------------------
suppressWarnings({
  inter_idx <- st_intersects(hex, hex)
})
pairs <- do.call(rbind, lapply(seq_along(inter_idx), function(i) {
  js <- inter_idx[[i]]
  js <- js[js > i]
  if (length(js) == 0) return(NULL)
  data.frame(i = i, j = js)
}))
pairs$state_i <- hex$state_code[pairs$i]
pairs$state_j <- hex$state_code[pairs$j]

# Keep only cross-state pairs with non-trivial area overlap
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

cross <- pairs[pairs$state_i != pairs$state_j, ]
cat(sprintf("Cross-state intersecting pairs: %d\n", nrow(cross)))

cross$overlap_area <- vapply(seq_len(nrow(cross)), function(k) {
  inter <- suppressWarnings(st_intersection(hex$geometry[cross$i[k]],
                                             hex$geometry[cross$j[k]]))
  if (length(inter) == 0) return(0)
  poly_area(inter)
}, numeric(1))
cross <- cross[cross$overlap_area > 1e-10, ]
cat(sprintf("Real-area overlap pairs: %d\n", nrow(cross)))

# --- Step 2: Per pair, decide loser (cell to be clipped) ---------------------
# Winner = smaller cell (denser state). Tie-break: alphabetically smaller state.
cross$w_i <- hex$cell_w[cross$i]
cross$w_j <- hex$cell_w[cross$j]
cross$loser <- ifelse(
  cross$w_i > cross$w_j, cross$i,
  ifelse(cross$w_j > cross$w_i, cross$j,
         ifelse(cross$state_i < cross$state_j, cross$j, cross$i))
)
cross$winner <- ifelse(cross$loser == cross$i, cross$j, cross$i)

# --- Step 3: For each loser cell, clip out union of all its winners ----------
losers_unique <- unique(cross$loser)
cat(sprintf("Cells to be clipped: %d\n", length(losers_unique)))

new_geoms <- hex$geometry  # start as a copy

for (l in losers_unique) {
  winner_ids <- cross$winner[cross$loser == l]
  winner_union <- suppressWarnings(
    st_union(hex$geometry[winner_ids])
  )
  clipped <- suppressWarnings(
    st_difference(hex$geometry[l], winner_union)
  )
  if (length(clipped) == 0 || st_is_empty(clipped)) {
    # Fully eaten — should not happen since winners are smaller; warn and keep
    # original geometry as fallback. Manual review needed for these.
    cat(sprintf("  WARNING: cell %d (region %s) fully eaten — keeping original\n",
                l, hex$region[l]))
    next
  }
  new_geoms[l] <- clipped
}

# --- Step 4: Build output sf -------------------------------------------------
hex_corrected <- st_sf(
  region     = hex$region,
  state_code = hex$state_code,
  geometry   = new_geoms,
  crs        = 4326
)

# Sanity: same number of municipios, same regions
stopifnot(nrow(hex_corrected) == nrow(hex))
stopifnot(all(hex_corrected$region == hex$region))

# Verify: re-run intersection check to confirm overlaps gone
suppressWarnings({
  check_idx <- st_intersects(hex_corrected, hex_corrected)
})
remaining <- 0
for (i in seq_along(check_idx)) {
  js <- check_idx[[i]]; js <- js[js > i]
  if (length(js) == 0) next
  for (j in js) {
    if (hex_corrected$state_code[i] == hex_corrected$state_code[j]) next
    inter <- suppressWarnings(
      st_intersection(hex_corrected$geometry[i], hex_corrected$geometry[j])
    )
    if (length(inter) > 0 && poly_area(inter) > 1e-10) remaining <- remaining + 1
  }
}
cat(sprintf("Remaining cross-state area overlaps after correction: %d\n", remaining))

# --- Step 5: Save ------------------------------------------------------------
saveRDS(hex_corrected, "data/mxmunicipio_hex_sf_corrected.rds")
cat(sprintf("Saved data/mxmunicipio_hex_sf_corrected.rds (%d hexes)\n",
            nrow(hex_corrected)))
