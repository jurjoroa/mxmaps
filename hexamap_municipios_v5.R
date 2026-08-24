suppressPackageStartupMessages({
  library(sf); library(dplyr); library(clue); library(mxmaps)
})
sf_use_s2(FALSE)

# v5: v3.1's framework (per-state grid, clip-to-boundary, full tessellation,
# Hungarian assignment) with AREA-BALANCED surplus absorption.
#
# Jorge's direction: keep cells inside the state boundary exactly as v3.1,
# but balance the space — surplus territory should flow to the tiniest
# (border-clipped) cells instead of piling onto the nearest one, so cell
# sizes are as even as possible inside each state. Sacrifice a little
# location precision to gain even space.
#
# Change vs v3.1 (hexamap_municipios_v3_1.R): unclaimed cells merge into the
# smallest-current-area ADJACENT assigned cell (space flows to small cells,
# evening sizes). v3.1 merged into the nearest cell with a mild load penalty,
# which still produced blob cells. Fallback for non-adjacent fragments
# (islands): nearest assigned cell.
#
# Output: data/mxmunicipio_hex_sf_v5.rds — same schema as v3/v3.1.

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

# Slack target can be passed: Rscript hexamap_municipios_v5.R 1.25
slack_args <- suppressWarnings(as.numeric(commandArgs(trailingOnly = TRUE)))
SLACK <- if (length(slack_args) >= 1 && !is.na(slack_args[1])) slack_args[1] else 1.2
cat(sprintf("Slack target: %.2fx\n", SLACK))

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

centroid_xy <- function(g) {
  cc <- st_coordinates(suppressWarnings(st_centroid(g)))
  c(cc[1, 1], cc[1, 2])
}

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

  grid_final <- NULL
  cs <- cellsize
  # Slack target (tunable): smaller cells clip less per cell and leave more
  # fragments to top up border cells — but too much slack makes multi-hex
  # unions the norm and sizes bimodal. Sweep externally to pick.
  slack_target <- SLACK * n_target
  for (factor in c(1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60)) {
    cs <- cellsize * factor
    suppressWarnings({
      grid <- st_make_grid(state_poly, cellsize = cs, square = FALSE)
      grid <- st_sf(geometry = grid, crs = 4326)
      grid <- grid[lengths(st_intersects(grid, state_poly)) > 0, ]
      grid_clip_geom <- st_intersection(grid$geometry, state_poly)
    })

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
    if (length(kept) >= slack_target) {
      grid_final <- st_sfc(kept, crs = 4326)
      break
    }
  }
  if (is.null(grid_final)) {
    stop(sprintf("State %s: could not generate enough cells (%d municipios)",
                 s, n_target))
  }

  # Assignment: Hungarian on centroid distances (same as v3.1)
  suppressWarnings({
    mun_xy  <- st_coordinates(st_centroid(state_muns))
    cell_xy <- st_coordinates(st_centroid(grid_final))
  })
  cost <- sqrt(outer(mun_xy[, 1], cell_xy[, 1], "-")^2 +
               outer(mun_xy[, 2], cell_xy[, 2], "-")^2)
  assigned <- as.integer(solve_LSAP(cost))

  # Area-balanced surplus absorption: each unclaimed cell goes to the
  # ADJACENT assigned cell whose final area lands closest to the ideal
  # equal share (state area / n). Space flows to small clipped cells
  # without overshooting; sizes even out inside the state.
  target_share <- st_area_deg2 / n_target
  target_for_cell <- integer(length(grid_final))
  target_for_cell[assigned] <- seq_len(n_target)
  unclaimed <- which(target_for_cell == 0)
  if (length(unclaimed) > 0) {
    cell_areas <- vapply(seq_along(grid_final),
                         function(i) poly_area(grid_final[i]), numeric(1))
    tgt_area <- cell_areas[assigned]
    hits_all <- st_intersects(grid_final, grid_final)
    for (uc in unclaimed[order(cell_areas[unclaimed])]) {
      hits <- hits_all[[uc]]
      hits <- hits[hits != uc & target_for_cell[hits] > 0]
      if (length(hits) > 0) {
        final_if <- tgt_area[target_for_cell[hits]] + cell_areas[uc]
        m <- hits[which.min(abs(final_if - target_share))]
        tgt <- target_for_cell[m]
      } else {
        # Non-adjacent fragment (island): distance-dominant, area tiebreak
        d <- sqrt((cell_xy[assigned, 1] - cell_xy[uc, 1])^2 +
                  (cell_xy[assigned, 2] - cell_xy[uc, 2])^2)
        final_if <- tgt_area + cell_areas[uc]
        score <- d / cs + 2 * abs(final_if - target_share) / target_share
        tgt <- which.min(score)
      }
      target_for_cell[uc] <- tgt
      tgt_area[tgt] <- tgt_area[tgt] + cell_areas[uc]
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

  out_areas <- vapply(seq_len(n_target),
                      function(m) poly_area(results[[s]]$geometry[m]), numeric(1))
  ratio <- out_areas / median(out_areas)
  disp <- sqrt((cell_xy[assigned, 1] - mun_xy[, 1])^2 +
               (cell_xy[assigned, 2] - mun_xy[, 2])^2)
  cat(sprintf(
    "State %s: %d mun, %d absorbed, area ratio p5=%.2f p95=%.2f, disp mean=%.4f\n",
    s, n_target, length(unclaimed),
    quantile(ratio, 0.05), quantile(ratio, 0.95), mean(disp)))
}

# --- 4. Combine + sanity-check ----------------------------------------------
hex_v5 <- do.call(rbind, results)
stopifnot(nrow(hex_v5) == nrow(mxmun_sf))
stopifnot(setequal(hex_v5$region, mxmun_sf$region))

suppressWarnings({ inter_idx <- st_intersects(hex_v5, hex_v5) })
cross <- 0
for (i in seq_along(inter_idx)) {
  js <- inter_idx[[i]]; js <- js[js > i]
  for (j in js) {
    if (hex_v5$state_code[i] != hex_v5$state_code[j]) {
      inter <- suppressWarnings(st_intersection(hex_v5$geometry[i],
                                                 hex_v5$geometry[j]))
      if (length(inter) > 0 && poly_area(inter) > 1e-8) cross <- cross + 1
    }
  }
}
cat(sprintf("Cross-state overlaps: %d\n", cross))

saveRDS(hex_v5, "data/mxmunicipio_hex_sf_v5.rds")
cat(sprintf("Saved data/mxmunicipio_hex_sf_v5.rds (%d cells)\n", nrow(hex_v5)))

# National evenness summary (cell area vs its state's ideal share)
hex_m <- st_transform(hex_v5, 6372)
areas <- as.numeric(st_area(hex_m))
share <- tapply(areas, hex_m$state_code, function(a) sum(a) / length(a))
ratio <- areas / share[hex_m$state_code]
cat(sprintf("NATIONAL evenness (slack %.2f): ratio p5=%.2f p25=%.2f median=%.2f p75=%.2f p95=%.2f\n",
            SLACK, quantile(ratio, 0.05), quantile(ratio, 0.25), median(ratio),
            quantile(ratio, 0.75), quantile(ratio, 0.95)))
