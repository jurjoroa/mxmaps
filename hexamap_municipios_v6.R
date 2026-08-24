suppressPackageStartupMessages({
  library(sf); library(dplyr); library(clue); library(mxmaps)
})
sf_use_s2(FALSE)

# v6: equal-area cells via weighted centroidal Voronoi (power diagram).
#
# Jorge's direction: maximize bin-size similarity inside each state — borrow
# space so all cells are even, while staying within the state boundary.
# Instead of clipping a hex lattice and absorbing surplus (v3.1/v5), v6
# REDRAWS the internal boundaries:
#   1. Assign municipios to hex-lattice positions (Hungarian, as before).
#   2. Treat the n assigned lattice centers as seeds; compute each seed's
#      power cell (convex, via half-plane clipping) intersected with the
#      state polygon.
#   3. Iterate: w_i += alpha * (target_area - area_i)  (area correction)
#               seed_i <- centroid(cell_i)             (Lloyd step)
#      Seeds spread until every cell ~ state_area / n.
# Interior cells stay near-hexagonal; border cells conform to the true
# boundary; nothing leaves the state; sizes are even by construction.
#
# Output: data/mxmunicipio_hex_sf_v6.rds — same schema.

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

# Optional: restrict to specific states, e.g. Rscript hexamap_municipios_v6.R 14 15 21
args <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
ONLY <- if (length(args) >= 1 && !any(is.na(args))) args else NULL

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

# Sutherland-Hodgman: clip convex polygon (matrix n x 2) by a*x <= b
clip_hp <- function(pts, a, b) {
  n <- nrow(pts)
  out <- vector("list", n + 1)
  m <- 0
  for (i in seq_len(n)) {
    cur <- pts[i, ]; nxt <- pts[(i %% n) + 1, ]
    dc <- sum(a * cur); dn <- sum(a * nxt)
    if (dc <= b) {
      m <- m + 1; out[[m]] <- cur
      if (dn > b) {
        t <- (b - dc) / (dn - dc)
        m <- m + 1; out[[m]] <- cur + t * (nxt - cur)
      }
    } else if (dn <= b) {
      t <- (b - dc) / (dn - dc)
      m <- m + 1; out[[m]] <- cur + t * (nxt - cur)
    }
  }
  if (m < 3) return(NULL)
  do.call(rbind, out[seq_len(m)])
}

state_codes <- sort(unique(mxmun_sf$state_code))
if (!is.null(ONLY)) state_codes <- ONLY

results <- vector("list", length(state_codes))
names(results) <- state_codes

for (s in state_codes) {
  state_poly <- state_sf$geometry[state_sf$state_code == s]
  state_muns <- mxmun_sf[mxmun_sf$state_code == s, ]
  n_target   <- nrow(state_muns)
  st_area_deg2 <- poly_area(state_poly)
  cellsize <- sqrt(st_area_deg2 / n_target / 0.866)

  # CVT clips against the largest polygon piece (mainland): islands make
  # centroid feedback oscillate, holes (Edomex/CDMX, Puebla/Tlaxcala) are
  # kept. Island municipios get mainland cells near their true position.
  if (length(state_poly[[1]]) > 1 || inherits(state_poly[[1]], "MULTIPOLYGON")) {
    pieces <- st_cast(state_poly, "POLYGON")
    p_areas <- vapply(seq_along(pieces), function(j) poly_area(pieces[j]),
                      numeric(1))
    state_clip <- st_sfc(pieces[[which.max(p_areas)]], crs = 4326)
  } else {
    state_clip <- st_sfc(state_poly[[1]], crs = 4326)
  }

  # 1. Lattice candidates: shrink until >= n centers near/inside the state
  seeds <- NULL
  cs <- cellsize
  for (factor in c(1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60)) {
    cs <- cellsize * factor
    grid <- st_make_grid(state_poly, cellsize = cs, square = FALSE)
    cents <- st_centroid(grid, warn = FALSE)
    suppressWarnings({
      d_poly <- as.numeric(st_distance(st_transform(cents, 6372),
                                       st_transform(state_poly, 6372)))
    })
    inside <- lengths(st_within(cents, state_poly)) > 0
    ord <- order(-as.integer(inside), d_poly)
    if (length(grid) >= n_target) {
      pick <- ord[seq_len(n_target)]
      seeds <- st_coordinates(cents)[pick, , drop = FALSE]
      break
    }
  }
  if (is.null(seeds)) stop(sprintf("State %s: no lattice fit", s))

  # 2. Assignment: Hungarian municipio -> lattice seed
  suppressWarnings(mun_xy <- st_coordinates(st_centroid(state_muns)))
  cost <- sqrt(outer(mun_xy[, 1], seeds[, 1], "-")^2 +
               outer(mun_xy[, 2], seeds[, 2], "-")^2)
  assigned <- as.integer(solve_LSAP(cost))
  # reorder seeds to municipio order
  P <- seeds[assigned, , drop = FALSE]
  rownames(P) <- NULL

  # 3. Weighted CVT: equalize cell areas inside the state polygon
  bb <- st_bbox(state_poly)
  pad <- 0.5 * max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"]) + cs
  sq <- matrix(c(bb["xmin"] - pad, bb["ymin"] - pad,
                 bb["xmax"] + pad, bb["ymin"] - pad,
                 bb["xmax"] + pad, bb["ymax"] + pad,
                 bb["xmin"] - pad, bb["ymax"] + pad),
               ncol = 2, byrow = TRUE)
  target <- st_area_deg2 / n_target
  w <- rep(0, n_target)
  alpha <- 0.6
  n_iter <- 60
  for (it in seq_len(n_iter)) {
    # fresh neighbor graph each iteration (seeds move)
    d_mat <- as.matrix(dist(P))
    knn <- apply(d_mat, 1, function(r) order(r)[2:min(ncol(d_mat), 17)])
    cells <- vector("list", n_target)
    areas <- numeric(n_target)
    cents <- matrix(NA_real_, n_target, 2)
    for (i in seq_len(n_target)) {
      pts <- sq
      for (j in knn[, i]) {
        a <- 2 * (P[j, ] - P[i, ])
        b2 <- sum(P[j, ]^2) - sum(P[i, ]^2) + w[i] - w[j]
        pts <- clip_hp(pts, a, b2)
        if (is.null(pts)) break
      }
      if (is.null(pts)) {
        cells[[i]] <- NULL; areas[i] <- 0; next
      }
      poly <- st_polygon(list(rbind(pts, pts[1, ])))
      clipped <- suppressWarnings(st_intersection(st_sfc(poly, crs = 4326),
                                                   state_clip))
      if (length(clipped) == 0) {
        cells[[i]] <- NULL; areas[i] <- 0; next
      }
      areas[i] <- poly_area(clipped)
      cc <- st_coordinates(suppressWarnings(st_centroid(clipped)))
      cents[i, ] <- c(cc[1, 1], cc[1, 2])
      cells[[i]] <- clipped[[1]]
    }
    # weight correction + Lloyd move (guard empty cells)
    ok <- areas > 0
    w[ok] <- w[ok] + alpha * (target - areas[ok])
    w <- w - mean(w)
    P[ok, ] <- cents[ok, ]
    if (it %% 10 == 0) {
      cat(sprintf("  state %s iter %d: area ratio in [%.2f, %.2f]\n",
                  s, it, min(areas[ok]) / target, max(areas[ok]) / target))
    }
    if (all(abs(areas[ok] / target - 1) < 0.02)) break
  }

  output_geoms <- st_sfc(lapply(seq_len(n_target), function(i) {
    if (!is.null(cells[[i]])) return(cells[[i]])
    # rare fallback: isolated seed never claimed area — tiny disk at seed
    st_buffer(st_point(P[i, ]), 0.1 * cs)
  }), crs = 4326)
  results[[s]] <- st_sf(
    region     = state_muns$region,
    state_code = s,
    geometry   = output_geoms,
    crs        = 4326
  )
  ratio <- areas[ok] / target
  disp <- sqrt((P[, 1] - mun_xy[, 1])^2 + (P[, 2] - mun_xy[, 2])^2)
  cat(sprintf("State %s: %d cells, ratio p5=%.2f p95=%.2f, disp mean=%.4f\n",
              s, n_target, quantile(ratio, 0.05), quantile(ratio, 0.95),
              mean(disp)))
}

hex_v6 <- do.call(rbind, results)
stopifnot(nrow(hex_v6) == length(unique(hex_v6$region)))

if (is.null(ONLY)) {
  stopifnot(nrow(hex_v6) == nrow(mxmun_sf))
  stopifnot(setequal(hex_v6$region, mxmun_sf$region))
  saveRDS(hex_v6, "data/mxmunicipio_hex_sf_v6.rds")
  cat(sprintf("Saved data/mxmunicipio_hex_sf_v6.rds (%d cells)\n", nrow(hex_v6)))
} else {
  saveRDS(hex_v6, sprintf("/tmp/v6_state_%s.rds", paste(ONLY, collapse = "_")))
  cat(sprintf("Saved single-state rds to /tmp (%d cells)\n", nrow(hex_v6)))
}
