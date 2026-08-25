suppressPackageStartupMessages({
  library(sf); library(dplyr); library(clue); library(mxmaps)
})
sf_use_s2(FALSE)

# v6-fix: rescue pass for states whose v6 run produced empty/fallback cells
# (14 Jalisco x1, 15 Edomex x3, 21 Puebla x13 — seeds swallowed by neighbors,
# usually near holes/coasts). Same weighted CVT as v6, plus:
#   - empty seed rescue: seed pulled to the nearest point of the state and
#     given a strong weight boost until it claims real territory
#   - Nayarit (18) included: one borderline small cell (0.53)
#
# Usage: Rscript hexamap_municipios_v6_fix.R 14 15 18 21
# Output: /tmp/v6fix_states.rds (only these states)

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

args <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)))
ONLY <- if (length(args) >= 1 && !any(is.na(args))) args else c(14, 15, 18, 21)

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

state_codes <- as.character(ONLY)
results <- vector("list", length(state_codes))
names(results) <- state_codes

for (s in state_codes) {
  state_poly <- state_sf$geometry[state_sf$state_code == s]
  state_muns <- mxmun_sf[mxmun_sf$state_code == s, ]
  n_target   <- nrow(state_muns)
  st_area_deg2 <- poly_area(state_poly)
  cellsize <- sqrt(st_area_deg2 / n_target / 0.866)

  if (length(state_poly[[1]]) > 1 || inherits(state_poly[[1]], "MULTIPOLYGON")) {
    pieces <- st_cast(state_poly, "POLYGON")
    p_areas <- vapply(seq_along(pieces), function(j) poly_area(pieces[j]),
                      numeric(1))
    state_clip <- st_sfc(pieces[[which.max(p_areas)]], crs = 4326)
  } else {
    state_clip <- st_sfc(state_poly[[1]], crs = 4326)
  }

  seeds <- NULL
  cs <- cellsize
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
  cost <- sqrt(outer(mun_xy[, 1], seeds[, 1], "-")^2 +
               outer(mun_xy[, 2], seeds[, 2], "-")^2)
  assigned <- as.integer(solve_LSAP(cost))
  P <- seeds[assigned, , drop = FALSE]
  rownames(P) <- NULL

  bb <- st_bbox(state_poly)
  pad <- 0.5 * max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"]) + cs
  sq <- matrix(c(bb["xmin"] - pad, bb["ymin"] - pad,
                 bb["xmax"] + pad, bb["ymin"] - pad,
                 bb["xmax"] + pad, bb["ymax"] + pad,
                 bb["xmin"] - pad, bb["ymax"] + pad),
               ncol = 2, byrow = TRUE)
  target <- st_area_deg2 / n_target
  w <- rep(0, n_target)
  alpha <- 0.15
  n_iter <- 250
  for (it in seq_len(n_iter)) {
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
    w[ok] <- w[ok] + alpha * (target - areas[ok])
    w <- w - mean(w)
    P[ok, ] <- P[ok, ] + 0.5 * (cents[ok, ] - P[ok, ])
    # RESCUE: empty seeds pulled onto the state + strong weight boost
    empty <- which(!ok)
    if (length(empty) > 0) {
      for (i in empty) {
        np <- st_nearest_points(st_sfc(st_point(P[i, ]), crs = 4326),
                                state_clip)
        cc <- st_coordinates(np)[2, ]
        P[i, ] <- c(cc[1], cc[2])
        w[i] <- w[i] + 0.3 * target
      }
      w <- w - mean(w)
    }
    if (it %% 20 == 0) {
      cat(sprintf("  state %s iter %d: ratio in [%.2f, %.2f], empty=%d\n",
                  s, it, min(areas[ok]) / target, max(areas[ok]) / target,
                  length(empty)))
    }
    if (length(empty) == 0 && all(abs(areas[ok] / target - 1) < 0.05)) break
  }

  output_geoms <- st_sfc(lapply(seq_len(n_target), function(i) {
    if (!is.null(cells[[i]])) return(cells[[i]])
    st_buffer(st_sfc(st_point(P[i, ]), crs = 4326), 0.1 * cs)[[1]]
  }), crs = 4326)

  results[[s]] <- st_sf(
    region     = state_muns$region,
    state_code = s,
    geometry   = output_geoms,
    crs        = 4326
  )
  fa <- as.numeric(st_area(st_transform(results[[s]], 6372)))
  ratio <- fa / (sum(fa) / n_target)
  cat(sprintf("State %s FIXED: ratio p5=%.2f p95=%.2f min=%.2f max=%.2f\n",
              s, quantile(ratio, 0.05), quantile(ratio, 0.95),
              min(ratio), max(ratio)))
}

fixed <- do.call(rbind, results)
stopifnot(nrow(fixed) == sum(mxmun_sf$state_code %in% state_codes))
saveRDS(fixed, "/tmp/v6fix_states.rds")
cat(sprintf("Saved /tmp/v6fix_states.rds (%d cells, states: %s)\n",
            nrow(fixed), paste(state_codes, collapse = ",")))
