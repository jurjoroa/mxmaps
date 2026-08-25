suppressPackageStartupMessages({
  library(sf); library(dplyr); library(clue); library(mxmaps)
})
sf_use_s2(FALSE)

# v8: reality-weighted cell sizes (per-state gamma) on top of v7.1.
#
# Jorge: sparse states have spare space — distribute it so cell sizes
# reflect true municipio sizes (Ensenada/Los Cabos/Mulege larger), while
# keeping balance. Dense states stay equal-area.
#
# Mechanism: weighted CVT. Target area per cell =
#   state_area * (true_area_i ^ gamma_s) / sum(true_area ^ gamma_s)
# gamma_s per state = auto rule from sparsity (cs_state / cs_national),
# overridden by the table below for named states. Final targets clamped
# to [0.45, 2.2] x equal share (renormalized) so bins stay balanced.
#
# Usage: Rscript hexamap_municipios_v8.R [state_code ...]
# Output: data/mxmunicipio_hex_sf_v8.rds (full) or /tmp/v8_states.rds.

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

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

# --- True municipio areas (km2, planar) --------------------------------------
mxmun_km2 <- mxmun_sf |>
  st_transform(6372) |>
  mutate(true_km2 = as.numeric(st_area(geometry)) / 1e6) |>
  st_drop_geometry()

# --- Per-state gamma: auto rule + named-state overrides ----------------------
state_areas <- mxmun_km2 |>
  group_by(state_code) |>
  summarise(state_km2 = sum(true_km2), n = dplyr::n(), .groups = "drop") |>
  mutate(cs_km = sqrt(state_km2 / n / 0.866))
cs_national <- sqrt(sum(state_areas$state_km2) / sum(state_areas$n) / 0.866)
state_areas$gamma_auto <- pmax(pmin(0.18 * (state_areas$cs_km / cs_national - 0.75),
                                    0.45), 0)
# Named-state overrides (Jorge's list: more space -> match reality harder)
gamma_override <- c("02" = 0.45,   # Baja California (Ensenada dominates)
                    "03" = 0.45,   # BCS (Mulege/Los Cabos)
                    "05" = 0.30,   # Coahuila (Ocampo/Sierra Mojada huge)
                    "08" = 0.30,   # Chihuahua
                    "26" = 0.25,   # Sonora
                    "10" = 0.25,   # Durango
                    "07" = 0.12,   # Chiapas (complex; mild)
                    "04" = 0.15,   # Campeche
                    "24" = 0.12,   # San Luis Potosi
                    "01" = 0.10,   # Aguascalientes
                    "09" = 0.10,   # CDMX (Milpa Alta vs Cuajimalpa)
                    "06" = 0.12)   # Colima (Manzanillo vs Minatitlan)
state_areas$gamma <- state_areas$gamma_auto
state_areas$gamma[match(names(gamma_override), state_areas$state_code)] <-
  gamma_override

cat("Per-state gamma:\n")
print(as.data.frame(state_areas[order(-state_areas$gamma),
      c("state_code", "n", "cs_km", "gamma_auto", "gamma")]))

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

local_search <- function(assigned, pos, st_pairs, mun_xy, lambda, touch_d,
                         max_pass = 25) {
  n <- length(assigned)
  disp_of <- function(k) sqrt(sum((mun_xy[k, ] - pos[assigned[k], ])^2))
  viol_of <- function(k) {
    if (length(st_pairs[[k]]) == 0) return(0L)
    sum(vapply(st_pairs[[k]], function(j)
      as.integer(sqrt(sum((pos[assigned[k], ] - pos[assigned[j], ])^2)) >
                   touch_d), 0L))
  }
  improved <- TRUE; pass <- 0
  while (improved && pass < max_pass) {
    improved <- FALSE; pass <- pass + 1
    for (a in sample(seq_len(n))) {
      if (length(st_pairs[[a]]) == 0) next
      cand <- unique(unlist(lapply(st_pairs[[a]], function(b) {
        nb <- assigned[b]
        dd <- sqrt((pos[, 1] - pos[nb, 1])^2 + (pos[, 2] - pos[nb, 2])^2)
        which(dd < 4 * touch_d)
      })))
      cand <- setdiff(cand, assigned[a])
      best_dJ <- -1e-9; best_c <- NA
      da_old <- disp_of(a); va_old <- viol_of(a)
      for (c in cand) {
        dc_old <- disp_of(c); vc_old <- viol_of(c)
        assigned[c(a, c)] <- c(assigned[c], assigned[a])
        dJ <- (disp_of(a) - da_old) + (disp_of(c) - dc_old) +
          lambda * (viol_of(a) - va_old + viol_of(c) - vc_old)
        assigned[c(a, c)] <- c(assigned[c], assigned[a])
        if (!is.na(dJ) && dJ < best_dJ) { best_dJ <- dJ; best_c <- c }
      }
      if (!is.na(best_c)) {
        cc <- best_c
        assigned[c(a, cc)] <- c(assigned[cc], assigned[a])
        improved <- TRUE
      }
    }
  }
  list(assigned = assigned,
       viol = sum(vapply(seq_len(n), viol_of, 0L)) %/% 2)
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

  # --- Weighted targets: true-area^gamma, clamped around equal share -------
  ta <- mxmun_km2$true_km2[mun_idx]
  w_t <- ta^gamma_s
  share <- st_area_deg2 / n_target
  tgt <- w_t / sum(w_t) * st_area_deg2
  for (r in 1:3) {
    tgt <- pmax(pmin(tgt, 2.2 * share), 0.45 * share)
    tgt <- tgt / sum(tgt) * st_area_deg2
  }
  tgt_ratio_range <- paste0(round(min(tgt / share), 2), "-",
                            round(max(tgt / share), 2))

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
  viol_base <- sum(vapply(seq_len(n_target), function(k) {
    if (length(st_pairs[[k]]) == 0) return(0L)
    sum(vapply(st_pairs[[k]], function(j)
      as.integer(sqrt(sum((seeds[assigned[k], ] - seeds[assigned[j], ])^2)) >
                   1.2 * cs), 0L))
  }, 0L)) %/% 2
  ls1 <- local_search(assigned, seeds, st_pairs, mun_xy, lambda, 1.2 * cs)
  assigned <- ls1$assigned
  P <- seeds[assigned, , drop = FALSE]
  rownames(P) <- NULL

  bb <- st_bbox(state_poly)
  pad <- 0.5 * max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"]) + cs
  sq <- matrix(c(bb["xmin"] - pad, bb["ymin"] - pad,
                 bb["xmax"] + pad, bb["ymin"] - pad,
                 bb["xmax"] + pad, bb["ymax"] + pad,
                 bb["xmin"] - pad, bb["ymax"] + pad), ncol = 2, byrow = TRUE)
  w <- rep(0, n_target)
  for (phase in list(list(0.6, 60, 1.0, 0.05),
                     list(0.05, 10, 0.4, 1.0),
                     list(0.3, 20, 0.6, 0.05))) {
    alpha <- phase[[1]]; n_iter <- phase[[2]]; damp <- phase[[3]]; tol <- phase[[4]]
    for (it in seq_len(n_iter)) {
      d_mat <- as.matrix(dist(P))
      knn <- apply(d_mat, 1, function(r) order(r)[2:min(ncol(d_mat), 17)])
      cells <- vector("list", n_target); areas <- numeric(n_target)
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
      w[ok] <- pmax(pmin(w[ok] + alpha * (tgt[ok] - areas[ok]), 1.5 * max(tgt)),
                    -1.5 * max(tgt))
      w <- w - mean(w)
      P[ok, ] <- P[ok, ] + damp * (cents[ok, ] - P[ok, ])
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
  }

  cell_pos <- cents
  na_rows <- which(!complete.cases(cell_pos))
  if (length(na_rows) > 0) cell_pos[na_rows, ] <- P[na_rows, ]
  g_cells <- st_transform(st_sfc(cells, crs = 4326), 6372)
  cadj <- st_is_within_distance(g_cells, g_cells, dist = 1, sparse = TRUE)
  spacing <- sqrt(mean(tgt))
  ls2 <- local_search(seq_len(n_target), cell_pos, st_pairs, mun_xy,
                      0.75 * spacing, 1.3 * spacing)
  assigned2 <- ls2$assigned

  largest_piece <- function(g) {
    if (!inherits(g[[1]], "MULTIPOLYGON")) return(g[[1]])
    pieces <- st_cast(g, "POLYGON")
    ar <- vapply(seq_along(pieces), function(j) poly_area(pieces[j]), numeric(1))
    pieces[[which.max(ar)]]
  }
  output_geoms <- st_sfc(lapply(seq_len(n_target), function(i) {
    ci <- assigned2[i]
    if (!is.null(cells[[ci]])) return(largest_piece(st_sfc(cells[[ci]], crs = 4326)))
    st_buffer(st_sfc(st_point(cell_pos[ci, ]), crs = 4326),
              0.05 * sqrt(tgt[ci]))[[1]]
  }), crs = 4326)

  results[[s]] <- st_sf(region = state_muns$region, state_code = s,
                        geometry = output_geoms, crs = 4326)

  fa <- as.numeric(st_area(st_transform(results[[s]], 6372)))
  ratio_eq <- fa / (sum(fa) / n_target)
  ratio_tg <- fa / (sum(fa) * tgt / sum(tgt))
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
  disp <- sqrt((cell_pos[assigned2, 1] - mun_xy[, 1])^2 +
               (cell_pos[assigned2, 2] - mun_xy[, 2])^2)
  cat(sprintf(
    "State %s: g=%.2f tgt[%s] | area-vs-target p5=%.2f p95=%.2f | adj %.0f%%->%.0f%% | disp=%.3f\n",
    s, gamma_s, tgt_ratio_range, quantile(ratio_tg, .05), quantile(ratio_tg, .95),
    100 * viol_base / max(tot, 1), 100 * pres / max(tot, 1), mean(disp)))
}

hex_v8 <- do.call(rbind, results)
stopifnot(nrow(hex_v8) == sum(mxmun_sf$state_code %in% state_codes))
stopifnot(inherits(hex_v8$geometry, "sfc"))

if (is.null(ONLY)) {
  stopifnot(nrow(hex_v8) == nrow(mxmun_sf))
  saveRDS(hex_v8, "data/mxmunicipio_hex_sf_v8.rds")
  cat(sprintf("Saved data/mxmunicipio_hex_sf_v8.rds (%d cells)\n", nrow(hex_v8)))
} else {
  saveRDS(hex_v8, "/tmp/v8_states.rds")
  cat(sprintf("Saved /tmp/v8_states.rds (%d cells)\n", nrow(hex_v8)))
}
