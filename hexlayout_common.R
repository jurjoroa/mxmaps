# ---------------------------------------------------------------------------
# Purpose: shared machinery for the municipio layout scripts.
#
#          WHY THIS FILE EXISTS. The mxmunicipio.map -> sf builder was
#          copy-pasted into 52 root-level scripts, and largest_piece() was
#          defined 11 times across 6 files (twice in some). Every bug in this
#          project so far survived because of that duplication: a cell index
#          used as a municipio index, a swap delta that ignored the neighbours
#          of the swapped pair, a site-distance proxy standing in for the true
#          cell adjacency, and strtoi("08") returning NA. Each had to be found
#          and fixed independently per copy. One definition, one fix.
#
#          These bodies are lifted VERBATIM from hexamap_municipios_v8_9.R --
#          the newest version, which carries all four fixes -- so behaviour is
#          unchanged, with one deliberate exception noted below.
#
#          DELIBERATE FIX: cvt_run() previously read `pull_mu` out of the
#          caller's environment by lexical scope leak (it was a loop variable
#          in the state loop, resolved at the global level). That works only by
#          accident of R's scoping and breaks silently if the caller renames or
#          drops the variable. It is now an explicit parameter defaulting to
#          0.25, the historical value.
#
# Author:  Claude (Opus 5) for Jorge Roa
# Version: shared
# Date:    26-August-2026
# Usage:   source("hexlayout_common.R")  from a root-level scratch script.
#          Do NOT source this from inside R/ -- root is scratch, R/ is package.
# ---------------------------------------------------------------------------

# Build the true municipio polygons from the package's ggplot-style long table.
# Returns an sf in EPSG:4326 with columns region, state_code.
build_mxmun_sf <- function() {
  utils::data("mxmunicipio.map", package = "mxmaps", envir = environment())
  df <- get("mxmunicipio.map", envir = environment())
  regions <- unique(df$region)
  polys <- lapply(regions, function(r) {
    sub <- df[df$region == r, ]
    groups <- unique(sub$group)
    rings <- lapply(groups, function(g) {
      pts <- sub[sub$group == g, ]
      pts <- pts[order(pts$order), ]
      co <- cbind(pts$long, pts$lat)
      if (nrow(co) < 4) return(NULL)
      if (!all(co[1, ] == co[nrow(co), ])) co <- rbind(co, co[1, ])
      list(coords = co, hole = any(pts$hole))
    })
    rings <- Filter(Negate(is.null), rings)
    if (length(rings) == 0) return(sf::st_polygon())
    sf::st_polygon(c(
      lapply(rings[!sapply(rings, `[[`, "hole")], `[[`, "coords"),
      lapply(rings[ sapply(rings, `[[`, "hole")], `[[`, "coords")))
  })
  out <- sf::st_make_valid(sf::st_sf(
    region = regions, geometry = sf::st_sfc(polys, crs = 4326)))
  out$state_code <- substr(out$region, 1, 2)
  out
}

# Per-state cell-size context: true areas, municipio counts, the theoretical
# cell size and its ratio to the national one. cs_ratio drives both the gamma
# rule and the position anchor in the v8 line.
state_size_table <- function(mxmun_sf) {
  km2 <- mxmun_sf |>
    sf::st_transform(6372) |>
    dplyr::mutate(true_km2 = as.numeric(sf::st_area(geometry)) / 1e6) |>
    sf::st_drop_geometry() |>
    dplyr::select(region, state_code, true_km2)
  tab <- km2 |>
    dplyr::group_by(state_code) |>
    dplyr::summarise(state_km2 = sum(true_km2), n = dplyr::n(),
                     .groups = "drop") |>
    dplyr::mutate(cs_km = sqrt(state_km2 / n / 0.866))
  cs_nat <- sqrt(sum(tab$state_km2) / sum(tab$n) / 0.866)
  tab$cs_ratio <- tab$cs_km / cs_nat
  # the v8.3+ derived gamma rule (hand overrides are NOT applied here; a caller
  # that wants them must say so explicitly, so they stay visible)
  tab$gamma_auto <- pmax(pmin(0.18 * (tab$cs_ratio - 0.75), 0.45), 0)
  list(per_municipio = km2, per_state = tab, cs_national = cs_nat)
}

# strtoi() defaults to base = 0L, which reads a leading zero as OCTAL, so
# strtoi("08") and strtoi("09") are NA. State codes here are zero-padded
# ("01".."32"), so ALWAYS use this instead. This exact trap crashed
# hexamap_municipios_v8_5.R at state 08 and is why that version cannot
# reproduce its own committed artifact.
state_code_int <- function(s) as.integer(s)

# ---- geometry primitives -------------------------------------------------

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

# ---- power diagram / CVT -------------------------------------------------

cvt_run <- function(P, w, tgt, alpha, n_iter, damp, tol,
                    state_clip, sq, mun_xy, anchor, true_nb = NULL,
                    pull_mu = 0.25) {
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

# ---- assignment search --------------------------------------------------

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

viol_base <- function(assigned, seeds, st_pairs, cs) {
  sum(vapply(seq_len(length(assigned)), function(k) {
    if (length(st_pairs[[k]]) == 0) return(0L)
    sum(vapply(st_pairs[[k]], function(j)
      as.integer(sqrt(sum((seeds[assigned[k], ] - seeds[assigned[j], ])^2)) >
                   1.2 * cs), 0L))
  }, 0L)) %/% 2
}
