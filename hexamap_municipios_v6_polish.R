suppressPackageStartupMessages({
  library(sf); library(dplyr); library(mxmaps)
})
sf_use_s2(FALSE)

# v6-polish: surgical fix for the last starved cells after CVT repair.
# For each cell with ratio < 0.6: shave its nearest healthy neighbors
# (bisection-calibrated half-plane cuts) and union the pieces into a new
# cell for the starved seed. All in-boundary, no overlaps, deg^2 units.
# Input:  /tmp/v6repair_states.rds (15, 21) + /tmp/v6fix_good.rds (14, 18)
# Output: /tmp/v6polish_states.rds

setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")

fixed <- readRDS("/tmp/v6fix_good.rds")
rep   <- readRDS("/tmp/v6repair_states.rds")

sanitize <- function(x) {
  if (!inherits(x$geometry, "sfc")) {
    g <- st_sfc(lapply(x$geometry, function(gg) if (inherits(gg, "sfc")) gg[[1]] else gg),
                crs = 4326)
    x <- st_as_sf(data.frame(region = x$region, state_code = x$state_code),
                  geometry = g)
  }
  x
}
cells <- rbind(sanitize(fixed), sanitize(rep))

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

shave_piece <- function(donor, toward_xy, want) {
  dc <- st_coordinates(st_centroid(donor))[1, ]
  u <- toward_xy - dc
  u <- u / sqrt(sum(u^2))
  lo <- 0; hi <- 3 * sqrt(poly_area(donor))
  best <- NULL; best_a <- 0
  for (k in 1:12) {
    d <- (lo + hi) / 2
    Q <- dc + d * u
    perp <- c(-u[2], u[1]); big <- 10
    quad <- rbind(Q - big * perp, Q + big * perp,
                  Q + big * perp + 10 * u, Q - big * perp + 10 * u,
                  Q - big * perp)
    hp <- st_polygon(list(quad))
    piece <- suppressWarnings(st_intersection(donor, st_sfc(hp, crs = 4326)))
    a <- if (length(piece) == 0) 0 else poly_area(piece)
    if (a > best_a) { best <- piece; best_a <- a }
    # piece area shrinks as d grows: too small -> move cut back (hi = d)
    if (a < want) hi <- d else lo <- d
    if (abs(a - want) / want < 0.05) break
  }
  best
}

polish <- function(sf_obj) {
  for (pass in 1:6) {
    areas <- vapply(seq_len(nrow(sf_obj)),
                    function(i) poly_area(sf_obj$geometry[i]), numeric(1))
    share <- sum(areas) / nrow(sf_obj)
    ratio <- areas / share
    worst <- which.min(ratio)
    if (ratio[worst] > 0.6) break
    need <- share - areas[worst]
    g <- sf_obj$geometry
    wcc <- st_coordinates(suppressWarnings(st_centroid(g[worst])))[1, ]
    d_cents <- st_coordinates(suppressWarnings(st_centroid(g)))
    dd <- sqrt((d_cents[, 1] - wcc[1])^2 + (d_cents[, 2] - wcc[2])^2)
    dd[worst] <- Inf
    healthy <- which(ratio > 0.75)
    donors <- healthy[order(dd[healthy])][seq_len(min(6, length(healthy)))]
    pieces <- list(); got <- 0
    for (dn in donors) {
      want <- min((need - got) / length(donors), poly_area(g[dn]) - 0.72 * share)
      if (want <= 0) next
      piece <- shave_piece(g[dn], wcc, want)
      if (is.null(piece) || length(piece) == 0) next
      ap <- poly_area(piece)
      if (ap <= 0) next
      pieces[[length(pieces) + 1]] <- piece[[1]]
      sf_obj$geometry[dn] <- st_sfc(suppressWarnings(st_difference(g[dn], piece)),
                                    crs = 4326)[[1]]
      got <- got + ap
      if (got >= 0.9 * need) break
    }
    if (length(pieces) == 0) break
    sf_obj$geometry[worst] <- st_sfc(suppressWarnings(st_union(st_as_sfc(pieces))),
                                     crs = 4326)[[1]]
    cat(sprintf("  polished: worst ratio %.2f -> fed by %d donor pieces\n",
                ratio[worst], length(pieces)))
  }
  sf_obj
}

out <- do.call(rbind, lapply(split(cells, cells$state_code), polish))
stopifnot(nrow(out) == nrow(cells))
stopifnot(inherits(out$geometry, "sfc"))
saveRDS(out, "/tmp/v6polish_states.rds")

for (s in unique(out$state_code)) {
  sub <- out[out$state_code == s, ]
  a <- vapply(seq_len(nrow(sub)), function(i) poly_area(sub$geometry[i]),
              numeric(1))
  r <- a / (sum(a) / length(a))
  cat(sprintf("State %s POLISHED: ratio p5=%.2f p95=%.2f min=%.2f max=%.2f\n",
              s, quantile(r, .05), quantile(r, .95), min(r), max(r)))
}
