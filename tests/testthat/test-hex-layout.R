# -----------------------------------------------------------------------
# Purpose : Property tests for the committed municipio hex-layout
#           artifacts in data/. These are the objects that can break
#           SILENTLY: a generator script is scratch code, but the .rds it
#           writes is what every plot downstream reads. Each test below
#           corresponds to a bug class that actually shipped in this
#           project (see the BUG: tags).
# Author  : Jorge Roa
# Version : 1.0.0
# Date    : 2026-08-26
# Usage   : testthat::test_file("tests/testthat/test-hex-layout.R")
#           or run as part of devtools::test() / R CMD check.
# Output  : testthat results. Missing artifacts skip(), they never fail,
#           because the artifacts are scratch output and are not
#           guaranteed to be present in a source checkout.
# -----------------------------------------------------------------------
#
# ALL TOLERANCES BELOW WERE MEASURED on
# data/mxmunicipio_hex_sf_v8_8.rds before being written down. The
# measured value is quoted in a comment next to every threshold. Two
# invariants from the ideal spec DO NOT hold exactly and are asserted at
# their real, measured tolerance instead of an invented one:
#   * "no cell overlaps another": 3522 pairs do overlap. Almost all are
#     numerical slivers (median 2.5e-14 of the smaller cell) but 12 pairs
#     exceed 0.1% and the worst reaches 75.3% of the smaller cell.
#   * "the union of a state's cells covers the state polygon": offshore
#     ISLAND parts of a state are never covered (Cozumel, Isla Mujeres,
#     Cedros, Angel de la Guarda...). Coverage is therefore asserted
#     against each state's MAINLAND part, with a looser separate bound on
#     the whole (island-inclusive) polygon.

testthat::skip_if_not_installed("sf")

library(sf)

suppressMessages(sf_use_s2(FALSE))

# Adjacency of hex cells must be measured from GEOMETRY, with a small
# snapping distance, not from cell-centre distance.
# BUG: a site-distance proxy once stood in for true cell adjacency and
# the optimiser happily improved the proxy. st_touches() is NOT usable
# here: ~46% of visually abutting cells share an infinitesimal sliver
# overlap, which makes st_touches() FALSE. Measured cohesion is flat at
# 0.7801 for snap distances from 0.001 m to 1000 m, so 1 m is safe.
adj_snap_m <- 1

# ---- fixtures ---------------------------------------------------------

hex_path <- function(file_name) {
  candidates <- c(
    test_path("..", "..", "data", file_name),
    file.path("data", file_name),
    system.file("data", file_name, package = "mxmaps")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

read_hex <- function(file_name) {
  path <- hex_path(file_name)
  if (is.na(path)) {
    skip(paste0("layout artifact not present: data/", file_name))
  }
  readRDS(path)
}

# Cache the expensive derived objects: 2469 true polygons, their
# state unions and the true intra-state adjacency list.
.cache <- new.env(parent = emptyenv())

true_municipios <- function() {
  if (!is.null(.cache$mx)) return(.cache$mx)
  skip_if_not_installed("mxmaps")
  df <- local({
    # mxmunicipio.map is LAZY-loaded, so it is not in the namespace directly;
    # asNamespace() finds it only if the package happens to be attached. Use
    # utils::data(), which works in a clean session and under R CMD check.
    e <- new.env()
    utils::data("mxmunicipio.map", package = "mxmaps", envir = e)
    get("mxmunicipio.map", envir = e)
  })
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
    if (length(rings) == 0) return(st_polygon())
    st_polygon(c(
      lapply(rings[!vapply(rings, `[[`, logical(1), "hole")],
             `[[`, "coords"),
      lapply(rings[vapply(rings, `[[`, logical(1), "hole")],
             `[[`, "coords")
    ))
  })
  mx <- st_make_valid(st_sf(
    region = regions,
    geometry = st_sfc(polys, crs = 4326)
  ))
  mx$state_code <- substr(mx$region, 1, 2)
  .cache$mx <- st_transform(mx, 6372)
  .cache$mx
}

# Per-state real outline, plus its largest connected part ("mainland").
state_outlines <- function() {
  if (!is.null(.cache$states)) return(.cache$states)
  mx <- true_municipios()
  codes <- sort(unique(mx$state_code))
  full <- list()
  main <- list()
  for (sc in codes) {
    poly <- st_union(st_geometry(mx)[mx$state_code == sc])
    parts <- st_cast(poly, "POLYGON")
    full[[sc]] <- poly
    main[[sc]] <- parts[which.max(as.numeric(st_area(parts)))]
  }
  .cache$states <- list(codes = codes, full = full, main = main)
  .cache$states
}

# Pairs (i, j), i < j, of TRUE municipio polygons that are adjacent and
# in the same state. Indices are into true_municipios().
true_intra_state_pairs <- function() {
  if (!is.null(.cache$pairs)) return(.cache$pairs)
  mx <- true_municipios()
  hits <- st_intersects(mx, mx)
  out <- do.call(rbind, lapply(seq_len(nrow(mx)), function(i) {
    js <- hits[[i]]
    js <- js[js > i & mx$state_code[js] == mx$state_code[i]]
    if (length(js) == 0) return(NULL)
    data.frame(i = i, j = js)
  }))
  .cache$pairs <- out
  out
}

# Cells projected to EPSG:6372, geometry repaired, area in km2.
hex_metric <- function(hex) {
  g <- st_make_valid(st_geometry(st_transform(hex, 6372)))
  list(g = g, area_km2 = as.numeric(st_area(g)) / 1e6)
}

overlap_table <- function(g, area_km2) {
  ov <- st_overlaps(g, g, sparse = TRUE)
  pr <- do.call(rbind, lapply(seq_along(ov), function(i) {
    js <- ov[[i]][ov[[i]] > i]
    if (length(js) == 0) return(NULL)
    data.frame(i = i, j = js)
  }))
  if (is.null(pr)) {
    return(data.frame(i = integer(), j = integer(),
                      km2 = numeric(), frac = numeric()))
  }
  pr$km2 <- vapply(seq_len(nrow(pr)), function(k) {
    as.numeric(st_area(st_intersection(g[[pr$i[k]]], g[[pr$j[k]]]))) / 1e6
  }, numeric(1))
  pr$frac <- pr$km2 / pmin(area_km2[pr$i], area_km2[pr$j])
  pr
}

# Share of true intra-state municipio adjacencies that survive as
# adjacent CELLS, measured from geometry at a given snap distance.
cohesion_at <- function(g, pairs, dist_m) {
  nb <- st_is_within_distance(g, g, dist = dist_m, sparse = TRUE)
  mean(mapply(function(i, j) j %in% nb[[i]], pairs$i, pairs$j))
}

# The artifact the artifact-level invariants are asserted against. This pointed
# at v8_8 for three versions after v8_8 stopped being the canonical layout, so
# "run the tests before shipping an artifact" was checking a superseded file --
# the one that renders 38 cells under a quarter of their state mean and 8
# fabricated placeholder discs, including Ecatepec (pop 1.65M) at 1.40 km2
# against its true 157.8 km2. A guard aimed at the wrong artifact is not a guard.
canonical <- "mxmunicipio_hex_sf_v11_3.rds"

n_municipios <- 2469L

# Every artifact that must at least satisfy the cheap structural and
# identity invariants. The invalid-geometry column is a RATCHET on
# measured baselines, not an ideal: v8_5 really does ship 23 invalid
# cells and v8_6 one, and a test that pretends otherwise would be a
# lie. What must never happen is that number going UP.
# min_rel_area is the same kind of ratchet for the "cell vanished to
# zero area" bug, expressed as the smallest cell area divided by the
# mean cell area of its own state. v8_6 ships a genuinely collapsed
# cell (21053, Puebla: 0.061 km2, 3.885e-04 of the Puebla mean) and
# v12 a smaller one relative to v8_8. Baselines are the measured
# values rounded down; they may only ever move UP.
known_artifacts <- data.frame(
  file = c("mxmunicipio_hex_sf_v8_5.rds", "mxmunicipio_hex_sf_v8_6.rds",
           "mxmunicipio_hex_sf_v8_8.rds", "mxmunicipio_hex_sf_v11.rds",
           "mxmunicipio_hex_sf_v12.rds", "mxmunicipio_hex_sf_v11_2.rds",
           "mxmunicipio_hex_sf_v11_3.rds"),
  max_invalid = c(23L, 1L, 0L, 0L, 0L, 0L, 0L),
  # measured: 0.03325, 3.885e-04, 0.007773, 0.4417, 0.002124, 0.449362,
  #           0.450817
  min_rel_area = c(0.03, 3e-4, 5e-3, 0.4, 2e-3, 0.44, 0.45),
  stringsAsFactors = FALSE
)

# ---- pure unit test: the strtoi trap ----------------------------------

test_that("strtoi(\"08\") is NA and as.integer(\"08\") is 8", {
  # BUG: set.seed(strtoi(state_code)) crashed on state "08"
  # (Chihuahua) and "09" (CDMX). strtoi() defaults to base = 0L, which
  # means "guess from the prefix", and a leading zero means OCTAL, so
  # "08" and "09" are not valid digits and come back NA. Every INEGI
  # state code in this package is zero-padded to two characters, so this
  # is a trap for 2 of the 32 states only, which is exactly why it
  # survived testing. Use as.integer(), or strtoi(x, base = 10L).
  expect_true(is.na(strtoi("08")))
  expect_true(is.na(strtoi("09")))
  expect_identical(as.integer("08"), 8L)
  expect_identical(as.integer("09"), 9L)
  expect_identical(strtoi("08", base = 10L), 8L)
  # And the reason it is not caught by the other 30 states:
  expect_identical(strtoi("07"), 7L)
})

# ---- structural + identity invariants, every artifact ----------------

for (row_i in seq_len(nrow(known_artifacts))) {
  art <- known_artifacts$file[row_i]
  max_invalid <- known_artifacts$max_invalid[row_i]
  min_rel_area <- known_artifacts$min_rel_area[row_i]

  test_that(paste(art, "has one valid, positive-area cell per municipio"), {
    hex <- read_hex(art)
    skip_if_not_installed("mxmaps")

    # Invariant 1: exactly one cell per municipio, regions unique and
    # equal to the true municipio set.
    expect_s3_class(hex, "sf")
    expect_identical(nrow(hex), n_municipios)
    expect_identical(anyDuplicated(hex$region), 0L)
    expect_true(all(nchar(hex$region) == 5L))
    mx <- true_municipios()
    expect_setequal(hex$region, mx$region)

    # Invariant 7, the index-confusion guard.
    # BUG: a cell index was used where a municipio index was meant,
    # silently permuting the assignment. A permutation that crossed a
    # state boundary is caught right here.
    expect_identical(hex$state_code, substr(hex$region, 1, 2))

    # Invariant 6: per-state cell count == true municipio count.
    true_n <- table(mx$state_code)
    hex_n <- table(hex$state_code)
    expect_setequal(names(hex_n), names(true_n))
    expect_equal(as.integer(hex_n[names(true_n)]),
                 as.integer(true_n))

    # Invariant 2: geometry non-empty, strictly positive area.
    # BUG: cells vanished to zero area during an iteration.
    expect_false(any(st_is_empty(hex)))
    m <- hex_metric(hex)
    expect_true(all(m$area_km2 > 0))
    # Measured smallest cell: v8_8 1.228 km2, v8_6 0.061 km2,
    # v12 0.330 km2, v8_5 8.161 km2, v11 55.84 km2.
    expect_gt(min(m$area_km2), 1e-3)
    # Relative floor, which is what "vanished" really means: no cell may
    # collapse relative to the mean cell in its own state. This is where
    # v8_6 fails a uniform 1e-3 floor, so the bound is the measured
    # per-artifact ratchet from known_artifacts, not one invented number.
    rel <- m$area_km2 / ave(m$area_km2, hex$state_code, FUN = mean)
    expect_gt(min(rel), min_rel_area)

    # Validity ratchet on the RAW stored geometry (before st_make_valid).
    expect_lte(sum(!st_is_valid(st_transform(hex, 6372))), max_invalid)
  })
}

# ---- geometric invariants, canonical artifact ------------------------

test_that("canonical cells do not overlap beyond the measured tolerance", {
  skip_on_cran()
  hex <- read_hex(canonical)
  m <- hex_metric(hex)
  ov <- overlap_table(m$g, m$area_km2)

  # Invariant 3. Cells still touch numerically -- 3811 pairs share a sliver of
  # area at floating-point scale -- but the REAL collisions are gone. The
  # v8_8-era defect this used to pin (pair 21083/29031, 75.3% of the smaller
  # cell covered, plus six more above 20%) does not exist in the canonical
  # artifact: the worst pair is 7e-04 and NO pair exceeds 1e-3. Bounds are the
  # v11_3 measurements plus a small multiple; they may only ever move DOWN.
  expect_lt(sum(ov$km2) / sum(m$area_km2), 2e-6)   # measured 4.813e-07
  expect_lt(median(ov$frac), 1e-12)                # measured 2.462e-14
  expect_identical(sum(ov$frac > 1e-3), 0L)        # measured 0 pairs
  expect_lt(max(ov$frac), 5e-3)                    # measured 7.0e-04
})

test_that("canonical keeps every cell inside its own state polygon", {
  skip_on_cran()
  hex <- read_hex(canonical)
  m <- hex_metric(hex)
  st <- state_outlines()

  # Invariant 4. Measured on v11_3: worst state spills 7.50e-06 of its own
  # area, national 1.07e-06 -- two orders tighter than the v8_8 figures this
  # used to pin (1.26e-04 worst, 6.92e-06 national).
  out_frac <- vapply(st$codes, function(sc) {
    cells <- st_union(m$g[hex$state_code == sc])
    outside <- st_difference(cells, st$full[[sc]])
    num <- if (length(outside) == 0) 0 else as.numeric(st_area(outside))
    num / as.numeric(st_area(st$full[[sc]]))
  }, numeric(1))

  expect_lt(max(out_frac), 5e-5)                   # measured 7.502e-06
  expect_lt(sum(out_frac * vapply(st$codes, function(sc)
    as.numeric(st_area(st$full[[sc]])), numeric(1))) /
      sum(vapply(st$codes, function(sc)
        as.numeric(st_area(st$full[[sc]])), numeric(1))), 5e-6)  # 1.065e-06
})

test_that("canonical covers each state's real outline, islands included", {
  skip_on_cran()
  hex <- read_hex(canonical)
  m <- hex_metric(hex)
  st <- state_outlines()

  # Invariant 5, the project's hard "real outline preserved" constraint.
  # THE ISLAND DEFECT IS CURED, and these bounds now say so. Under v8_8 the
  # whole-polygon miss was 2.0e-02 and entirely offshore: state 23 lost Cozumel
  # (475.8 km2) and Isla Mujeres (45.6), state 02 lost Cedros (926.1) and Angel
  # de la Guarda (350.3), state 03 had 12 of its 13 parts uncovered. The v11.1
  # fix -- a largest_piece() call at OUTPUT was discarding island fragments --
  # removed all of it. Mainland and whole-polygon figures now agree to three
  # significant figures because there is no island gap left to separate them.
  # Holding the old 0.03 tolerance would let every island vanish again unseen.
  areas <- vapply(st$codes, function(sc)
    as.numeric(st_area(st$main[[sc]])), numeric(1))
  unc <- vapply(st$codes, function(sc) {
    cells <- st_union(m$g[hex$state_code == sc])
    gap <- st_difference(st$main[[sc]], cells)
    if (length(gap) == 0) 0 else as.numeric(st_area(gap))
  }, numeric(1))

  expect_lt(max(unc / areas), 5e-5)            # measured 9.100e-06 (26)
  expect_lt(median(unc / areas), 2e-6)         # measured 2.463e-07
  expect_lt(sum(unc) / sum(areas), 1e-5)       # measured 1.820e-06

  areas_full <- vapply(st$codes, function(sc)
    as.numeric(st_area(st$full[[sc]])), numeric(1))
  unc_full <- vapply(st$codes, function(sc) {
    cells <- st_union(m$g[hex$state_code == sc])
    gap <- st_difference(st$full[[sc]], cells)
    if (length(gap) == 0) 0 else as.numeric(st_area(gap))
  }, numeric(1))
  expect_lt(max(unc_full / areas_full), 5e-5)  # measured 9.038e-06 (26)
  expect_lt(sum(unc_full) / sum(areas_full), 1e-5)  # measured 1.818e-06
})

# ---- the assignment is the right assignment -------------------------

test_that("canonical assignment is not a permutation of the cells", {
  skip_on_cran()
  hex <- read_hex(canonical)
  skip_if_not_installed("mxmaps")
  m <- hex_metric(hex)
  mx <- true_municipios()

  # BUG: a cell index used where a municipio index was meant silently
  # PERMUTED the assignment. Invariant 7 only catches permutations that
  # cross a state boundary; a within-state shuffle passes it untouched
  # and every area/coverage invariant above as well, because the set of
  # cells is unchanged. The only thing that detects it is displacement:
  # each cell must sit near ITS OWN municipio. Compared against a
  # seeded within-state random shuffle of the same cells, which is
  # exactly the bug's signature.
  k <- match(hex$region, mx$region)
  expect_false(anyNA(k))
  true_xy <- st_coordinates(st_centroid(st_geometry(mx)))[k, ]
  cell_xy <- st_coordinates(st_centroid(m$g))
  disp_km <- sqrt(rowSums((cell_xy - true_xy)^2)) / 1000

  set.seed(1L)
  perm <- seq_len(nrow(hex))
  for (sc in unique(hex$state_code)) {
    idx <- which(hex$state_code == sc)
    # NOT sample(idx): sample() on a length-1 vector permutes 1:idx
    # instead of the vector, the same shape of trap as strtoi("08").
    perm[idx] <- idx[sample.int(length(idx))]
  }
  disp_perm_km <- sqrt(rowSums((cell_xy[perm, ] - true_xy)^2)) / 1000

  # Displacement is the one term allowed to move between versions -- the pair
  # repair moves sites deliberately -- so headroom here is wider than on the
  # geometric invariants, but not 2x wider.
  expect_lt(mean(disp_km), 30)          # measured 24.75 km
  expect_lt(median(disp_km), 24)        # measured 18.61 km
  # Measured ratio 0.179: the real assignment is 5.6x better than a
  # within-state shuffle. Approaching 0.30 means it has been scrambled.
  expect_lt(mean(disp_km) / mean(disp_perm_km), 0.30)
})

test_that("canonical cohesion is measured from geometry, not from centres", {
  skip_on_cran()
  hex <- read_hex(canonical)
  skip_if_not_installed("mxmaps")
  m <- hex_metric(hex)
  pairs <- true_intra_state_pairs()
  expect_identical(nrow(pairs), 6288L)

  # BUG: a swap delta ignored the neighbours of the swapped pair, so the
  # optimiser improved a number that was not cohesion. A floor on
  # geometry-derived cohesion is what would have caught that: the score
  # can only be reported as improved if the real adjacency improved.
  # THE LOAD-BEARING GUARD. This floor sat at 0.74 -- a v8_8 measurement left
  # behind when canonical moved to v11_3 -- against an actual 0.8290, so the
  # suite permitted an 8.9 pp regression: twice the whole v11.2 -> v11.3 gain
  # (+4.39 pp), and enough to accept a silent revert to v11.2 (0.7851).
  coh <- cohesion_at(m$g, pairs, adj_snap_m)
  expect_gt(coh, 0.82)                  # measured 0.8290

  # BUG: a site-distance proxy stood in for the true cell adjacency.
  # Two guards. First, the geometric number must be a PLATEAU in the
  # snap distance: a proxy masquerading as adjacency would drift with
  # the threshold. Measured 0.7801 at 0.001 m and at 1 m, identical.
  coh_tight <- cohesion_at(m$g, pairs, 0.001)
  expect_lt(abs(coh - coh_tight), 0.01)

  # Second, a cell-centre distance rule is NOT equivalent to geometric
  # adjacency, and this asserts they disagree on a material number of
  # pairs so nobody "simplifies" the scorer into the proxy again.
  xy <- st_coordinates(st_centroid(m$g))
  thr <- sqrt(mean(m$area_km2) * 1e6) * 1.3
  proxy <- sqrt(rowSums((xy[pairs$i, ] - xy[pairs$j, ])^2)) < thr
  nb <- st_is_within_distance(m$g, m$g, dist = adj_snap_m, sparse = TRUE)
  geom <- mapply(function(i, j) j %in% nb[[i]], pairs$i, pairs$j)
  expect_gt(sum(proxy != geom), 100L)
})
