# ---------------------------------------------------------------------------
# Purpose: Sweep the area-equalisation exponent gamma and measure the
#          LEGIBILITY <-> FIDELITY frontier the project has never plotted.
#
#          Cell target area is proportional to true_area^gamma. gamma = 0 gives
#          perfectly equal cells (maximum legibility, the point of the layout);
#          gamma = 1 gives true areas (a plain map, no legibility gain). Every
#          released version so far picked one gamma from a formula. Nobody has
#          measured what the exchange rate actually is, so nobody could say
#          whether the formula sits on the good part of the curve.
#
#          Engine: ot_municipios_v11.R, unmodified, driven through its
#          V11_GAMMA / V11_KNN env vars and its positional state list. It pins
#          sites at the true municipio centroids and solves cell areas by
#          damped Newton on the semi-discrete optimal-transport dual, so gamma
#          is the ONLY thing that moves between runs -- no CVT, no assignment
#          search, no restarts. That is what makes a clean sweep possible.
#
# LEGIBILITY -- the definition used here, and why.
#   Primary   leg_ratio = area_p5 / area_p95 within the state.
#             The reader's problem is the SMALLEST cells: a municipio rendered
#             far below its neighbours cannot be colour-matched to a legend,
#             hovered, or pointed at. leg_ratio is the robust smallest-to-
#             largest area ratio; 1.00 = perfectly equal cells, 0 = one cell
#             invisible next to another. p5/p95 rather than min/max because a
#             single sliver clipped by the state outline should not define the
#             legibility of 217 cells -- min/max is reported too (leg_minmax).
#   Tail      tiny_frac = share of cells under 0.25x their state's mean cell
#             area. Same threshold evaluate_states.R now uses for its
#             "legibility tail", kept identical so the two scorers agree.
#   Absolute  leg_print = share of cells whose equal-area-circle diameter is
#             at least 2 mm when the state is printed with its longer bbox
#             side at 120 mm (a one-state page figure), and leg_print_nat, the
#             same test when the WHOLE COUNTRY is printed 180 mm wide. 2 mm is
#             about the floor at which a filled patch can still be matched to a
#             legend swatch by eye and hit by a cursor. These are the only
#             measures with real units; leg_ratio is the one used for the
#             frontier because it is scale-free and comparable across states.
#   Diagnostic leg_ratio_req is the legibility the raw gamma rule ASKS for,
#             before the engine's [0.45, 2.2] x share target clamp. Comparing
#             it with leg_ratio_clamp (what the clamp leaves) and leg_ratio
#             (what came out) shows how much of the sweep the clamp is eating;
#             clamp_bind_frac is the share of cells pinned to a bound.
#   Costs     cohesion_pct (share of true intra-state municipio adjacencies
#             surviving as touching cells, the project's existing metric) and
#             disp_km (mean cell-centroid displacement from the true municipio
#             centroid). Both computed here, not read from the engine, so every
#             gamma is scored by identical code.
#
#   VALIDATION Because "legibility" has never been validated here, the script
#             also measures what gamma does to READING REAL DATA. A choropleth's
#             perceived average is the INK-weighted (cell-area-weighted) mean of
#             the variable; the honest per-municipio average is the unweighted
#             mean. gamma = 0 makes those identical by construction; gamma = 1
#             hands the visual average to whichever municipios happen to be
#             large. The bias is reported for four genuinely different 2020
#             census variables (indigenous-language share, Afro-Mexican share,
#             sex ratio, log10 population density).
#
#   KNEE      Kneedle / maximum-distance-to-chord on the (legibility,
#             cohesion) curve after min-max normalising both over the swept
#             range: the gamma whose point lies furthest from the straight line
#             joining the two extreme gammas. Objective, no tuned constant.
#
# Author:  Claude (Opus 5) for Jorge Roa
# ISS:     -
# Version: gamma-frontier v1
# Date:    2026-08-26
#
# Usage:   Rscript gamma_frontier.R [sweep|metrics|plot|all]
#          GF_CACHE=<dir>   per-(gamma,state) engine output cache
#                           (default /tmp/gamma_frontier_cache)
#          GF_KNN=<k>       starting neighbour count handed to the engine
#                           (default 60 -- see the note on Puebla below)
# Output:  data/gamma_frontier.csv          tidy, one row per (gamma, state)
#          state_pages/gamma_frontier.png   the frontier plot
#          stdout                           panel aggregate + knee statement
#
# NOTE ON V11_KNN. The engine's default V11_KNN = 40 does not converge on
# Puebla: at gamma >= 0.15 damped Newton stalls on iteration 1 and returns the
# unweighted Voronoi with max area error 502-878%. V11_KNN = 60 converges in 5
# iterations to 0.1%. The sweep therefore runs at 60 and RETRIES any state whose
# reported maxerr exceeds 5% at 90 then 120, recording knn_used and maxerr_pct
# per row so an unconverged cell is visible rather than silently averaged in.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(ggplot2); library(scales)
  library(patchwork); library(mxmaps)
})
sf_use_s2(FALSE)
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
options(dplyr.summarise.inform = FALSE)

STAGE <- { a <- commandArgs(trailingOnly = TRUE)
           if (length(a) >= 1) a[1] else "all" }
CACHE <- { v <- Sys.getenv("GF_CACHE")
           if (nzchar(v)) v else "/tmp/gamma_frontier_cache" }
KNN0  <- { v <- Sys.getenv("GF_KNN"); if (nzchar(v)) as.integer(v) else 60L }
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
dir.create("state_pages", showWarnings = FALSE)

GAMMAS <- c(0, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0)
# panel spanning the whole density range, 03/06/13/21/30 mandatory;
# 26 SON and 09 CDMX bracket the cell-size extremes, 17 MOR fills the middle.
# 20 OAX (n = 570) deliberately excluded: too slow for a 7-point sweep.
PANEL <- c("03", "06", "09", "13", "17", "21", "26", "30")
PRINT_MM     <- 120   # longer side of the printed one-state figure
PRINT_MM_NAT <- 180   # longer side of a printed national figure
NAT_SPAN_KM  <- 2900  # Mexico's E-W extent, for the national print scale
READ_MM      <- 2     # smallest legible filled-patch diameter
TINY_FRAC  <- 0.25  # "legibility tail" threshold, matches evaluate_states.R

# ---- true municipio polygons (the construction every script here uses) ----
build_true <- function() {
  data("mxmunicipio.map")
  df <- mxmunicipio.map
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
    st_polygon(c(lapply(rings[!sapply(rings, `[[`, "hole")], `[[`, "coords"),
                 lapply(rings[ sapply(rings, `[[`, "hole")], `[[`, "coords")))
  })
  mx <- st_make_valid(st_sf(region = regions,
                            geometry = st_sfc(polys, crs = 4326)))
  mx$state_code <- substr(mx$region, 1, 2)
  mx
}

# ---- stage 1: drive the engine once per (gamma, state) --------------------
cache_path <- function(g, s) file.path(CACHE, sprintf("g%s_s%s.rds",
                                                     formatC(g, format = "f",
                                                             digits = 2), s))
meta_path  <- function(g, s) sub("\\.rds$", "_meta.rds", cache_path(g, s))

run_one <- function(g, s) {
  cp <- cache_path(g, s); mp <- meta_path(g, s)
  if (file.exists(cp) && file.exists(mp)) return(readRDS(mp))
  knn <- KNN0
  repeat {
    out <- suppressWarnings(system2(
      "/usr/local/bin/Rscript",
      c("ot_municipios_v11.R", as.character(as.integer(s))),
      env = c(sprintf("V11_GAMMA=%g", g), sprintf("V11_KNN=%d", knn)),
      stdout = TRUE, stderr = TRUE))
    ln <- grep("^state ", out, value = TRUE)
    if (!length(ln))
      stop(sprintf("engine gave no state line for g=%g s=%s", g, s))
    num <- function(pat) {
      m <- regmatches(ln[1], regexpr(pat, ln[1]))
      as.numeric(gsub("[^0-9.]", "", m))
    }
    meta <- data.frame(
      gamma = g, state_code = s, knn_used = knn,
      newton_it  = num("newton it=\\s*[0-9]+"),
      maxerr_pct = num("maxerr=\\s*[0-9.]+"),
      eng_cohesion = num("COHESION\\s+[0-9.]+"))
    if (meta$maxerr_pct <= 5 || knn >= 120) break
    knn <- if (knn < 90) 90L else 120L
    cat(sprintf("    retry g=%.2f s=%s at KNN=%d (maxerr %.1f%%)\n",
                g, s, knn, meta$maxerr_pct))
  }
  file.copy("/tmp/v11_states.rds", cp, overwrite = TRUE)
  saveRDS(meta, mp)
  meta
}

if (STAGE %in% c("sweep", "all")) {
  # The engine always writes /tmp/v11_states.rds, so two sweeps running at
  # once silently race on it and can cache a layout under the wrong gamma.
  # Refuse to start rather than produce numbers nobody can trust.
  lock <- file.path(CACHE, "sweep.lock")
  if (file.exists(lock)) {
    pid <- suppressWarnings(readLines(lock, warn = FALSE)[1])
    alive <- length(suppressWarnings(system2(
      "ps", c("-p", pid, "-o", "pid="), stdout = TRUE, stderr = FALSE))) > 0
    if (isTRUE(alive))
      stop(sprintf("a sweep is already running (pid %s)", pid))
    cat(sprintf("  stale lock from dead pid %s - taking over\n", pid))
  }
  writeLines(as.character(Sys.getpid()), lock)
  cat("== stage 1: sweep ==\n")
  for (g in GAMMAS) for (s in PANEL) {
    t0 <- Sys.time()
    m <- run_one(g, s)
    cat(sprintf(
      "  g=%.2f s=%s knn=%d it=%d maxerr=%6.2f%% eng_coh=%5.1f%% (%.0fs)\n",
      g, s, m$knn_used, m$newton_it, m$maxerr_pct, m$eng_cohesion,
      as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  unlink(lock)
}

# ---- stage 2: score every cached layout with identical code ---------------
if (STAGE %in% c("metrics", "all")) {
  cat("== stage 2: metrics ==\n")
  mx <- build_true()
  mx_m <- st_transform(mx, 6372)
  true_km2 <- as.numeric(st_area(mx_m)) / 1e6
  tcn <- st_coordinates(suppressWarnings(
    st_centroid(mx_m, of_largest_polygon = TRUE)))
  adj <- st_intersects(mx, mx)

  # four genuinely different real variables, for the ink-bias validation
  vt <- df_mxmunicipio_2020
  vt <- data.frame(
    region = vt$region,
    pct_indig = 100 * vt$indigenous_language / vt$pop,
    pct_afro  = 100 * vt$afromexican / vt$pop,
    sex_ratio = 100 * vt$pop_male / vt$pop_female,
    log_dens  = log10(vt$pop / true_km2[match(vt$region, mx$region)]))
  VARS <- c("pct_indig", "pct_afro", "sex_ratio", "log_dens")

  score_state <- function(g, s) {
    cp <- cache_path(g, s)
    if (!file.exists(cp)) return(NULL)
    hex <- readRDS(cp)
    hex <- hex[hex$state_code == s, ]
    hm <- st_transform(hex, 6372)
    a_km2 <- as.numeric(st_area(hm)) / 1e6
    n <- nrow(hex)
    share <- sum(a_km2) / n
    r_eq <- a_km2 / share
    pos <- st_coordinates(suppressWarnings(st_centroid(hm)))
    idx_true <- match(hex$region, mx$region)
    disp <- sqrt((pos[, 1] - tcn[idx_true, 1])^2 +
                 (pos[, 2] - tcn[idx_true, 2])^2) / 1000

    # cohesion: true intra-state adjacencies surviving as touching cells
    touch <- st_is_within_distance(hm, hm, dist = 1, sparse = TRUE)
    pres <- 0L; tot <- 0L
    for (a in seq_len(n)) {
      js <- adj[[idx_true[a]]]
      js <- js[mx$state_code[js] == s]
      js <- match(mx$region[js], hex$region)
      js <- js[!is.na(js) & js > a]
      for (b in js) {
        tot <- tot + 1L
        if (b %in% touch[[a]]) pres <- pres + 1L
      }
    }

    # absolute print legibility: state on its own page, and the same cells on
    # a national page (where the 223x between-state cell-size range, which
    # gamma cannot touch, is what actually decides readability)
    bb <- st_bbox(hm)
    span_km <- max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"]) / 1000
    mm_per_km <- PRINT_MM / span_km
    dia_mm <- 2 * sqrt(a_km2 * mm_per_km^2 / pi)
    dia_mm_nat <- 2 * sqrt(a_km2 * (PRINT_MM_NAT / NAT_SPAN_KM)^2 / pi)

    # what the raw gamma rule asked for, and what the engine's clamp allowed
    wt <- pmax(true_km2[idx_true], 1e-12)^g
    req <- wt / sum(wt) * n
    cl <- req; bind <- rep(FALSE, n)
    for (r in 1:3) {
      cl2 <- pmax(pmin(cl, 2.2), 0.45)
      bind <- cl2 != cl              # test BEFORE the rescale moves them off
      cl <- cl2 / sum(cl2) * n
    }
    qr <- quantile(req, c(.05, .95))
    qc <- quantile(cl, c(.05, .95))

    iq <- vapply(seq_len(n), function(i) {
      co <- st_coordinates(hm$geometry[i])
      co <- co[!is.na(co[, 1]), , drop = FALSE]
      if (nrow(co) < 3) return(NA_real_)
      x <- co[, 1]; y <- co[, 2]
      per <- sum(sqrt((x - c(x[-1], x[1]))^2 + (y - c(y[-1], y[1]))^2))
      4 * pi * (a_km2[i] * 1e6) / per^2
    }, 0)

    # ink bias: area-weighted (what the reader sees) vs unweighted mean
    vv <- vt[match(hex$region, vt$region), VARS, drop = FALSE]
    ink <- as.data.frame(lapply(VARS, function(v) {
      x <- vv[[v]]; ok <- is.finite(x)
      c(sum(a_km2[ok] * x[ok]), sum(a_km2[ok]), sum(x[ok]), sum(ok))
    }))
    names(ink) <- VARS
    ink <- setNames(
      as.list(unlist(lapply(VARS, function(v) ink[[v]]))),
      as.vector(outer(c("wsum", "asum", "vsum", "vn"), VARS,
                      paste, sep = "_")))

    q <- quantile(r_eq, c(.05, .95))
    data.frame(
      gamma = g, state_code = s, n = n,
      true_km2_per_cell = round(sum(true_km2[idx_true]) / n, 1),
      cohesion_pct = 100 * pres / max(tot, 1), pairs = tot, pairs_kept = pres,
      disp_km = mean(disp), disp_km_p95 = as.numeric(quantile(disp, .95)),
      area_p5 = as.numeric(q[1]), area_p95 = as.numeric(q[2]),
      area_min = min(r_eq), area_max = max(r_eq),
      leg_ratio = as.numeric(q[1] / q[2]),
      leg_minmax = min(r_eq) / max(r_eq),
      tiny_frac = mean(r_eq < TINY_FRAC),
      leg_print = mean(dia_mm >= READ_MM),
      dia_mm_p5 = as.numeric(quantile(dia_mm, .05)),
      leg_print_nat = mean(dia_mm_nat >= READ_MM),
      dia_mm_nat_p5 = as.numeric(quantile(dia_mm_nat, .05)),
      leg_ratio_req = as.numeric(qr[1] / qr[2]),
      leg_ratio_clamp = as.numeric(qc[1] / qc[2]),
      clamp_bind_frac = mean(bind),
      iq_mean = mean(iq, na.rm = TRUE),
      as.data.frame(ink))
  }

  rows <- list(); k <- 0L
  for (g in GAMMAS) for (s in PANEL) {
    r <- score_state(g, s)
    if (is.null(r)) next
    k <- k + 1L; rows[[k]] <- r
    cat(sprintf("  scored g=%.2f s=%s coh=%5.1f%% leg=%.3f disp=%5.1fkm\n",
                g, s, r$cohesion_pct, r$leg_ratio, r$disp_km))
  }
  res <- do.call(rbind, rows)
  meta <- do.call(rbind, lapply(GAMMAS, function(g)
    do.call(rbind, lapply(PANEL, function(s) {
      mp <- meta_path(g, s); if (file.exists(mp)) readRDS(mp) else NULL
    }))))
  res <- left_join(res, meta, by = c("gamma", "state_code"))
  st_nm <- df_mxmunicipio_2020 |>
    distinct(state_code, state_abbr)
  res <- left_join(res, st_nm, by = "state_code") |>
    mutate(converged = maxerr_pct <= 5) |>
    arrange(state_code, gamma)
  write.csv(res, "data/gamma_frontier.csv", row.names = FALSE)
  cat("wrote data/gamma_frontier.csv\n")
}

# ---- stage 3: aggregate, knee, plot --------------------------------------
if (STAGE %in% c("plot", "all")) {
  res <- read.csv("data/gamma_frontier.csv",
                  colClasses = c(state_code = "character"))

  # panel aggregate: cohesion pair-weighted (the honest weighting), the rest
  # cell-weighted
  agg <- res |>
    group_by(gamma) |>
    summarise(cohesion_pct = 100 * sum(pairs_kept) / sum(pairs),
              disp_km   = sum(disp_km * n) / sum(n),
              leg_ratio = sum(leg_ratio * n) / sum(n),
              tiny_frac = sum(tiny_frac * n) / sum(n),
              leg_print = sum(leg_print * n) / sum(n),
              leg_print_nat = sum(leg_print_nat * n) / sum(n),
              leg_req   = sum(leg_ratio_req * n) / sum(n),
              leg_clamp = sum(leg_ratio_clamp * n) / sum(n),
              clamped   = sum(clamp_bind_frac * n) / sum(n),
              iq_mean   = sum(iq_mean * n) / sum(n),
              n_states  = dplyr::n(), .groups = "drop") |>
    arrange(gamma)

  # Kneedle: normalise both axes over the swept range, order by legibility,
  # take the point furthest from the chord joining the extremes.
  knee_of <- function(leg, coh, gam) {
    ok <- is.finite(leg) & is.finite(coh)
    leg <- leg[ok]; coh <- coh[ok]; gam <- gam[ok]
    if (length(unique(leg)) < 3) return(NA_real_)
    if (diff(range(leg)) == 0 || diff(range(coh)) == 0) return(NA_real_)
    x <- (leg - min(leg)) / diff(range(leg))
    y <- (coh - min(coh)) / diff(range(coh))
    o <- order(x); x <- x[o]; y <- y[o]; gg <- gam[o]
    x1 <- x[1]; y1 <- y[1]; x2 <- x[length(x)]; y2 <- y[length(y)]
    d <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) /
      sqrt((y2 - y1)^2 + (x2 - x1)^2)
    if (!any(is.finite(d))) return(NA_real_)
    gg[which.max(d)]
  }
  knee_panel <- knee_of(agg$leg_ratio, agg$cohesion_pct, agg$gamma)
  knee_state <- res |>
    group_by(state_code, state_abbr) |>
    summarise(knee_gamma = knee_of(leg_ratio, cohesion_pct, gamma),
              .groups = "drop")

  # exchange rate walking DOWN gamma (the direction the project cares about):
  # cohesion points paid per 0.01 of legibility bought
  ex <- agg |> arrange(desc(gamma)) |>
    mutate(d_leg = leg_ratio - lag(leg_ratio),
           d_coh = cohesion_pct - lag(cohesion_pct),
           d_disp = disp_km - lag(disp_km),
           coh_per_leg01 = -d_coh / (d_leg / 0.01),
           km_per_leg01 = d_disp / (d_leg / 0.01),
           step = paste0(lag(sprintf("%.2f", gamma)), "->",
                         sprintf("%.2f", gamma)))

  cat("\n================ PANEL AGGREGATE ================\n")
  print(as.data.frame(agg |>
    mutate(across(c(cohesion_pct, disp_km), ~round(.x, 1)),
           across(c(leg_ratio, tiny_frac, leg_print, leg_print_nat, leg_req,
                    leg_clamp, clamped, iq_mean), ~round(.x, 3)))),
    row.names = FALSE)
  cat("\n===== EXCHANGE RATE, walking gamma DOWNWARD =====\n")
  cat("(cohesion points surrendered per 0.01 legibility gained)\n")
  print(as.data.frame(ex |> filter(!is.na(d_leg)) |>
    transmute(step, d_leg = round(d_leg, 3), d_coh = round(d_coh, 2),
              d_disp_km = round(d_disp, 2),
              coh_per_leg01 = round(coh_per_leg01, 3),
              km_per_leg01 = round(km_per_leg01, 3))), row.names = FALSE)
  cat("\n===== PER-STATE KNEE (Kneedle on legibility/cohesion) =====\n")
  print(as.data.frame(knee_state), row.names = FALSE)
  if (is.finite(knee_panel)) {
    cat(sprintf("\nPANEL KNEE: gamma = %.2f\n", knee_panel))
  } else {
    cat("\nPANEL KNEE: undefined (a frontier axis did not vary)\n")
  }
  # ---- the verdict the owner needs -------------------------------------
  at <- function(g, col) agg[[col]][which.min(abs(agg$gamma - g))]
  seg <- function(g_hi, g_lo) sprintf(
    paste("  gamma %.2f -> %.2f : legibility %+.3f | cohesion %+.1f pts |",
          "disp %+.1f km | %.3f cohesion pts per 0.01 legibility"),
    g_hi, g_lo, at(g_lo, "leg_ratio") - at(g_hi, "leg_ratio"),
    at(g_lo, "cohesion_pct") - at(g_hi, "cohesion_pct"),
    at(g_lo, "disp_km") - at(g_hi, "disp_km"),
    (at(g_hi, "cohesion_pct") - at(g_lo, "cohesion_pct")) /
      ((at(g_lo, "leg_ratio") - at(g_hi, "leg_ratio")) / 0.01))
  cat("\n===== VERDICT: three regimes, walking gamma DOWNWARD =====\n")
  cat("DEAD ZONE\n"); cat(seg(1.00, 0.60), "\n")
  cat("EXPENSIVE STRETCH\n"); cat(seg(0.60, 0.30), "\n")
  cat("CHEAP STRETCH\n"); cat(seg(0.30, 0.00), "\n")
  cat(sprintf(
    paste0("\nThe clamp, not gamma, sets the floor: requested p5/p95 falls to ",
           "%.3f at\ngamma = 1 but the engine's [0.45, 2.2] x share target ",
           "clamp only allows\n%.3f, and the layout delivers %.3f. Anything ",
           "above gamma ~= 0.60 is inert.\n"),
    at(1.0, "leg_req"), at(1.0, "leg_clamp"), at(1.0, "leg_ratio")))
  cat(sprintf(
    paste0("Legibility tail is empty at EVERY gamma (%d of %d rows with any ",
           "cell under\n%.2fx its state mean) -- the clamp makes an invisible ",
           "cell impossible here.\n"),
    sum(res$tiny_frac > 0), nrow(res), TINY_FRAC))
  cat(sprintf(
    paste0("On a national %d mm page only %.0f%% of panel cells clear %d mm ",
           "at ANY gamma:\nthe 223x between-state cell-size range is outside ",
           "gamma's reach.\n"),
    PRINT_MM_NAT, 100 * at(0.30, "leg_print_nat"), READ_MM))
  VARS <- c("pct_indig", "pct_afro", "sex_ratio", "log_dens")
  if (all(paste0("wsum_", VARS) %in% names(res))) {
    ink <- do.call(rbind, lapply(VARS, function(v) {
      d <- res |> group_by(gamma) |>
        summarise(visual = sum(.data[[paste0("wsum_", v)]]) /
                    sum(.data[[paste0("asum_", v)]]),
                  plain = sum(.data[[paste0("vsum_", v)]]) /
                    sum(.data[[paste0("vn_", v)]]), .groups = "drop")
      data.frame(variable = v, gamma = d$gamma,
                 ink_weighted = round(d$visual, 3),
                 per_municipio = round(d$plain, 3),
                 bias = round(d$visual - d$plain, 3))
    }))
    cat("\n===== READING REAL DATA: ink-weighted vs per-municipio mean =====\n")
    cat("(panel-wide: what the map's average LOOKS like, vs what it IS)\n")
    print(as.data.frame(ink), row.names = FALSE)

    # WITHIN-state bias is the part gamma actually controls: the panel-wide
    # figure above stays biased at gamma = 0 because equal cells are equal
    # inside a state, and the 223x between-state cell-size range survives.
    wb <- do.call(rbind, lapply(VARS, function(v) {
      d <- res
      vis <- d[[paste0("wsum_", v)]] / d[[paste0("asum_", v)]]
      pla <- d[[paste0("vsum_", v)]] / d[[paste0("vn_", v)]]
      agg2 <- data.frame(gamma = d$gamma, n = d$n, b = abs(vis - pla)) |>
        group_by(gamma) |>
        summarise(mab = sum(b * n) / sum(n), .groups = "drop")
      data.frame(variable = v, gamma = agg2$gamma,
                 within_state_abs_bias = round(agg2$mab, 4))
    }))
    cat("\n== the part gamma controls: mean |within-state ink bias| ==\n")
    print(as.data.frame(
      do.call(cbind, c(list(gamma = unique(wb$gamma)),
        setNames(lapply(VARS, function(v)
          wb$within_state_abs_bias[wb$variable == v]), VARS)))),
      row.names = FALSE)
  }

  cat("\n===== COST OF GOING FROM THE KNEE (0.30) TO FULLY EQUAL (0) =====\n")
  d0 <- res |> filter(gamma %in% c(0, 0.3)) |>
    select(state_code, state_abbr, n, gamma, cohesion_pct, disp_km,
           leg_ratio) |>
    arrange(state_code, gamma) |>
    group_by(state_code, state_abbr, n) |>
    summarise(d_cohesion = cohesion_pct[gamma == 0] -
                cohesion_pct[gamma == 0.3],
              d_disp_km = disp_km[gamma == 0] - disp_km[gamma == 0.3],
              d_legib = leg_ratio[gamma == 0] - leg_ratio[gamma == 0.3],
              .groups = "drop")
  print(as.data.frame(d0 |> mutate(across(where(is.numeric), ~round(.x, 2)))),
        row.names = FALSE)

  nonconv <- sum(!res$converged)
  cat(sprintf("Rows failing the 5%% area-convergence gate: %d / %d\n",
              nonconv, nrow(res)))
  if (nonconv > 0)
    print(as.data.frame(res[!res$converged,
      c("gamma", "state_code", "knn_used", "maxerr_pct")]), row.names = FALSE)

  # ---- plot -------------------------------------------------------------
  res$lab <- sprintf("%02d", as.integer(res$state_code))
  res$facet <- sprintf("%s %s (n=%d)", res$lab, res$state_abbr, res$n)
  kn <- left_join(knee_state, res, by = c("state_code", "state_abbr")) |>
    filter(is.finite(knee_gamma) & abs(gamma - knee_gamma) < 1e-9)

  p1 <- ggplot(res, aes(leg_ratio, cohesion_pct,
                        colour = state_abbr, group = state_abbr)) +
    geom_path(linewidth = .6, alpha = .85) +
    geom_point(aes(size = gamma), alpha = .85) +
    geom_point(data = kn, shape = 21, size = 6, stroke = 1.1,
               fill = NA, colour = "black") +
    geom_text(data = res |> filter(gamma %in% c(0, 1)),
              aes(label = sprintf("g=%.1f", gamma)), size = 2.6,
              vjust = -1.1, show.legend = FALSE) +
    scale_size_continuous(range = c(1.2, 3.6)) +
    scale_x_continuous(labels = label_number(accuracy = .01)) +
    labs(title = "Pareto frontier: legibility vs adjacency, per state",
         subtitle = paste("x = area p5/p95 within state (1 = equal cells).",
                          "Ringed point = Kneedle knee."),
         x = "legibility  (area p5 / p95)", y = "cohesion (%)",
         colour = "state", size = "gamma") +
    theme_minimal(base_size = 10)

  p2 <- ggplot(agg, aes(leg_ratio, cohesion_pct)) +
    geom_path(aes(x = leg_req), linewidth = .5, linetype = 2,
              colour = "grey65") +
    geom_path(linewidth = .9, colour = "grey30") +
    geom_point(aes(fill = gamma), shape = 21, size = 4, colour = "grey20") +
    geom_text(aes(label = sprintf("%.2f", gamma)), size = 2.7, vjust = -1.3) +
    geom_point(data = agg |>
                 filter(is.finite(knee_panel) & abs(gamma - knee_panel) < 1e-9),
               shape = 21, size = 8, stroke = 1.2, fill = NA,
               colour = "firebrick") +
    scale_fill_viridis_c(option = "C") +
    labs(title = sprintf("Panel aggregate frontier - KNEE at gamma = %.2f",
                         knee_panel),
         subtitle = paste("solid = delivered, dashed = what gamma asked for",
                          "before the engine's area clamp"),
         x = "legibility  (area p5 / p95, cell-weighted)",
         y = "cohesion (%), pair-weighted") +
    theme_minimal(base_size = 10)

  long <- do.call(rbind, lapply(
    c("legibility", "cohesion_pct_scaled", "disp_km_scaled", "tiny_frac"),
    function(v) {
      d <- res
      val <- switch(v,
        legibility = d$leg_ratio,
        cohesion_pct_scaled = d$cohesion_pct / 100,
        disp_km_scaled = d$disp_km / max(res$disp_km),
        tiny_frac = d$tiny_frac)
      data.frame(gamma = d$gamma, facet = d$facet, metric = v, value = val)
    }))
  p3 <- ggplot(long, aes(gamma, value, colour = metric)) +
    geom_line(linewidth = .7) + geom_point(size = 1.1) +
    facet_wrap(~facet, nrow = 2) +
    scale_colour_manual(values = c(legibility = "#1b7837",
                                   cohesion_pct_scaled = "#2166ac",
                                   disp_km_scaled = "#b2182b",
                                   tiny_frac = "#f1a340")) +
    labs(title = "Each metric against gamma, per state (all scaled to 0-1)",
         subtitle = paste("disp_km divided by the panel max",
                          round(max(res$disp_km), 1), "km"),
         x = "gamma", y = "scaled value", colour = NULL) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom")

  pp <- (p1 | p2) / p3 + plot_layout(heights = c(1, 1.05))
  ggsave("state_pages/gamma_frontier.png", pp, width = 15, height = 10,
         dpi = 150, bg = "white")
  cat("wrote state_pages/gamma_frontier.png\n")
}
