# ---------------------------------------------------------------------------
# Purpose: check two claims about the canonical artifact before either is
#          written into the project's documentation.
#
#          CLAIM 1. "v11.2's cells are all convex." CLAUDE.md and several
#          commit messages say so, and the reasoning behind rejecting non-convex
#          cell families rests on it. A power diagram's cells ARE convex, but
#          this artifact clips them to the real state outline, which is not.
#
#          CLAIM 2. "56-64% of municipios' cells do not intersect their own true
#          territory at all." If true, mean displacement in km understates how
#          far the layout moves things and every per-state km figure is
#          misleading. Striking enough that it must not be repeated unchecked.
#
# Author:  Claude (Opus 5) for Jorge Roa
# ISS:     -
# Version: verify v1
# Date:    27-August-2026
#
# Usage:   Rscript verify_v11_2_claims.R [rds_path]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(mxmaps)
})
sf_use_s2(FALSE)
setwd("/Users/jorgeroa/Documents/GitHub/mxmaps")
source("hexlayout_common.R")

MX_CRS <- 6372

args <- commandArgs(trailingOnly = TRUE)
rds <- if (length(args) >= 1) args[1] else "data/mxmunicipio_hex_sf_v11_2.rds"
hex <- readRDS(rds)
cat("Artifact:", rds, "\n\n")

hex_m <- st_transform(hex, MX_CRS)
mx <- build_mxmun_sf()
mx_m <- st_transform(mx, MX_CRS)

# ---- CLAIM 1: convexity --------------------------------------------------
a_cell <- as.numeric(st_area(hex_m))
a_hull <- as.numeric(st_area(st_convex_hull(hex_m)))
solidity <- a_cell / a_hull
# 1e-9 relative tolerance: a clipped cell whose outline is a straight state
# border is convex to within floating point, and should count as convex.
is_convex <- solidity > 1 - 1e-9

cat("== CLAIM 1: are v11.2's cells convex? ==\n")
cat(sprintf("strictly convex cells: %d / %d = %.1f%%\n",
            sum(is_convex), length(is_convex), 100 * mean(is_convex)))
cat(sprintf("mean solidity: %.4f | min %.4f | p5 %.4f\n",
            mean(solidity), min(solidity), quantile(solidity, 0.05)))
cat(sprintf("cells with solidity < 0.90: %d (%.1f%%)\n",
            sum(solidity < 0.90), 100 * mean(solidity < 0.90)))
cat(sprintf("multi-part cells: %d\n",
            sum(vapply(st_geometry(hex_m),
                       function(g) length(st_cast(st_sfc(g), "POLYGON")) > 1,
                       TRUE))))

# ---- CLAIM 2: does a cell cover any of its own municipio? ----------------
same <- match(hex$region, mx$region)
stopifnot(!any(is.na(same)))
hits <- st_intersects(hex_m, mx_m[same, ], sparse = FALSE)
own <- diag(hits)

# Overlap AREA, not just a boolean: a cell can clip a corner of its own
# territory and still be essentially elsewhere.
ov_frac <- vapply(seq_len(nrow(hex_m)), function(i) {
  g <- suppressWarnings(tryCatch(
    st_intersection(hex_m$geometry[i], mx_m$geometry[same[i]]),
    error = function(e) NULL))
  if (is.null(g) || length(g) == 0) return(0)
  p <- suppressWarnings(tryCatch(st_collection_extract(g, "POLYGON"),
                                 error = function(e) NULL))
  if (is.null(p) || length(p) == 0) return(0)
  sum(as.numeric(st_area(p))) / a_cell[i]
}, 0)

cat("\n== CLAIM 2: does each cell intersect its own true municipio? ==\n")
cat(sprintf("cells intersecting own territory at all: %d / %d = %.1f%%\n",
            sum(own), length(own), 100 * mean(own)))
cat(sprintf("cells NOT intersecting own territory:    %d = %.1f%%\n",
            sum(!own), 100 * mean(!own)))
cat(sprintf("overlap as a fraction of cell area: mean %.3f median %.3f\n",
            mean(ov_frac), median(ov_frac)))
cat(sprintf("cells with < 1%% of their area on own territory: %d (%.1f%%)\n",
            sum(ov_frac < 0.01), 100 * mean(ov_frac < 0.01)))
cat(sprintf("cells with < 50%% of their area on own territory: %d (%.1f%%)\n",
            sum(ov_frac < 0.50), 100 * mean(ov_frac < 0.50)))
