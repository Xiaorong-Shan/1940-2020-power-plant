#!/usr/bin/env Rscript
# ============================================================
# PM2.5 (annual) for 1940/50/60/70/80/90 using 1999 met proxy
# Minimal + robust:
#   1) Pre-aggregate HyADS to enforce uniqueness:
#      - by-unit: sum(hyads) by (x, y, uID)
#      - total : sum(hyads) by (x, y)
#      (Prevents data.table::dcast defaulting to `length`, which yields near-zeros)
#   2) Always use mask.use = mask.usa (CONUS, exclude AK/HI)
#   3) “Fake 1999” inputs -> run with 1999 met -> rename outputs back to target year
#   4) Build annual TOTAL from the by-unit output (summing unit columns), robust to schema
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
  library(sf)
  library(USAboundaries)
})

# ----------------------------- CONFIG ---------------------------------
YEARS          <- c(1940, 1950, 1960, 1970, 1980, 1990)  # decades to run (annual)
MET_PROXY_YEAR <- 1999L                                  # force meteorology year to 1999
NAME_X         <- "hyads"                                # exposure column expected by the model function

# Input/Output paths (adjust if needed)
EXPO_DIR   <- "/scratch/xshan2/R_Code/disperseR/main/output/exposures_proxy1999_shifted_1940_1990"
OUT_DIR    <- "/scratch/xshan2/R_Code/disperseR/main/output/pm25_decades_model.lm.cv_single_poly_proxy1999met"

# Model & function
FX_PATH    <- "/scratch/xshan2/lucas/hyads_to_pm25_functions_update.R"
MODEL_RDA  <- "/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/hyads_to_cmaq_models20230516.RData"
MODEL_OBJ  <- "preds.ann.hyads06w05"                     # object name inside MODEL_RDA
MODEL_NAME <- "model.lm.cv_single_poly"                  # model selector used by the function
MET_DEST   <- "/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/met"

# Grid CRS (must match your HyADS grid; here: CONUS Albers)
P4S <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +ellps=GRS80 +datum=NAD83 +units=m"

# Ensure output directory exists
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------- LOAD FUNCTIONS & MODELS ----------------------
# Use UTC to avoid met time alignment issues (e.g., PBL dates)
Sys.setenv(TZ = "UTC")

# Source conversion function (provides hyads_to_pm25_unit)
stopifnot(file.exists(FX_PATH), file.exists(MODEL_RDA))
source(FX_PATH)

# Load annual model objects and select the one we want
load(MODEL_RDA)
stopifnot(exists(MODEL_OBJ, inherits = TRUE))
model.dataset <- get(MODEL_OBJ, inherits = TRUE)

# ----------------------------- MASK -----------------------------------
# Build a CONUS mask (exclude AK/HI) in the same projection as HyADS
us_states <- st_transform(USAboundaries::us_states(), P4S)
mask.usa  <- sf::as_Spatial(us_states)[us_states$state_abbr %in% setdiff(state.abb, c("HI","AK")),]

# --------------------------- UTILITIES --------------------------------
# Safe/soft reader: return NULL if file is absent or unreadable
read_fst_or_null <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  tryCatch(read_fst(fp, as.data.table = TRUE), error = function(e) NULL)
}

# Copy a 1999-named file and rename it to the target year (preserve directory)
copy_as_year <- function(fp1999, yr, out_dir = OUT_DIR) {
  stopifnot(file.exists(fp1999))
  bn  <- basename(fp1999)
  ext <- tools::file_ext(bn)
  # Prefer replacing "_1999." with "_<yr>."
  bn2 <- sub(paste0("_", MET_PROXY_YEAR, "\\."), paste0("_", yr, "."), bn)
  # Fallback: append _<yr> if pattern not found
  if (identical(bn2, bn)) bn2 <- paste0(tools::file_path_sans_ext(bn), "_", yr, ".", ext)
  out <- file.path(out_dir, bn2)
  file.copy(fp1999, out, overwrite = TRUE)
  out
}

# --------------------------- CORE RUNNER ------------------------------
# Run one target year (annual): fake inputs as 1999 -> run -> rename
run_year <- function(YR) {
  message("== YEAR ", YR, " | MET proxy = ", MET_PROXY_YEAR, " ==")

  # ---- 1) Read HyADS inputs for the target year ----
  f_byu_in <- file.path(EXPO_DIR, sprintf("grids_exposures_byunit_%d.fst", YR))
  f_tot_in <- file.path(EXPO_DIR, sprintf("grids_exposures_total_%d.fst",  YR))

  byu <- read_fst_or_null(f_byu_in)
  tot <- read_fst_or_null(f_tot_in)
  if (is.null(byu) || is.null(tot)) {
    message("  ✖ Missing HyADS inputs for ", YR, ":\n    ", f_byu_in, "\n    ", f_tot_in)
    return(invisible(NULL))
  }

  # ---- 2) PRE-AGGREGATE HyADS to ensure uniqueness (prevents dcast 'length') ----
  # by-unit table needs columns: x, y, uID, hyads
  if (!(NAME_X %in% names(byu))) {
    # Try to identify the HyADS column automatically if it's not named 'hyads'
    cand <- intersect(names(byu), c("hyads","hyads_weighted","N","counts","count","value","val"))
    stopifnot(length(cand) >= 1)
    setnames(byu, cand[1], NAME_X)
  }
  stopifnot(all(c("x","y","uID", NAME_X) %in% names(byu)))
  byu <- byu[, .(hyads = sum(get(NAME_X), na.rm = TRUE)), by = .(x, y, uID)]

  # total table needs columns: x, y, hyads
  if (!(NAME_X %in% names(tot))) {
    cand <- intersect(names(tot), c("hyads","hyads_weighted","N","counts","count","value","val"))
    stopifnot(length(cand) >= 1)
    setnames(tot, cand[1], NAME_X)
  }
  stopifnot(all(c("x","y", NAME_X) %in% names(tot)))
  tot <- tot[, .(hyads = sum(get(NAME_X), na.rm = TRUE)), by = .(x, y)]

  # ---- 3) Write FAKE 1999 inputs into a temporary working folder ----
  tmp_dir <- file.path(OUT_DIR, sprintf(".fake_%d_as_%d", YR, MET_PROXY_YEAR))
  run_dir <- file.path(tmp_dir, "out")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

  byu[, year := MET_PROXY_YEAR][]
  tot[, year := MET_PROXY_YEAR][]
  fst::write_fst(byu, file.path(tmp_dir, sprintf("grids_exposures_byunit_%d.fst", MET_PROXY_YEAR)))
  fst::write_fst(tot, file.path(tmp_dir, sprintf("grids_exposures_total_%d.fst",  MET_PROXY_YEAR)))

  # ---- 4) Run annual PM2.5 (by-unit path), using 1999 met + CONUS mask ----
  hyads_to_pm25_unit(
    year.m        = MET_PROXY_YEAR,                                # forces 1999 meteorology usage
    month.n       = NULL,
    fstart        = file.path(tmp_dir, "grids_exposures_byunit_"), # points to our fake-1999 by-unit HyADS
    fstart.total  = file.path(tmp_dir, "grids_exposures_total_"),  # supplied for API compatibility
    fstart_out    = file.path(run_dir, "grids_pm25_byunit_"),      # by-unit PM out (1999 naming)
    model.dataset = model.dataset,                                 # preds.ann.hyads06w05
    model.name    = MODEL_NAME,                                    # "model.lm.cv_single_poly"
    name.x        = NAME_X,                                        # "hyads"
    met.dest      = MET_DEST,                                      # meteorology root
    mask.use      = mask.usa,                                      # ★ required: apply CONUS mask
    p4s           = P4S,
    total         = FALSE                                          # by-unit path (not total)
  )

  # ---- 5) Locate the 1999 by-unit PM output (fst preferred; rds fallback) ----
  out1999 <- file.path(run_dir, sprintf("grids_pm25_byunit_%d.fst", MET_PROXY_YEAR))
  if (!file.exists(out1999)) {
    alt <- file.path(run_dir, sprintf("grids_pm25_byunit_%d.rds", MET_PROXY_YEAR))
    if (file.exists(alt)) out1999 <- alt
  }
  stopifnot(file.exists(out1999))

  # ---- 6) Copy & rename the 1999 by-unit output to the target year ----
  out_byu <- copy_as_year(out1999, YR, OUT_DIR)
  message("  ✓ BY-UNIT PM written: ", out_byu)

  # ---- 7) Construct annual TOTAL from the by-unit output (robust to schema) ----
  #    Wide schema (many numeric unit columns): total = rowSums(unit columns)
  #    Already-a-total schema (has 'vals.out'): just keep (x,y,vals.out)
  #    Long schema (x,y,uID,value): sum by (x,y)
  res_byu <- if (grepl("\\.rds$", out_byu, ignore.case = TRUE)) {
    x <- readRDS(out_byu); if (!is.data.table(x)) as.data.table(x) else x
  } else {
    fst::read_fst(out_byu, as.data.table = TRUE)
  }

  cols     <- names(res_byu)
  num_cols <- cols[vapply(res_byu, is.numeric, logical(1))]
  drop_nc  <- intersect(cols, c("x","y","year","hyads","E","emiss","so2","SO2","vals.out"))
  unit_nc  <- setdiff(num_cols, drop_nc)

  if (length(unit_nc) >= 1) {
    # WIDE: sum across unit columns
    res_tot <- res_byu[, .(x, y)]
    res_tot[, vals.out := rowSums(res_byu[, ..unit_nc], na.rm = TRUE)]
  } else if ("vals.out" %in% cols) {
    # ALREADY TOTAL: keep as-is
    res_tot <- res_byu[, .(x, y, vals.out)]
  } else if ("uID" %in% cols) {
    # LONG: sum a numeric value column by (x, y)
    vcol <- setdiff(num_cols, c("x","y"))[1]
    res_tot <- res_byu[, .(vals.out = sum(get(vcol), na.rm = TRUE)), by = .(x, y)]
  } else {
    stop("Cannot derive total from by-unit output. Please inspect columns of: ", out_byu)
  }

  # ---- 8) Write annual TOTAL to disk and print quick QC ----
  out_tot <- file.path(OUT_DIR, sprintf("grids_pm25_total_%d.fst", YR))
  fst::write_fst(res_tot, out_tot)
  message("  ✓ TOTAL PM written:  ", out_tot)

  if (nrow(res_tot)) {
    rng <- range(res_tot$vals.out, na.rm = TRUE)
    q99 <- quantile(res_tot$vals.out, c(.9,.99,.999), na.rm = TRUE)
    message(sprintf("  QC: range[%.4g, %.4g], q90=%.4g q99=%.4g q99.9=%.4g",
                    rng[1], rng[2], q99[1], q99[2], q99[3]))
  } else {
    message("  ⚠ TOTAL has 0 rows (empty).")
  }
}

# ----------------------------- RUN ALL --------------------------------
for (y in YEARS) {
  try(run_year(y), silent = FALSE)
}
message("Done: decades ", paste(YEARS, collapse = ", "),
        " (all run with meteorology proxy year = ", MET_PROXY_YEAR, ").")
