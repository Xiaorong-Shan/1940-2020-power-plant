###############################################################
#  run_sensitivity_spatial_thresholds_hyads_mapping.R
#
#  PURPOSE
#  -------
#  Sensitivity analysis for the "spatial boundary" (maximum
#  nearest-neighbor distance, MAX_KM) used when mapping power
#  plants (from the historical emissions CSV) to HyADS source
#  units (uIDs) derived from the 2022 per-ton kernel.
#
#  WHAT IT TESTS (Sensitivity)
#  ---------------------------
#  • How results change when we tighten or relax the spatial
#    mapping boundary: MAX_KM ∈ {Inf, 300, 200, 100} km by default.
#  • Impact on: kept mappings, distance distribution, retained
#    SO2 emissions share, and number of usable years.
#
#  MAIN OUTPUTS (per MAX_KM)
#  -------------------------
#  /.../weighted_fst_MAXKM<k>/
#    ├─ grids_exposures_byunit_<YYYY>.fst      (annual, weighted)
#    ├─ grids_exposures_total_<YYYY>.fst       (annual, weighted)
#    ├─ plant_year_hyads_at_plant.csv          (per-plant QA, includes weighted hyads_at_plant_w)
#    └─ uID_mapping_QA_full.csv                (all mappings with distance and far-flag)
#
#  MASTER SUMMARY
#  --------------
#  /.../sensitivity_summary.csv                (side-by-side metrics for all thresholds)
#
#  ASSUMPTIONS
#  -----------
#  • HyADS kernel (2022) is a per-ton, by-unit grid with columns: x, y, uID, hyads
#  • Historical emissions CSV includes: Year, Plant_ID, Generator_ID,
#    convert_fuel, Longitude, Latitude, SO2_tons (or *_kg)
#  • Coal-only analysis by default (edit filter if you want Oil/Gas).
#
#  AUTHOR / DATE
#  -------------
#  [Your Name], [Date]
###############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
  library(sf)
  library(fs)
})

options(stringsAsFactors = FALSE)

# =================== USER CONFIG =================== #
# Per-ton HyADS kernel (2022)
HYADS_KERNEL_FST <- "/home/xshan2/HAQ_LAB/lhennem/disperseR/main/output/exp/coal_hyads_no_weight/grids_exposures_byunit_2022.fst"

# Historical emissions (1940–1990)
EMISS_CSV        <- "/scratch/xshan2/R_Code/powerplant/data/so2_facility_emissions_1940_1990.csv"

# Base output directory (each threshold creates a subfolder)
BASE_OUT_DIR     <- "/scratch/xshan2/R_Code/disperseR"

# Years to export weighted FSTs for
YEARS_ANN        <- 1940:1990

# Spatial reference and I/O settings
P4S_ALBERS       <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
FST_COMPRESS     <- 50

# Sensitivity set: spatial boundary thresholds (km)
MAX_KM_LIST      <- c(Inf, 300, 200, 100)
# =================================================== #


cat("\n[INIT] Load HyADS kernel (2022) ...\n")
kernel <- read_fst(HYADS_KERNEL_FST, as.data.table = TRUE)
stopifnot(all(c("x","y","uID","hyads") %in% names(kernel)))
kernel <- kernel[is.finite(x) & is.finite(y) & is.finite(hyads)]
kernel[, uID := as.character(uID)]
setkey(kernel, uID)
setorder(kernel, uID, -hyads)

# uID peaks (Albers meters) for nearest-neighbor mapping
uID_peak <- kernel[, .SD[1L, .(uID, src_x_m = x, src_y_m = y, hyads_peak = hyads)], by = uID]
uID_sf   <- st_as_sf(uID_peak, coords = c("src_x_m","src_y_m"), crs = P4S_ALBERS, remove = FALSE)
uID_m    <- st_transform(uID_sf, 3857)  # meters, for NN search

# Split kernel by uID for fast weighting
kernel_small <- kernel[, .(x, y, uID, hyads)]
kernel_small[, uID := as.character(uID)]
kernel_split <- split(kernel_small, f = kernel_small$uID, drop = TRUE)
names(kernel_split) <- as.character(unique(kernel_small$uID))

cat("Total uIDs in kernel:", length(kernel_split), "\n")

# =================== Emissions (Coal only) =================== #
cat("\n[INIT] Load emissions CSV (Coal only) ...\n")
emiss <- fread(EMISS_CSV)

# Normalize column names
if ("year" %in% names(emiss))         setnames(emiss, "year", "Year")
if ("plant_id" %in% names(emiss))     setnames(emiss, "plant_id", "Plant_ID")
if ("generator_id" %in% names(emiss)) setnames(emiss, "generator_id", "Generator_ID")
if ("convert_f" %in% names(emiss))    setnames(emiss, "convert_f", "convert_fuel")

stopifnot(all(c("Year","Plant_ID","Generator_ID","convert_fuel","Longitude","Latitude") %in% names(emiss)))

# Fuel filter — Coal only (extend to Oil/Gas if needed)
emiss <- emiss[convert_fuel %chin% c("Coal","COAL","coal")]

# Coordinates and units
emiss[, plant_lon := suppressWarnings(as.numeric(Longitude))]
emiss[, plant_lat := suppressWarnings(as.numeric(Latitude))]
emiss <- emiss[is.finite(plant_lon) & is.finite(plant_lat)]

if (!"SO2_tons" %in% names(emiss)) {
  if ("so2_kg" %in% names(emiss)) emiss[, SO2_tons := so2_kg / 907.185]
  if ("SO2_kg" %in% names(emiss)) emiss[, SO2_tons := SO2_kg / 907.185]
}
stopifnot("SO2_tons" %in% names(emiss))

setkey(emiss, Year, Plant_ID)
years_use <- sort(unique(emiss$Year))
cat("Years present in the emissions CSV:", paste(head(years_use, 8), collapse = ","), " ...\n")


# =================== Worker: run one threshold =================== #
run_one_threshold <- function(MAX_KM) {
  tag <- if (is.infinite(MAX_KM)) "Inf" else as.character(MAX_KM)
  out_dir <- file.path(BASE_OUT_DIR, sprintf("weighted_fst_MAXKM%s", tag))
  dir_create(out_dir, recurse = TRUE)

  cat("\n================= RUN MAX_KM =", tag, "=================\n")

  # B) Map (Year, Plant_ID) → nearest uID (no empty placeholders)
  map_year_plant_uid <- rbindlist(lapply(years_use, function(yr) {
    sub <- unique(emiss[Year == yr, .(Plant_ID, plant_lon, plant_lat)])
    if (nrow(sub) == 0L) return(data.table())
    pl_sf <- st_as_sf(sub, coords = c("plant_lon","plant_lat"), crs = 4326, remove = FALSE)
    pl_m  <- st_transform(pl_sf, 3857)
    idx   <- st_nearest_feature(pl_m, uID_m)
    d_km  <- as.numeric(st_distance(pl_m, uID_m[idx, ], by_element = TRUE)) / 1000
    data.table(Year = yr, Plant_ID = sub$Plant_ID, uID = uID_peak$uID[idx], nn_km = d_km)
  }), fill = TRUE)

  map_year_plant_uid[, uID := as.character(uID)]
  setkey(map_year_plant_uid, Year, Plant_ID)

  # QA — all mappings (before filtering)
  qa_full <- copy(map_year_plant_uid)
  qa_full[, flag_far := if (is.infinite(MAX_KM)) FALSE else nn_km > MAX_KM]
  fwrite(qa_full, file.path(out_dir, "uID_mapping_QA_full.csv"))

  # Distance filtering (the “spatial boundary”)
  map_year_plant_uid_filt <- if (is.infinite(MAX_KM)) {
    map_year_plant_uid[!is.na(uID)]
  } else {
    map_year_plant_uid[!is.na(uID) & nn_km <= MAX_KM]
  }

  # C) Aggregate emissions to (Year, uID)
  emiss_plant_year <- emiss[, .(SO2_tons = sum(SO2_tons, na.rm = TRUE)), by = .(Year, Plant_ID)]
  emiss_uid_year   <- map_year_plant_uid_filt[emiss_plant_year, on = .(Year, Plant_ID), nomatch = 0L][
    , .(SO2_tons = sum(SO2_tons, na.rm = TRUE)), by = .(Year, uID)]
  emiss_uid_year[, uID := as.character(uID)]
  setkey(emiss_uid_year, Year, uID)

  write_fst(emiss_uid_year, file.path(out_dir, "emiss_uid_year.fst"), compress = FST_COMPRESS)

  # D1) Compute hyads_at_plant and its SO2-weighted version
  # Reattach plant coordinates and project to Albers
  plant_xy <- unique(emiss[, .(Year, Plant_ID, plant_lon, plant_lat)])
  map_year_plant_uid_filt <- map_year_plant_uid_filt[plant_xy, on = .(Year, Plant_ID), nomatch = 0L]
  map_year_plant_uid_filt <- map_year_plant_uid_filt[is.finite(plant_lon) & is.finite(plant_lat)]
  pl_sf  <- st_as_sf(map_year_plant_uid_filt, coords = c("plant_lon","plant_lat"), crs = 4326, remove = FALSE)
  pl_alb <- st_transform(pl_sf, P4S_ALBERS)
  pl_xy  <- st_coordinates(pl_alb)
  map_year_plant_uid_filt[, `:=`(plant_x = pl_xy[,1], plant_y = pl_xy[,2],
                                 grid_x = NA_real_, grid_y = NA_real_, hyads_at_plant = NA_real_)]

  # Nearest grid cell within each uID kernel
  valid_coord <- is.finite(map_year_plant_uid_filt$plant_x) & is.finite(map_year_plant_uid_filt$plant_y)
  for (uid in unique(map_year_plant_uid_filt$uID)) {
    uid_chr <- as.character(uid)
    ki <- kernel_split[[uid_chr]]
    if (is.null(ki) || !nrow(ki)) next
    ki <- ki[is.finite(x) & is.finite(y) & is.finite(hyads)]
    if (!nrow(ki)) next
    idx <- which(map_year_plant_uid_filt$uID == uid_chr & valid_coord)
    if (!length(idx)) next
    for (i in idx) {
      px <- map_year_plant_uid_filt$plant_x[i]
      py <- map_year_plant_uid_filt$plant_y[i]
      d2 <- (ki$x - px)^2 + (ki$y - py)^2
      if (!any(is.finite(d2))) next
      d2[!is.finite(d2)] <- Inf
      j <- which.min(d2)
      map_year_plant_uid_filt$grid_x[i]         <- ki$x[j]
      map_year_plant_uid_filt$grid_y[i]         <- ki$y[j]
      map_year_plant_uid_filt$hyads_at_plant[i] <- ki$hyads[j]
    }
  }

  plant_map_emiss <- map_year_plant_uid_filt[
    emiss_plant_year, on = .(Year, Plant_ID), nomatch = 0L]
  plant_map_emiss[, hyads_at_plant_w := hyads_at_plant * SO2_tons]

  fwrite(plant_map_emiss[, .(Year, Plant_ID, uID, nn_km,
                             plant_lon, plant_lat, plant_x, plant_y,
                             grid_x, grid_y, SO2_tons,
                             hyads_at_plant, hyads_at_plant_w)],
         file.path(out_dir, "plant_year_hyads_at_plant.csv"))

  # D2) Write annual weighted FSTs (by-unit / total)
  BYUNIT_PREFIX <- "grids_exposures_byunit_"
  TOTAL_PREFIX  <- "grids_exposures_total_"

  for (yy in YEARS_ANN) {
    eyy <- emiss_uid_year[Year == yy]
    if (nrow(eyy) == 0L) {
      message("MAX_KM=", tag, " | Year ", yy, ": no valid plants → skip.")
      next
    }
    parts <- vector("list", nrow(eyy))
    for (i in seq_len(nrow(eyy))) {
      uid_i <- as.character(eyy$uID[i]); so2_i <- eyy$SO2_tons[i]
      ki <- kernel_split[[uid_i]]
      if (is.null(ki)) next
      ki_w <- copy(ki)
      ki_w[, `:=`(hyads = hyads * so2_i, uID = uid_i)]
      parts[[i]] <- ki_w
    }
    byunit_out <- rbindlist(parts, use.names = TRUE, fill = TRUE)
    if (nrow(byunit_out) == 0L) next
    byunit_out[, year := yy]
    data.table::setcolorder(byunit_out, c("x","y","uID","year","hyads"))
    fst::write_fst(byunit_out,
                   file.path(out_dir, sprintf("%s%d.fst", BYUNIT_PREFIX, yy)),
                   compress = FST_COMPRESS)

    total_out <- byunit_out[, .(hyads = sum(hyads, na.rm = TRUE)), by = .(x, y)]
    total_out[, year := yy]
    data.table::setcolorder(total_out, c("x","y","year","hyads"))
    fst::write_fst(total_out,
                   file.path(out_dir, sprintf("%s%d.fst", TOTAL_PREFIX, yy)),
                   compress = FST_COMPRESS)

    message("MAX_KM=", tag, " | Year ", yy, " written.")
  }

  # ===== Metrics for the master summary =====
  n_map_all   <- nrow(map_year_plant_uid)
  n_map_keep  <- nrow(map_year_plant_uid_filt)

  nn_stats <- map_year_plant_uid_filt[, .(
    n = .N,
    nn_p50 = quantile(nn_km, 0.5, na.rm = TRUE),
    nn_p90 = quantile(nn_km, 0.9, na.rm = TRUE),
    nn_p95 = quantile(nn_km, 0.95, na.rm = TRUE),
    nn_max = max(nn_km, na.rm = TRUE)
  )]

  total_so2_all  <- emiss[, sum(SO2_tons, na.rm = TRUE)]
  total_so2_keep <- plant_map_emiss[, sum(SO2_tons, na.rm = TRUE)]
  years_with_any <- length(unique(plant_map_emiss$Year))

  data.table(
    MAX_KM = tag,
    mapped_all = n_map_all,
    mapped_kept = n_map_keep,
    mapped_kept_ratio = ifelse(n_map_all == 0, NA_real_, n_map_keep / n_map_all),
    nn_p50 = nn_stats$nn_p50,
    nn_p90 = nn_stats$nn_p90,
    nn_p95 = nn_stats$nn_p95,
    nn_max = nn_stats$nn_max,
    so2_total_all = total_so2_all,
    so2_total_kept = total_so2_keep,
    so2_keep_ratio = ifelse(total_so2_all == 0, NA_real_, total_so2_keep / total_so2_all),
    years_with_any = years_with_any,
    out_dir = out_dir
  )
}

# =================== Run all thresholds =================== #
summary_list <- lapply(MAX_KM_LIST, run_one_threshold)
sens_summary <- rbindlist(summary_list, fill = TRUE)
fwrite(sens_summary, file.path(BASE_OUT_DIR, "sensitivity_summary.csv"))

cat("\n[DONE] Sensitivity finished.\n")
cat("  Master summary:\n    ", file.path(BASE_OUT_DIR, "sensitivity_summary.csv"), "\n")
cat("  Per-threshold outputs under:\n    ", file.path(BASE_OUT_DIR, "weighted_fst_MAXKM*"), "\n")
