###############################################################
#  FULL PIPELINE: Steps A–F  (Weighted HyADS Generation)
#
#  A. Build HyADS uID source points from 2022 per-ton kernel
#  B. Map (Year, Plant_ID) → nearest HyADS uID using your emission CSV
#  C. Aggregate emissions to (Year, uID)
#  D. Compute weighted HyADS at plant (hyads_at_plant_w)
#  E2. Compute grid-aligned translation offsets (plume alignment)
#  F. Generate decadal per-ton by-unit & total FSTs (1940–1990)
#
#  Author: Xiaorong Shan
#  Updated: October 15, 2025
###############################################################


suppressPackageStartupMessages({
  library(data.table)
  library(fst)
  library(sf)
  library(fs)
})

options(stringsAsFactors = FALSE)

# ===================== USER CONFIG ===================== #
HYADS_KERNEL_FST <- "/home/xshan2/HAQ_LAB/lhennem/disperseR/main/output/exp/coal_hyads_no_weight/grids_exposures_byunit_2022.fst"
EMISS_CSV        <- "/scratch/xshan2/R_Code/powerplant/data/so2_facility_emissions_1940_1990.csv"

# Output directories
WEIGHTED_DIR <- "/scratch/xshan2/R_Code/disperseR/weighted_fst"
PER_TON_DIR  <- "/scratch/xshan2/R_Code/disperseR/byunit_fst_perton"
dir_create(WEIGHTED_DIR, recurse = TRUE)
dir_create(PER_TON_DIR,  recurse = TRUE)

YEARS_ANN  <- 1940:1990            # Annual weighted output
YEARS_DEC  <- seq(1940, 1990, 10)  # Decadal per-ton output
P4S_ALBERS <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
GRID_SPACING_M <- 36000            # 36 km grid spacing
MAX_KM <- 300                      # Main mapping threshold (km)
FST_COMPRESS <- 50
OUT_BYUNIT_PREFIX <- "grids_exposures_byunit_"
OUT_TOTAL_PREFIX  <- "grids_exposures_total_"
# ======================================================= #


# ======================= STEP A ============================
# Build HyADS source points (uID peaks) from 2022 kernel
# ===========================================================
cat("\n[STEP A] Loading 2022 HyADS kernel and building uID peaks...\n")
kernel <- read_fst(HYADS_KERNEL_FST, as.data.table = TRUE)
stopifnot(all(c("x","y","uID","hyads") %in% names(kernel)))
kernel <- kernel[is.finite(x) & is.finite(y) & is.finite(hyads)]
kernel[, uID := as.character(uID)]
setkey(kernel, uID)
setorder(kernel, uID, -hyads)

# uID peak coordinates (Albers meters)
uID_peak <- kernel[, .SD[1L, .(uID, src_x_m = x, src_y_m = y, hyads_peak = hyads)], by = uID]
uID_sf   <- st_as_sf(uID_peak, coords = c("src_x_m","src_y_m"), crs = P4S_ALBERS, remove = FALSE)
uID_m    <- st_transform(uID_sf, 3857)

# Split kernel by uID for faster weighting
kernel_small <- kernel[, .(x, y, uID, hyads)]
kernel_small[, uID := as.character(uID)]
kernel_split <- split(kernel_small, f = kernel_small$uID, drop = TRUE)
names(kernel_split) <- as.character(unique(kernel_small$uID))
cat("Total uIDs in kernel:", length(kernel_split), "\n")


# ======================= STEP B ============================
# Map (Year, Plant_ID) → nearest HyADS uID  (Coal only)
# ===========================================================
cat("\n[STEP B] Reading emission CSV and mapping Plant_ID → uID (Coal only)...\n")
emiss <- fread(EMISS_CSV)
if ("year" %in% names(emiss))         setnames(emiss, "year", "Year")
if ("plant_id" %in% names(emiss))     setnames(emiss, "plant_id", "Plant_ID")
if ("generator_id" %in% names(emiss)) setnames(emiss, "generator_id", "Generator_ID")
if ("convert_f" %in% names(emiss))    setnames(emiss, "convert_f", "convert_fuel")

stopifnot(all(c("Year","Plant_ID","Generator_ID","convert_fuel","Longitude","Latitude") %in% names(emiss)))
emiss <- emiss[convert_fuel %chin% c("Coal","COAL","coal")]
emiss[, plant_lon := as.numeric(Longitude)]
emiss[, plant_lat := as.numeric(Latitude)]
emiss <- emiss[is.finite(plant_lon) & is.finite(plant_lat)]

# Ensure SO2_tons exists
if (!"SO2_tons" %in% names(emiss)) {
  if ("so2_kg" %in% names(emiss)) emiss[, SO2_tons := so2_kg / 907.185]
  if ("SO2_kg" %in% names(emiss)) emiss[, SO2_tons := SO2_kg / 907.185]
}
stopifnot("SO2_tons" %in% names(emiss))
setkey(emiss, Year, Plant_ID)

years_use <- sort(unique(emiss$Year))
cat("Years in emission file:", paste(head(years_use, 8), collapse = ","), "...\n")

# Map plant → nearest uID
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

qa_full <- copy(map_year_plant_uid)
qa_full[, flag_far := nn_km > MAX_KM]
fwrite(qa_full, file.path(WEIGHTED_DIR, sprintf("uID_mapping_QA_full_MAXKM%d.csv", MAX_KM)))

map_year_plant_uid_filt <- map_year_plant_uid[!is.na(uID) & nn_km <= MAX_KM]


# ======================= STEP C ============================
# Aggregate emissions to (Year, uID)
# ===========================================================
cat("\n[STEP C] Aggregating emissions to (Year, uID)...\n")

# --- aggregate to annual total tons ---
emiss_plant_year <- emiss[, .(SO2_tons = sum(SO2_tons, na.rm = TRUE)), by = .(Year, Plant_ID)]

# Aggregate to uID
emiss_uid_year <- map_year_plant_uid_filt[emiss_plant_year,
                                          on = .(Year, Plant_ID), nomatch = 0L][
  , .(SO2_tons = sum(SO2_tons, na.rm = TRUE)), by = .(Year, uID)]
emiss_uid_year[, uID := as.character(uID)]
setkey(emiss_uid_year, Year, uID)

write_fst(emiss_uid_year, file.path(WEIGHTED_DIR, "emiss_uid_year.fst"), compress = FST_COMPRESS)
cat("Saved emiss_uid_year.fst with", nrow(emiss_uid_year), "rows (annual tons).\n")


# ======================= STEP D ============================
# Compute weighted HyADS at plant (hyads_at_plant_w)
# and write annual weighted by-unit / total FSTs
# ===========================================================
cat("\n[STEP D] Computing hyads_at_plant_w and writing annual weighted FSTs...\n")

plant_xy <- unique(emiss[, .(Year, Plant_ID, plant_lon, plant_lat)])
map_year_plant_uid_filt <- map_year_plant_uid_filt[plant_xy, on = .(Year, Plant_ID), nomatch = 0L]
map_year_plant_uid_filt <- map_year_plant_uid_filt[is.finite(plant_lon) & is.finite(plant_lat)]

# Project plant coordinates to Albers
pl_sf  <- st_as_sf(map_year_plant_uid_filt, coords = c("plant_lon","plant_lat"), crs = 4326, remove = FALSE)
pl_alb <- st_transform(pl_sf, P4S_ALBERS)
pl_xy  <- st_coordinates(pl_alb)
map_year_plant_uid_filt[, `:=`(plant_x = pl_xy[,1], plant_y = pl_xy[,2])]

map_year_plant_uid_filt[, `:=`(grid_x = NA_real_, grid_y = NA_real_, hyads_at_plant = NA_real_)]
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

# Save plant-level table
fwrite(plant_map_emiss[, .(Year, Plant_ID, uID, nn_km,
                           plant_lon, plant_lat, plant_x, plant_y,
                           grid_x, grid_y, SO2_tons,
                           hyads_at_plant, hyads_at_plant_w)],
       file.path(WEIGHTED_DIR, sprintf("plant_year_hyads_at_plant_MAXKM%d.csv", MAX_KM)))

# Annual weighted by-unit / total FST
for (yy in YEARS_ANN) {
  eyy <- emiss_uid_year[Year == yy]
  if (nrow(eyy) == 0L) {
    message("Year ", yy, ": no valid plants within ", MAX_KM, " km → skip.")
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
  setcolorder(byunit_out, c("x","y","uID","year","hyads"))
  write_fst(byunit_out,
            file.path(WEIGHTED_DIR, sprintf("%s%d.fst", OUT_BYUNIT_PREFIX, yy)),
            compress = FST_COMPRESS)

  total_out <- byunit_out[, .(hyads = sum(hyads, na.rm = TRUE)), by = .(x, y)]
  total_out[, year := yy]
  setcolorder(total_out, c("x","y","year","hyads"))
  write_fst(total_out,
            file.path(WEIGHTED_DIR, sprintf("%s%d.fst", OUT_TOTAL_PREFIX, yy)),
            compress = FST_COMPRESS)

  message("Weighted FST written for ", yy)
}


# ======================= STEP E2 ===========================
# Compute grid-aligned translation offsets (plume alignment)
# ===========================================================
cat("\n[STEP E2] Computing grid-aligned offsets...\n")
src_xy <- unique(uID_peak[, .(uID, src_x_m, src_y_m)])
plant_map_emiss <- plant_map_emiss[src_xy, on = "uID"]
plant_map_emiss[, `:=`(
  ngrid_dx   = as.integer(round((plant_x - src_x_m) / GRID_SPACING_M)),
  ngrid_dy   = as.integer(round((plant_y - src_y_m) / GRID_SPACING_M)),
  shift_dx_m = ngrid_dx * GRID_SPACING_M,
  shift_dy_m = ngrid_dy * GRID_SPACING_M
)]
cat("Offset summary (km):\n")
print(summary(plant_map_emiss$shift_dx_m / 1000))
print(summary(plant_map_emiss$shift_dy_m / 1000))
fwrite(plant_map_emiss,
       file.path(WEIGHTED_DIR,
                 sprintf("plant_year_hyads_at_plant_MAXKM%d_with_shift.csv", MAX_KM)))


# ======================= STEP F ============================
# Generate decadal per-ton by-unit & total FSTs (1940–1990)
# Schema identical to 2022: x, y, uID, hyads, year.E, year.D
# ===========================================================
cat("\n[STEP F] Writing decadal per-ton by-unit & total FSTs...\n")

for (yy in YEARS_DEC) {
  active_uids <- emiss_uid_year[Year == yy & SO2_tons > 0, unique(uID)]

  dt <- copy(kernel_small)
  dt[, hyads := 0.0]
  if (length(active_uids))
    dt[uID %chin% active_uids, hyads := kernel_small[uID %chin% active_uids, hyads]]

  dt[, `:=`(year.E = yy, year.D = yy)]
  out_byunit <- file.path(PER_TON_DIR, sprintf("%s%d.fst", OUT_BYUNIT_PREFIX, yy))
  write_fst(dt[, .(x, y, uID, hyads, year.E, year.D)], out_byunit, compress = FST_COMPRESS)

  total_dt <- dt[, .(hyads = sum(hyads, na.rm = TRUE)), by = .(x, y)]
  total_dt[, `:=`(year.E = yy, year.D = yy)]
  out_total <- file.path(PER_TON_DIR, sprintf("%s%d.fst", OUT_TOTAL_PREFIX, yy))
  write_fst(total_dt[, .(x, y, hyads, year.E, year.D)], out_total, compress = FST_COMPRESS)

  message("Per-ton decadal FST written for ", yy)
}

cat("\nAll tasks completed successfully.\n")
cat("- Annual weighted FSTs & plant tables → ", WEIGHTED_DIR, "\n")
cat("- Decadal per-ton FSTs → ", PER_TON_DIR, "\n")
