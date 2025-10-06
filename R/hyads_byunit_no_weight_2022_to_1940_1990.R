###############################################################
# Map coal plants (1940–1990 by decade) to nearest HyADS-2022 uID
# Output: one tidy table with year, plant info, plant lon/lat,
#         nearest HyADS uID, distance (km), hyads_peak, hyads_sum,
#         and HyADS source lon/lat.
# Author: You
# Dependencies: data.table, readxl, fst, sf
###############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(fst)
  library(sf)
})

# ----------------------- USER PATHS -----------------------
EIA860_PATH <- "/home/xshan2/HAQ_LAB/xshan2/R_Code/powerplant/december_generator2022.xlsx"
HYADS_BYUNIT_2022_FST <- "/home/xshan2/HAQ_LAB/lhennem/disperseR/main/output/exp/coal_hyads_no_weight/grids_exposures_byunit_2022.fst"
P4S_ALBERS <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
YEARS <- seq(1940, 1990, by = 10)  # decades
# ----------------------------------------------------------


# ======================= STEP A ============================
# Read & prepare HyADS-2022 (per-ton by-unit)
# - build per-uID peak point (lon/lat) and total plume strength
# ===========================================================

hyads_2022 <- read_fst(HYADS_BYUNIT_2022_FST, as.data.table = TRUE)

req_cols <- c("x","y","uID","hyads")
if (!all(req_cols %in% names(hyads_2022))) {
  stop("HyADS file missing columns: ", paste(setdiff(req_cols, names(hyads_2022)), collapse = ", "))
}

# strict NA filter and standardize types
hyads_clean <- hyads_2022[is.finite(x) & is.finite(y) & is.finite(hyads)]
hyads_clean[, uID := as.character(uID)]

# peak (max hyads) per uID for source point & total plume strength
hy_sum  <- hyads_clean[, .(hyads_sum  = sum(hyads, na.rm = TRUE)), by = uID]
hy_pos  <- hyads_clean[is.finite(hyads) & hyads >= 0]
setorder(hy_pos, uID, -hyads)
hy_peak_xy <- hy_pos[, .SD[1L, .(x, y, hyads_peak = hyads)], by = uID]

# project HyADS source points to lon/lat for reporting
hy_src_sf  <- st_as_sf(hy_peak_xy, coords = c("x","y"), crs = P4S_ALBERS, remove = FALSE)
hy_src_ll  <- st_transform(hy_src_sf, 4326)
src_coords <- st_coordinates(hy_src_ll)

hy_src <- data.table(
  uID        = hy_peak_xy$uID,
  src_x      = hy_peak_xy$x,
  src_y      = hy_peak_xy$y,
  src_lon    = src_coords[,1],
  src_lat    = src_coords[,2],
  hyads_peak = hy_peak_xy$hyads_peak
)[hy_sum, on = "uID"]  # add hyads_sum

# sf objects for fast nearest-neighbor distance
hy_src_m <- st_transform(hy_src_ll, 3857)  # meters


# ======================= STEP B ============================
# Read EIA-860 once
# ===========================================================

generator_act <- data.table(read_excel(EIA860_PATH, sheet = 1, skip = 2))
generator_ret <- data.table(read_excel(EIA860_PATH, sheet = 3, skip = 2))

cols_keep <- c("Plant ID","Plant Name","Plant State","County",
               "Generator ID","Energy Source Code",
               "Operating Year","Retirement Year",
               "Nameplate Capacity (MW)","Latitude","Longitude")

generator_act <- generator_act[, intersect(cols_keep, names(generator_act)), with = FALSE]
generator_ret <- generator_ret[, intersect(cols_keep, names(generator_ret)), with = FALSE]
if (!"Retirement Year" %in% names(generator_ret)) {
  stop("The 'Retired' sheet is missing the 'Retirement Year' column.")
}


# ======================= STEP C ============================
# Loop decades → build mapping table for each year → rbind
# ===========================================================

match_tbl_all <- rbindlist(lapply(YEARS, function(yr) {
  cat("\n================ YEAR:", yr, "================\n")

  # --- filter generators active in this year ---
  act_yr <- generator_act[`Operating Year` <= yr]
  ret_yr <- generator_ret[`Operating Year` <= yr & `Retirement Year` > yr]
  gen_yr <- rbind(act_yr, ret_yr, fill = TRUE)

  # --- keep coal only ---
  gen_yr[, convert_fuel := fifelse(`Energy Source Code` %in% c("BIT","SUB","LIG"), "Coal", "NonCoal")]
  coal_yr <- gen_yr[convert_fuel == "Coal"]

  # --- clean coords & deduplicate plants ---
  coal_yr[, plant_id := suppressWarnings(as.integer(`Plant ID`))]
  coal_yr[, Latitude  := suppressWarnings(as.numeric(Latitude))]
  coal_yr[, Longitude := suppressWarnings(as.numeric(Longitude))]
  coal_yr <- coal_yr[is.finite(Latitude) & is.finite(Longitude)]

  plants_yr <- unique(coal_yr[, .(
    plant_id, `Plant Name`, `Plant State`, County, Latitude, Longitude
  )])

  cat("Coal generators (rows): ", nrow(coal_yr), "\n")
  cat("Unique coal plants     : ", plants_yr[, uniqueN(plant_id)], "\n")

  if (nrow(plants_yr) == 0L) {
    return(data.table(
      year = integer(), plant_id = integer(), `Plant Name` = character(),
      `Plant State` = character(), County = character(),
      plant_lon = numeric(), plant_lat = numeric(),
      hyads_uID_nn = character(), nn_dist_km = numeric(),
      hyads_peak = numeric(), hyads_sum = numeric(),
      src_lon = numeric(), src_lat = numeric()
    ))
  }

  # --- build plant sf and compute nearest HyADS uID ---
  pl_sf <- st_as_sf(plants_yr, coords = c("Longitude","Latitude"), crs = 4326, remove = FALSE)
  pl_m  <- st_transform(pl_sf, 3857)

  nn_idx  <- st_nearest_feature(pl_m, hy_src_m)
  nn_dist <- st_distance(pl_m, hy_src_m[nn_idx, ], by_element = TRUE)  # units: m
  nn_km   <- as.numeric(nn_dist) / 1000

  # assemble tidy mapping with lon/lat & year
  out <- cbind(
    data.table(
      year       = yr,
      plant_id   = plants_yr$plant_id,
      `Plant Name`  = plants_yr$`Plant Name`,
      `Plant State` = plants_yr$`Plant State`,
      County        = plants_yr$County,
      plant_lon  = plants_yr$Longitude,
      plant_lat  = plants_yr$Latitude,
      hyads_uID_nn = hy_src$uID[nn_idx],
      nn_dist_km   = nn_km,
      hyads_peak   = hy_src$hyads_peak[nn_idx],
      hyads_sum    = hy_src$hyads_sum[nn_idx],
      src_lon      = hy_src$src_lon[nn_idx],
      src_lat      = hy_src$src_lat[nn_idx]
    )
  )

  # basic preview per year
  cat("Distance (km) summary:\n"); print(summary(out$nn_dist_km))
  cat("Preview (top 5 by proximity):\n")
  print(out[order(nn_dist_km)][1:min(5, .N),
        .(year, plant_id, `Plant Name`, nn_dist_km, hyads_uID_nn, hyads_peak, hyads_sum,
          plant_lon, plant_lat, src_lon, src_lat)])

  out
}), fill = TRUE, use.names = TRUE)


# ======================= STEP D ============================
# Append per-ton HYADS at the plant location (nearest grid cell)
# ===========================================================

# 1) Project plant lon/lat to the HyADS CRS (Albers) to get plant_x/plant_y
pl_all_sf  <- st_as_sf(match_tbl_all, coords = c("plant_lon","plant_lat"), crs = 4326, remove = FALSE)
pl_all_alb <- st_transform(pl_all_sf, P4S_ALBERS)
pl_xy      <- st_coordinates(pl_all_alb)

mt <- as.data.table(match_tbl_all)
mt[, `:=`(plant_x = pl_xy[,1], plant_y = pl_xy[,2])]

# 2) Prep HyADS grid (only needed cols; ensure character uID)
hc <- hyads_clean[, .(uID = as.character(uID), x, y, hyads)]

# 3) For each row: within its matched uID, find nearest (x,y) and fetch hyads
mt[, `:=`(hyads_at_plant = NA_real_, grid_x = NA_real_, grid_y = NA_real_)]
uids_needed <- unique(mt$hyads_uID_nn)

for (uid in uids_needed) {
  sub <- hc[uID == uid]
  idx <- which(mt$hyads_uID_nn == uid)
  if (nrow(sub) == 0L || length(idx) == 0L) next
  # per row, pick nearest grid cell
  for (i in idx) {
    d2 <- (sub$x - mt$plant_x[i])^2 + (sub$y - mt$plant_y[i])^2
    j  <- which.min(d2)
    mt$hyads_at_plant[i] <- sub$hyads[j]
    mt$grid_x[i]         <- sub$x[j]
    mt$grid_y[i]         <- sub$y[j]
  }
}

# 4) Replace the master table
match_tbl_all <- mt[]

# ======================= STEP E ============================
# Final preview & (optional) write
# ===========================================================

cat("\n================ FINAL PREVIEW (top 10) ================\n")
print(match_tbl_all[order(year, nn_dist_km)][1:min(10, .N),
      .(year, plant_id, `Plant Name`, hyads_uID_nn, nn_dist_km,
        hyads_at_plant,  # <-- per-ton hyads at plant location
        hyads_peak, hyads_sum,
        plant_lon, plant_lat, src_lon, src_lat,
        grid_x, grid_y)])

cat("\nRows in final table:", nrow(match_tbl_all), "\n")
cat("Years covered:", paste(sort(unique(match_tbl_all$year)), collapse = ", "), "\n")

# Optional: write to CSV
data.table::fwrite(match_tbl_all, "/scratch/xshan2/R_Code/disperseR/coal_plant_to_hyads_mapping_1940_1990.csv")
