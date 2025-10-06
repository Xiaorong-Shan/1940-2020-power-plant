###############################################################
# Map 1940 coal plants to nearest HyADS-2022 uID
# Author: Xiaorong Shan
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
# ----------------------------------------------------------


# ======================= STEP 1 ============================
# Read and filter 1940 coal facilities (generator level)
# ===========================================================

generator_act <- data.table(read_excel(EIA860_PATH, sheet = 1, skip = 2))
generator_ret <- data.table(read_excel(EIA860_PATH, sheet = 3, skip = 2))

cols_keep <- c("Plant ID","Plant Name","Plant State","County",
               "Generator ID","Energy Source Code",
               "Operating Year","Retirement Year",
               "Nameplate Capacity (MW)","Latitude","Longitude")

generator_act <- generator_act[, intersect(cols_keep, names(generator_act)), with = FALSE]
generator_ret <- generator_ret[, intersect(cols_keep, names(generator_ret)), with = FALSE]

act40 <- generator_act[`Operating Year` <= 1940]
if (!"Retirement Year" %in% names(generator_ret)) stop("Missing 'Retirement Year' in retired sheet.")
ret40 <- generator_ret[`Operating Year` <= 1940 & `Retirement Year` > 1940]

generator1940 <- rbind(act40, ret40, fill = TRUE)
generator1940[, convert_fuel := fifelse(`Energy Source Code` %in% c("BIT","SUB","LIG"), "Coal", "NonCoal")]
coal1940 <- generator1940[convert_fuel == "Coal"]

coal1940[, plant_id := suppressWarnings(as.integer(`Plant ID`))]
coal1940[, Latitude  := suppressWarnings(as.numeric(Latitude))]
coal1940[, Longitude := suppressWarnings(as.numeric(Longitude))]
coal1940 <- coal1940[is.finite(Latitude) & is.finite(Longitude)]

plants1940 <- unique(coal1940[, .(
  plant_id, `Plant Name`, `Plant State`, County, Latitude, Longitude
)])

cat("Rows (generator-level, coal, operating in 1940): ", nrow(coal1940), "\n")
cat("Unique plants: ", plants1940[, uniqueN(plant_id)], "\n")
print(coal1940[, .N, by = `Energy Source Code`][order(-N)])
cat("Plant-level preview:\n")
print(head(plants1940, 10))


# ======================= STEP 2 ============================
# Read and clean HyADS 2022 per-ton by-unit grid
# ===========================================================

hyads_2022 <- read_fst(HYADS_BYUNIT_2022_FST, as.data.table = TRUE)

req_cols <- c("x", "y", "uID", "hyads")
if (!all(req_cols %in% names(hyads_2022))) {
  stop("Missing required columns in HyADS file: ",
       paste(setdiff(req_cols, names(hyads_2022)), collapse = ", "))
}

hyads_clean <- hyads_2022[is.finite(x) & is.finite(y) & is.finite(hyads)]
hyads_clean[, uID := as.character(uID)]

cat("\nRows before cleaning:", nrow(hyads_2022),
    "\nRows after cleaning :", nrow(hyads_clean), "\n")
cat("HyADS value summary:\n"); print(summary(hyads_clean$hyads))
cat("HyADS preview:\n"); print(head(hyads_clean, 5))


# ======================= STEP 3 ============================
# Convert HyADS Albers coordinates to WGS84 lon/lat
# ===========================================================

hy_sf  <- st_as_sf(hyads_clean, coords = c("x","y"), crs = P4S_ALBERS, remove = FALSE)
hy_lls <- st_transform(hy_sf, 4326)
coords <- st_coordinates(hy_lls)
hyads_ll <- as.data.table(hy_lls)
hyads_ll[, `:=`(lon = coords[,1], lat = coords[,2])]

cat("\nLongitude range:\n"); print(summary(hyads_ll$lon))
cat("Latitude range:\n"); print(summary(hyads_ll$lat))
cat("HyADS lon/lat preview:\n"); print(hyads_ll[1:5, .(uID, hyads, lon, lat)])


# ======================= STEP 4 ============================
# Extract one peak (max hyads) grid cell per uID
# ===========================================================

hy_pos <- hyads_ll[is.finite(hyads) & hyads >= 0 & is.finite(lon) & is.finite(lat)]
setorder(hy_pos, uID, -hyads)
hy_peak <- hy_pos[, .SD[1L, .(lon, lat, hyads_peak = hyads)], by = uID]

hy_sum <- hyads_clean[, .(hyads_sum = sum(hyads, na.rm = TRUE)), by = uID]
hy_src <- merge(hy_peak, hy_sum, by = "uID", all.x = TRUE)

cat("\nuIDs with a peak point:", nrow(hy_src), "\n")
cat("Peak/strength preview:\n")
print(head(hy_src))


# ======================= STEP 5 ============================
# Nearest-neighbor mapping: 1940 plant → closest HyADS uID
# ===========================================================

hy_peak_sf <- st_as_sf(hy_src, coords = c("lon","lat"), crs = 4326, remove = FALSE)
plants1940_sf <- st_as_sf(
  plants1940[is.finite(Latitude) & is.finite(Longitude)],
  coords = c("Longitude","Latitude"), crs = 4326, remove = FALSE
)

# Transform to meters (Web Mercator) for distance
hy_m <- st_transform(hy_peak_sf, 3857)
pl_m <- st_transform(plants1940_sf, 3857)

nn_idx  <- st_nearest_feature(pl_m, hy_m)
nn_dist <- st_distance(pl_m, hy_m[nn_idx, ], by_element = TRUE)
nn_km   <- as.numeric(nn_dist) / 1000

match_tbl <- cbind(
  st_drop_geometry(plants1940_sf),
  data.frame(
    hyads_uID_nn = hy_src$uID[nn_idx],
    nn_dist_km   = nn_km
  )
)

match_tbl <- merge(
  match_tbl,
  hy_src[, .(uID, hyads_peak, hyads_sum)],
  by.x = "hyads_uID_nn", by.y = "uID", all.x = TRUE
)

cat("\nMatched plants:", nrow(match_tbl),
    "(should equal number of input plants)\n")
cat("Distance (km) summary:\n"); print(summary(match_tbl$nn_dist_km))
cat("\nNearest mapping preview (top 10 by proximity):\n")
print(match_tbl[order(nn_dist_km)][1:10, .(
  plant_id, `Plant Name`, hyads_uID_nn, nn_dist_km, hyads_peak, hyads_sum
)])

# ---------------------- OPTIONAL OUTPUTS --------------------
# fwrite(match_tbl, "map1940_to_hyads2022_nearest.csv")
# fwrite(hy_src,    "hyads2022_uID_peak_and_sum.csv")
# ------------------------------------------------------------

cat("\nDone.\nYou now have:\n")
cat(" - `match_tbl`: each 1940 coal plant matched to its nearest HyADS-2022 uID\n")
cat(" - `hy_src`: per-uID plume peak and total hyads (for diagnostics or weighting)\n")
