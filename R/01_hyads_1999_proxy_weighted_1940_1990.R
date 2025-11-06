#!/usr/bin/env Rscript
# =============================================================================
# Exposures for decadal years only: 1940, 1950, 1960, 1970, 1980, 1990
# Proxied by 1999 monthly kernels (HyADS counts) with 36-km grid-aligned shifts.
# Consistency rules:
#   - If a target year has no emissions rows: write empty by-unit & total FSTs.
#   - If some (uID_src, MM) kernel files are missing: those months contribute 0;
#     outputs are still written for that year to keep structure consistent.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(fst); library(sf)
})

# ========= USER PATHS (EDIT AS NEEDED) =========
GL_1999_DIR <- "/scratch/xshan2/R_Code/disperseR/gridlinks_1999_unpacked/gridlinks1999_2018"
UNITS_RDA   <- "/projects/HAQ_LAB/mrasel/R/ampd-raw-data-processing/data/units_coal_1997_2021.rda"
EMISS_CSV   <- "/scratch/xshan2/R_Code/powerplant/data/so2_facility_emissions_1940_1990.csv"

SCRATCH_ROOT <- "/scratch/xshan2/R_Code/disperseR"
XWALK_FST    <- file.path(SCRATCH_ROOT, "crosswalk/hist_to_1999_crosswalk.fst")
OUT_DIR_EXP  <- file.path(SCRATCH_ROOT, "main/output/exposures_proxy1999_shifted_decades")

dir.create(dirname(XWALK_FST), recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_EXP,        recursive = TRUE, showWarnings = FALSE)

# ========= GRID / CRS =========
GRID_SPACING_M <- 36000
CRS_ALBERS     <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +ellps=GRS80 +datum=NAD83 +units=m"
MAX_KM         <- 300

# ========= HELPERS =========
yearmons <- function(yy) sprintf("%04d%02d", yy, 1:12)
id_from_bn <- function(bn) sub("^gridlinks_([^_]+)_.*$", "\\1", bn)
ym_from_bn <- function(bn) sub("^gridlinks_[^_]+_([0-9]{4})-([0-9]{2})-.*$", "\\1\\2", bn)

write_empty_outputs <- function(year, out_dir) {
  byu <- data.table(x=numeric(0), y=numeric(0), uID=character(0), year=integer(0), hyads=numeric(0))
  tot <- data.table(x=numeric(0), y=numeric(0), year=integer(0), hyads=numeric(0))
  f_byu <- file.path(out_dir, sprintf("grids_exposures_byunit_%d.fst", year))
  f_tot <- file.path(out_dir, sprintf("grids_exposures_total_%d.fst", year))
  write_fst(byu, f_byu, compress = 50)
  write_fst(tot, f_tot, compress = 50)
  message(sprintf("⚠️  Year %d has no usable data; wrote empty outputs:\n  %s\n  %s", year, f_byu, f_tot))
}

# ========= 1) BUILD / LOAD CROSSWALK =========
if (!file.exists(XWALK_FST)) {
  message("🔧 Building crosswalk (historical facility → 1999.ID)...")

  em <- fread(EMISS_CSV)
  setnames(em, names(em), tolower(names(em)))
  stopifnot("year" %in% names(em))

  if (!"so2_tons" %in% names(em)) {
    if ("so2_short_tons" %in% names(em)) setnames(em, "so2_short_tons", "so2_tons")
    if ("so2_kg" %in% names(em)) em[, so2_tons := so2_kg / 907.185]
  }
  stopifnot("so2_tons" %in% names(em))

  if (!("longitude" %in% names(em) && "latitude" %in% names(em))) {
    cand_lon <- intersect(names(em), c("lon","lng","plant_lon","longitude"))
    cand_lat <- intersect(names(em), c("lat","plant_lat","latitude"))
    if (length(cand_lon) >= 1 && length(cand_lat) >= 1) {
      setnames(em, cand_lon[1], "longitude"); setnames(em, cand_lat[1], "latitude")
    } else stop("Historical emissions are missing longitude/latitude.")
  }
  em[, `:=`(longitude = as.numeric(longitude), latitude = as.numeric(latitude))]
  em <- em[is.finite(longitude) & is.finite(latitude)]

  nm <- names(em)
  pid <- nm[grepl("^plant", nm)][1]
  gid <- nm[grepl("generator|unit", nm)][1]
  stopifnot(!is.na(pid)); setnames(em, pid, "plant_id")
  if (!is.na(gid)) setnames(em, gid, "generator_id")
  if ("generator_id" %in% names(em)) {
    em[, uID_hist := ifelse(is.na(generator_id) | generator_id=="",
                            as.character(plant_id),
                            paste0(plant_id, ".", generator_id))]
  } else {
    em[, uID_hist := as.character(plant_id)]
  }

  load(UNITS_RDA)
  if (!exists("units_updated")) {
    obj <- ls(); obj <- obj[grepl("units", obj, ignore.case = TRUE)][1]
    stopifnot(length(obj) == 1); units_updated <- get(obj)
  }
  u <- as.data.table(units_updated)
  u <- u[year == 1999 & is.finite(Latitude) & is.finite(Longitude)]
  stopifnot(nrow(u) > 0)

  if (!("ID" %in% names(u))) {
    cand_unit <- intersect(names(u), c("unit","Generator","generator","GenID","UnitID"))
    if (length(cand_unit)) {
      u[, ID := paste0(PlantId, ".", get(cand_unit[1]))]
    } else stop("No 'ID' column in 1999 unit metadata and cannot infer one.")
  }
  u[, ID := as.character(ID)]

  u_sf   <- st_as_sf(u, coords = c("Longitude","Latitude"), crs = 4326, remove = FALSE)
  u_aea  <- st_transform(u_sf, CRS_ALBERS)
  u_xy   <- st_coordinates(u_aea)
  u[, `:=`(src_x_m = u_xy[,1], src_y_m = u_xy[,2])]
  u_m    <- st_transform(u_sf, 3857)

  hist_sites <- unique(em[, .(uID_hist, longitude, latitude)], by = "uID_hist")
  hs_sf  <- st_as_sf(hist_sites, coords = c("longitude","latitude"), crs = 4326, remove = FALSE)
  hs_m   <- st_transform(hs_sf, 3857)
  hs_aea <- st_transform(hs_sf, CRS_ALBERS)
  hs_xy  <- st_coordinates(hs_aea)
  hist_sites[, `:=`(plant_x = hs_xy[,1], plant_y = hs_xy[,2])]

  idx  <- st_nearest_feature(hs_m, u_m)
  d_km <- as.numeric(st_distance(hs_m, u_m[idx,], by_element = TRUE)) / 1000

  xwalk <- data.table(
    uID_hist = hist_sites$uID_hist,
    uID_src  = u$ID[idx],
    nn_km    = d_km,
    plant_x  = hist_sites$plant_x,
    plant_y  = hist_sites$plant_y
  )[nn_km <= MAX_KM]

  xwalk <- xwalk[u[, .(uID_src = ID, src_x_m, src_y_m)], on="uID_src", nomatch=0L]
  write_fst(xwalk, XWALK_FST, compress = 50)
  message("✅ Crosswalk saved: ", XWALK_FST, "  rows=", nrow(xwalk))
} else {
  message("ℹ️ Reusing existing crosswalk: ", XWALK_FST)
  xwalk <- as.data.table(read_fst(XWALK_FST))
}

# ========= 2) INDEX 1999 MONTHLY KERNELS (BY MM) =========
files_1999 <- list.files(GL_1999_DIR, pattern="\\.fst$", full.names=TRUE)
stopifnot(length(files_1999) > 0)

idx_dt <- data.table(
  file   = files_1999,
  ID     = id_from_bn(basename(files_1999)),
  YYYYMM = ym_from_bn(basename(files_1999))
)
idx_dt <- idx_dt[substr(YYYYMM,1,4) == "1999"]
idx_dt[, MM := substr(YYYYMM,5,6)]
idx_dt[, YYYYMM := NULL]

# ========= 3) PRECOMPUTE GRID-ALIGNED SHIFTS =========
xw <- copy(xwalk)
xw[, `:=`(
  shift_dx = as.integer(round((plant_x - src_x_m) / GRID_SPACING_M)) * GRID_SPACING_M,
  shift_dy = as.integer(round((plant_y - src_y_m) / GRID_SPACING_M)) * GRID_SPACING_M
)]

# ========= 4) LOOP OVER DECADAL YEARS ONLY =========
EM_ALL <- fread(EMISS_CSV)
setnames(EM_ALL, names(EM_ALL), tolower(names(EM_ALL)))

if (!"so2_tons" %in% names(EM_ALL)) {
  if ("so2_short_tons" %in% names(EM_ALL)) setnames(EM_ALL, "so2_short_tons", "so2_tons")
  if ("so2_kg" %in% names(EM_ALL)) EM_ALL[, so2_tons := so2_kg / 907.185]
}
stopifnot("so2_tons" %in% names(EM_ALL))

if (!("longitude" %in% names(EM_ALL) && "latitude" %in% names(EM_ALL))) {
  cand_lon <- intersect(names(EM_ALL), c("lon","lng","plant_lon","longitude"))
  cand_lat <- intersect(names(EM_ALL), c("lat","plant_lat","latitude"))
  if (length(cand_lon) >= 1 && length(cand_lat) >= 1) {
    setnames(EM_ALL, cand_lon[1], "longitude"); setnames(EM_ALL, cand_lat[1], "latitude")
  } else stop("Historical emissions are missing longitude/latitude.")
}
EM_ALL[, `:=`(longitude = as.numeric(longitude), latitude = as.numeric(latitude))]

nm_all  <- names(EM_ALL)
pid_all <- nm_all[grepl("^plant", nm_all)][1]
gid_all <- nm_all[grepl("generator|unit", nm_all)][1]
stopifnot(!is.na(pid_all)); setnames(EM_ALL, pid_all, "plant_id")
if (!is.na(gid_all)) setnames(EM_ALL, gid_all, "generator_id")
if ("generator_id" %in% names(EM_ALL)) {
  EM_ALL[, uID_hist := ifelse(is.na(generator_id) | generator_id=="",
                              as.character(plant_id),
                              paste0(plant_id, ".", generator_id))]
} else {
  EM_ALL[, uID_hist := as.character(plant_id)]
}

YEARS <- seq(1940, 1990, by = 10)  # <-- only decadal years

for (YR in YEARS) {
  message("\n====== YEAR ", YR, " ======")
  ea <- EM_ALL[year == YR & is.finite(longitude) & is.finite(latitude)]
  if (!nrow(ea)) { write_empty_outputs(YR, OUT_DIR_EXP); next }

  ym <- yearmons(YR)
  em_month <- ea[, .(uID_hist,
                     YYYYMM    = rep(ym, each = .N),
                     so2_month = rep(so2_tons / 12, each = 12))]
  em_month[, MM := substr(YYYYMM, 5, 6)]

  emx <- merge(em_month, xw, by = "uID_hist", all.x = TRUE)
  emx <- emx[!is.na(uID_src)]
  if (!nrow(emx)) { write_empty_outputs(YR, OUT_DIR_EXP); next }

  emx <- merge(emx, idx_dt, by.x = c("uID_src","MM"), by.y = c("ID","MM"), all.x = TRUE)
  emx_nonNA <- emx[!is.na(file)]
  setorder(emx_nonNA, uID_src, MM, file)
  emx_nonNA <- emx_nonNA[!duplicated(emx_nonNA, by = c("uID_src","MM","uID_hist"))]

  OUT_BYUNIT <- file.path(OUT_DIR_EXP, sprintf("grids_exposures_byunit_%d.fst", YR))
  OUT_TOTAL  <- file.path(OUT_DIR_EXP, sprintf("grids_exposures_total_%d.fst",  YR))

  if (!nrow(emx_nonNA)) { write_empty_outputs(YR, OUT_DIR_EXP); next }

  parts  <- vector("list", nrow(emx_nonNA))
  n_used <- 0L
  for (i in seq_len(nrow(emx_nonNA))) {
    fp   <- emx_nonNA$file[i]
    uidH <- emx_nonNA$uID_hist[i]
    tons <- emx_nonNA$so2_month[i]
    sdx  <- emx_nonNA$shift_dx[i]
    sdy  <- emx_nonNA$shift_dy[i]

    dt <- tryCatch(as.data.table(read_fst(fp)), error = function(e) NULL)
    if (is.null(dt) || !nrow(dt)) next

    nm <- names(dt)
    if (!("x" %in% nm && "y" %in% nm)) next
    wc <- intersect(nm, c("N","n","counts","count","weight","w","hits","val","value"))
    if (!length(wc)) next
    setnames(dt, wc[1], "N")

    dt <- dt[is.finite(x) & is.finite(y) & is.finite(N)]
    if (!nrow(dt)) next

    dt[, `:=`(
      x     = x + sdx,
      y     = y + sdy,
      uID   = uidH,
      hyads = N * tons
    )][, N := NULL]

    parts[[i]] <- dt; n_used <- n_used + 1L
    if (n_used %% 400 == 0) message("  processed ", n_used, " / ", nrow(emx_nonNA))
  }

  if (n_used == 0) { write_empty_outputs(YR, OUT_DIR_EXP); next }

  byunit <- rbindlist(parts, use.names = TRUE, fill = TRUE)
  byunit <- byunit[is.finite(x) & is.finite(y) & is.finite(hyads)]
  if (!nrow(byunit)) { write_empty_outputs(YR, OUT_DIR_EXP); next }

  byunit[, year := YR]
  setcolorder(byunit, c("x","y","uID","year","hyads"))
  write_fst(byunit, OUT_BYUNIT, compress = 50)

  total <- byunit[, .(hyads = sum(hyads, na.rm = TRUE)), by = .(x, y)]
  total <- total[is.finite(x) & is.finite(y) & is.finite(hyads)]
  total[, year := YR]
  write_fst(total[, .(x, y, year, hyads)], OUT_TOTAL, compress = 50)

  qv <- tryCatch(quantile(total$hyads, probs=c(0,.25,.5,.9,.99,1), na.rm=TRUE),
                 error = function(e) rep(NA_real_, 6))
  cat(sprintf("✅ %d done:\n  %s\n  %s\n", YR, OUT_BYUNIT, OUT_TOTAL))
  if (all(is.finite(qv))) {
    cat(sprintf("  Range: min=%.3g p25=%.3g p50=%.3g p90=%.3g p99=%.3g max=%.3g\n",
                qv[1], qv[2], qv[3], qv[4], qv[5], qv[6]))
  } else cat("  (empty or non-finite QC stats)\n")
}
