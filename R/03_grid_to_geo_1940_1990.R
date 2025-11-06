#!/usr/bin/env Rscript
# ============================================================
# Convert grid-level PM2.5 FST → ZIP-level (area-weighted)
#   Source: Sherry Shan (2025)
#   Years: 1940, 1950, 1960, 1970, 1980, 1990
#   Input:  /scratch/xshan2/R_Code/disperseR/main/output/pm25_decades_model.lm.cv_single_poly_proxy1999met
#   Output: ZIP-level area-weighted PM2.5 per year (FST)
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(raster)
  library(data.table)
  library(fst)
  library(areal)
  library(USAboundaries)
})

# ============================================================
# 1. Projection (must match PM2.5 grids)
# ============================================================
p4s <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

# ============================================================
# 2. ZIP shapefile + crosswalk loader
# ============================================================
zip_sf_reader <- function(
  shapefile_dir = "/projects/HAQ_LAB/lhennem/disperseR/main/input/zcta_500k/"
) {
  zcta_shp <- file.path(shapefile_dir, "cb_2017_us_zcta510_500k.shp")

  # The crosswalk should come from disperseR package; make sure it’s available
  cw <- disperseR::crosswalk

  zips <- st_read(zcta_shp, quiet = TRUE)
  zips <- st_transform(zips, p4s)
  if ("ZCTA5CE10" %in% names(zips)) names(zips)[names(zips) == "ZCTA5CE10"] <- "ZCTA"
  zips <- merge(zips, cw, by = "ZCTA", all = FALSE, allow.cartesian = TRUE)
  zips$ZIP <- formatC(zips$ZIP, width = 5, format = "d", flag = "0")
  zips[, "ZIP"]
}

# ============================================================
# 3. Core converter: grids → ZIPs (area-weighted)
# ============================================================
grids_to_zips <- function(file.in, path.out, zips, p4s) {
  message("\nProcessing: ", file.in)
  out_base <- gsub("^grids_", "zips_", basename(file.in))
  out_path <- file.path(path.out, out_base)

  # ---- Read grid PM2.5 ----
  g <- read_fst(file.in, as.data.table = TRUE)
  if (!nrow(g)) { message("Empty grid file, skipped."); return(invisible(NULL)) }

  # Replace NAs with zeros
  for (nm in names(g)) set(g, which(is.na(g[[nm]])), nm, 0)

  # Identify value columns (PM2.5)
  value.cols <- setdiff(names(g), c("x","y"))
  if (!length(value.cols)) {
    warning("No numeric value columns found in: ", file.in)
    return(invisible(NULL))
  }

  # Skip if all-zero
  keep <- value.cols[colSums(as.matrix(g[, ..value.cols])) > 0]
  if (!length(keep)) { message("All-zero file, skipped."); return(invisible(NULL)) }
  g.trim <- g[, c("x","y", keep), with = FALSE]

  # ---- Convert grid points to polygons ----
  r <- rasterFromXYZ(g.trim)
  crs(r) <- p4s
  names(r) <- keep
  rp_sf <- st_as_sf(rasterToPolygons(r))
  rp_sf <- st_transform(rp_sf, p4s)
  rp_sf$GID <- seq_len(nrow(rp_sf))

  # ---- Reshape to long table ----
  rp_dt <- data.table(rp_sf)[, geometry := NULL]
  ncin_m <- melt(rp_dt, id.vars = "GID", variable.name = "uID", value.name = "val")

  # ---- ZIP polygons & intersection weights ----
  zips <- st_transform(zips, p4s)
  ncin_train <- rp_sf[, "GID"]
  weights <- areal::aw_intersect(zips, source = ncin_train, areaVar = "area")
  wdt <- data.table(weights)

  zips.a <- data.table(ZIP = zips$ZIP, ZIP.area = as.vector(st_area(zips)))
  wdt <- merge(wdt, zips.a, by = "ZIP")
  wdt[, areaWeight := area / ZIP.area]

  # ---- Apply weights and aggregate ----
  mm <- merge(wdt, ncin_m, by = "GID", allow.cartesian = TRUE)
  mm[, val_aw := areaWeight * val]
  zip_sums <- mm[, .(pm25 = sum(val_aw)), by = .(ZIP, uID)]

  # ---- Write ZIP-level FST ----
  out_wide <- dcast(zip_sums, ZIP ~ uID, value.var = "pm25")
  dir.create(path.out, recursive = TRUE, showWarnings = FALSE)
  write_fst(out_wide, out_path)
  message("Wrote ZIP-level FST: ", out_path)
  invisible(out_path)
}

# ============================================================
# 4. Main runner
# ============================================================
# Input directory = where your grids_pm25_total_YYYY.fst live
input_dir  <- "/scratch/xshan2/R_Code/disperseR/main/output/pm25_decades_model.lm.cv_single_poly_proxy1999met"
output_dir <- file.path(input_dir, "zips_exposures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

years <- c(1940,1950,1960,1970,1980,1990)
zips  <- zip_sf_reader()

cat("\n[RUN] Input directory :  ", input_dir,
    "\n[RUN] Output directory: ", output_dir, "\n", sep = "")

for (yr in years) {
  file.in <- file.path(input_dir, sprintf("grids_pm25_total_%d.fst", yr))
  if (!file.exists(file.in)) {
    warning("[RUN] Missing file: ", file.in, call. = FALSE)
    next
  }
  grids_to_zips(file.in, output_dir, zips, p4s)
}

cat("\n[ALL DONE] ZIP-level PM2.5 exposures written to:\n  ", output_dir, "\n")
