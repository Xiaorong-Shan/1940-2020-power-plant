###############################################################################
# Script: national_PWE_decadal.R
# Author: Xiaorong Shan
# Purpose:
#   Calculate national and county-level population-weighted PM2.5 exposure (PWE)
#   for the United States from 1940 to 1990 (decadal intervals).
#
# Description:
#   - Input: 
#       (1) Decadal HyADS-modeled PM2.5 concentration grids (.fst files, 36 km)
#           Columns: x, y, vals.out (in Albers Equal Area projection)
#       (2) NHGIS county-level population time series (A00AA1940–A00AA1990)
#   - Process:
#       For each decade, the script:
#         1. Loads PM2.5 raster for that year.
#         2. Aggregates 36-km grid cells to county averages (area-weighted).
#         3. Merges county population for the same year.
#         4. Computes national PWE = Σ(P × E) / Σ(P).
#         5. Saves county-level results and prints the national PWE table.
#   - Output:
#       • county_pm25_pop_YYYY.csv : County-level PM2.5 and population data
#       • PWE_trend_1940_1990.png  : National PWE trend plot
#       • Printed PWE table in console
###############################################################################

# ================== LOAD PACKAGES ==================
library(raster)
library(sf)
library(USAboundaries)
library(data.table)
library(dplyr)
library(ggplot2)
library(fst)

# ================== PATHS & PARAMETERS ==================
# HyADS PM2.5 decadal outputs (.fst) — 36 km grid, Albers projection
input_dir <- "/scratch/xshan2/R_Code/disperseR/main/output/pm25_decades_model.lm.cv_single_poly_proxy1999met"

# NHGIS population time series (contains A00AA1940 … A00AA1990)
pop_file  <- "/scratch/xshan2/R_Code/powerplant/data/population/nhgis0012_csv/nhgis0012_ts_nominal_county.csv"

# Output directory
out_dir   <- file.path(input_dir, "pwe_outputs_decades")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Decades to process
years <- c(1940, 1950, 1960, 1970, 1980, 1990)

# Albers Equal-Area projection (same as your grid)
p4s <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

# ================== 1) LOAD COUNTY BOUNDARIES ==================
# Get U.S. county shapefile (low-resolution is fine)
counties_sf <- us_counties()

# Remove duplicate columns (some usboundaries versions duplicate names)
counties_sf <- counties_sf[, !duplicated(names(counties_sf))]

# Keep only essential columns
counties_sf <- counties_sf %>%
  select(geoid, statefp, countyfp, state_name, state_abbr, name, geometry)

# Exclude non-CONUS regions (AK, HI, PR)
counties_sf <- counties_sf %>% filter(!(state_abbr %in% c("AK", "HI", "PR")))

# Reproject to match PM2.5 grid CRS
counties_sf <- st_transform(counties_sf, crs = p4s)

# ================== 2) LOAD NHGIS POPULATION DATA ==================
# Read NHGIS population time series (wide format)
pop_ts <- fread(pop_file)

# Expected columns — adjust if different in your file
cols_needed <- c("STATEFP","COUNTYFP","A00AA1940","A00AA1950",
                 "A00AA1960","A00AA1970","A00AA1980","A00AA1990")
stopifnot(all(cols_needed %in% names(pop_ts)))

# Convert from wide to long format → geoid, YEAR, pop
pop_long <- melt(
  pop_ts[, ..cols_needed],
  id.vars = c("STATEFP","COUNTYFP"),
  variable.name = "VAR", value.name = "pop"
)[, `:=`(
  geoid = paste0(sprintf("%02d", as.integer(STATEFP)),
                 sprintf("%03d", as.integer(COUNTYFP))),
  YEAR  = as.integer(sub("A00AA", "", VAR))
)][, .(geoid, YEAR, pop)]

# Keep only target decades
pop_dt <- pop_long[YEAR %in% years]

# ================== 3) HELPER FUNCTION: FST → RASTER ==================
# Convert .fst (x, y, vals.out) into a RasterLayer in Albers CRS
make_raster_albers <- function(fst_path, p4s) {
  dt <- as.data.table(read_fst(fst_path))
  if (!all(c("x","y","vals.out") %in% names(dt))) {
    stop("FST columns must include x, y, vals.out. Found: ", paste(names(dt), collapse=", "))
  }
  rasterFromXYZ(dt[, .(x, y, vals.out)], crs = CRS(p4s))
}

# ================== 4) MAIN LOOP ==================
pwe_list   <- list()
county_all <- list()

for (yy in years) {
  message("====== Processing year: ", yy, " ======")
  fpath <- file.path(input_dir, sprintf("grids_pm25_total_%d.fst", yy))
  if (!file.exists(fpath)) {
    warning("Missing ", fpath, ", skipping.")
    pwe_list[[as.character(yy)]] <- data.frame(year = yy, PWE = NA_real_)
    next
  }

  # 4.1 Read PM2.5 raster
  r_yy <- make_raster_albers(fpath, p4s)

  # 4.2 Compute area-weighted county means
  v <- raster::extract(
    r_yy, counties_sf,
    fun = mean, na.rm = TRUE,
    weights = TRUE, exact = FALSE
  )

  # 4.3 Build county table
  county_out <- counties_sf %>%
    st_drop_geometry() %>%
    select(geoid, state_name, name) %>%
    mutate(year = yy, pm25 = as.numeric(v))

  # 4.4 Merge population for the same year
  pop_year_dt <- pop_dt[YEAR == yy]
  county_out  <- left_join(county_out, pop_year_dt, by = "geoid") %>%
                 rename(pop_year = pop)

  # 4.5 Compute national population-weighted exposure
  pwe_val <- if (nrow(pop_year_dt)) {
    sum(county_out$pm25 * county_out$pop_year, na.rm = TRUE) /
      sum(county_out$pop_year, na.rm = TRUE)
  } else NA_real_

  # 4.6 Save county-level file
  fwrite(county_out, file.path(out_dir, sprintf("county_pm25_pop_%d.csv", yy)))

  # Store results
  pwe_list[[as.character(yy)]]   <- data.frame(year = yy, PWE = pwe_val)
  county_all[[as.character(yy)]] <- county_out

  message("Year ", yy, " PWE = ", signif(pwe_val, 4))
}

# ================== 5) PRINT NATIONAL RESULTS ==================
pwe_df <- rbindlist(pwe_list)

cat("\n========== National Population-Weighted PM2.5 (1940–1990) ==========\n")
print(format(pwe_df, digits = 6, scientific = TRUE))
cat("====================================================================\n")
county_all_df <- rbindlist(county_all, fill = TRUE)
fwrite(county_all_df, file.path(out_dir, "county_pm25_pop_all_years.csv"))

# ================== 6) PLOT NATIONAL TREND ==================
p <- ggplot(pwe_df, aes(x = year, y = PWE)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  theme_minimal(base_size = 13) +
  labs(title = "National Population-Weighted PM2.5 Exposure (1940–1990)",
       x = "Year", y = expression(PWE~(mu*g/m^3)))

ggsave(file.path(out_dir, "PWE_trend_1940_1990.pdf"), p, width = 7.2, height = 4.2, dpi = 300)

message("✅ Done. Results stored in: ", out_dir)
