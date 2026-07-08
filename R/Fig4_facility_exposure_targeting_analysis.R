#!/usr/bin/env Rscript

###############################################################################
# Full first-run script:
# Facility-level exposure concentration and targeting analysis
#
# This script:
#   1. Loads by-unit PM2.5 gridded FST files.
#   2. Area-weights gridded PM2.5 to counties.
#   3. Applies county population weights.
#   4. Computes unit-year exposure contributions.
#   5. Saves the intermediate unit-year contribution file.
#   6. Loads SO2 emissions.
#   7. Aggregates unit-level contributions to facility level.
#   8. Generates Fig. 4A, Fig. 4B, and Fig. 4C as separate PDFs.
#
# Important:
#   - This is the complete first-run script.
#   - It will run the slow by-unit extraction step.
#   - If the intermediate CSV already exists and reuse_intermediate = TRUE,
#     the script will skip the slow extraction and continue from the saved file.
#   - All comments are written in English.
#   - Figures do not include plot titles.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
  library(sf)
  library(raster)
  library(USAboundaries)
  library(ggplot2)
  library(scales)
  library(methods)
})

# ===================== USER SETTINGS =====================

years <- c(1940, 1950, 1960, 1970, 1980, 1990)

pm25_dir <- "/scratch/xshan2/R_Code/disperseR/main/output/pm25_decades_model.lm.cv_single_poly_proxy1999met"

pop_file <- "/scratch/xshan2/R_Code/powerplant/data/population/nhgis0013_ts_nominal_county.csv"

so2_file <- "/scratch/xshan2/R_Code/disperseR/pp/pp_emissions_1940_1990/so2_facility_emissions_1940_1990.csv"

out_dir <- file.path(pm25_dir, "facility_concentration_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

unit_year_file <- file.path(
  out_dir,
  "unit_year_exposure_contribution_1940_1990.csv"
)

p4s <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

# If TRUE and the intermediate file already exists, skip the slow extraction.
reuse_intermediate <- TRUE

# Number of unit columns processed per batch.
batch_size <- 100

# Number of facilities shown in Fig. 4B.
top_n_facilities <- 20

# Facility targeting fractions used in Fig. 4C.
target_cutoffs <- c(0.01, 0.05, 0.10, 0.20, 0.30)

cat("Output directory:\n", out_dir, "\n\n")
cat("Population file:\n", pop_file, "\n\n")
cat("SO2 file:\n", so2_file, "\n\n")

if (!file.exists(pop_file)) stop("Population file does not exist: ", pop_file)
if (!file.exists(so2_file)) stop("SO2 file does not exist: ", so2_file)

# ===================== HELPER FUNCTIONS =====================

fmt_geoid5 <- function(x) {
  x <- as.character(x)
  x <- gsub("\\D", "", x)
  x <- sprintf("%05d", as.integer(x))
  x
}

safe_numeric <- function(x) {
  as.numeric(gsub(",", "", as.character(x)))
}

make_uid <- function(plant_id, generator_id) {
  make.names(paste0("X", plant_id, ".", generator_id))
}

first_nonmissing_chr <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

first_nonmissing_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(as.numeric(NA))
  x[1]
}

get_unit_cols <- function(dt) {
  setdiff(names(dt), c("x", "y", "year", "vals.out", "hyads"))
}

weighted_county_means_for_batch <- function(grid_dt, unit_cols_batch, counties_sf, p4s) {

  sub_dt <- grid_dt[, c("x", "y", unit_cols_batch), with = FALSE]

  r <- raster::rasterFromXYZ(sub_dt)
  raster::crs(r) <- p4s
  names(r) <- unit_cols_batch

  counties_sp <- methods::as(counties_sf, "Spatial")

  extracted <- raster::extract(
    r,
    counties_sp,
    weights = TRUE,
    normalizeWeights = TRUE,
    na.rm = FALSE
  )

  out <- matrix(
    NA_real_,
    nrow = nrow(counties_sf),
    ncol = length(unit_cols_batch),
    dimnames = list(counties_sf$geoid, unit_cols_batch)
  )

  for (i in seq_along(extracted)) {

    m <- extracted[[i]]

    if (is.null(m)) next

    m <- as.data.frame(m)

    if (nrow(m) == 0) next

    weight_col <- intersect(names(m), c("weight", "weights"))[1]

    if (is.na(weight_col)) {
      weight_col <- names(m)[ncol(m)]
    }

    w <- as.numeric(m[[weight_col]])
    layer_cols <- intersect(unit_cols_batch, names(m))

    if (length(layer_cols) == 0) next

    vals <- as.matrix(m[, layer_cols, drop = FALSE])
    vals[!is.finite(vals)] <- NA_real_

    denom <- sum(w[is.finite(w)], na.rm = TRUE)

    if (!is.finite(denom) || denom <= 0) next

    weighted_vals <- colSums(vals * w, na.rm = TRUE) / denom

    out[i, layer_cols] <- weighted_vals
  }

  out
}

# ===================== 1) LOAD COUNTY GEOMETRY =====================

cat("Loading county geometry...\n")

counties_sf <- USAboundaries::us_counties()
counties_sf <- counties_sf[, !duplicated(names(counties_sf))]

keep_cols <- intersect(
  c("geoid", "state_abbr", "stusps", "state_name", "name", "geometry"),
  names(counties_sf)
)

counties_sf <- counties_sf[, keep_cols]

if (!("state_abbr" %in% names(counties_sf)) && ("stusps" %in% names(counties_sf))) {
  counties_sf$state_abbr <- counties_sf$stusps
}

if (!("state_name" %in% names(counties_sf))) {
  counties_sf$state_name <- NA_character_
}

if (!("name" %in% names(counties_sf))) {
  counties_sf$name <- NA_character_
}

# Keep CONUS only.
counties_sf <- counties_sf[!(counties_sf$state_abbr %in% c("AK", "HI", "PR")), ]
counties_sf$geoid <- fmt_geoid5(counties_sf$geoid)
counties_sf <- sf::st_transform(counties_sf, crs = p4s)

cat("CONUS counties:", nrow(counties_sf), "\n\n")

# ===================== 2) LOAD COUNTY POPULATION =====================

cat("Loading NHGIS county population...\n")

pop_ts <- fread(pop_file)

pop_cols <- paste0("A00AA", years)
cols_needed <- c("STATEFP", "COUNTYFP", pop_cols)

missing_cols <- setdiff(cols_needed, names(pop_ts))

if (length(missing_cols) > 0) {
  stop(
    "Population file is missing required columns:\n",
    paste(missing_cols, collapse = ", ")
  )
}

pop_long <- melt(
  pop_ts[, ..cols_needed],
  id.vars = c("STATEFP", "COUNTYFP"),
  variable.name = "VAR",
  value.name = "pop_year"
)

pop_long[, year := as.integer(sub("A00AA", "", VAR))]
pop_long[, geoid := paste0(
  sprintf("%02d", as.integer(STATEFP)),
  sprintf("%03d", as.integer(COUNTYFP))
)]
pop_long[, geoid := fmt_geoid5(geoid)]
pop_long[, pop_year := safe_numeric(pop_year)]

pop_dt <- pop_long[
  year %in% years,
  .(geoid, year, pop_year)
]

pop_dt <- pop_dt[
  geoid %in% counties_sf$geoid &
    !is.na(pop_year) &
    is.finite(pop_year) &
    pop_year > 0
]

cat("Population rows:", nrow(pop_dt), "\n\n")

# ===================== 3) COMPUTE OR LOAD UNIT-YEAR EXPOSURE CONTRIBUTION =====================

if (reuse_intermediate && file.exists(unit_year_file)) {

  cat("Existing intermediate file found. Loading:\n", unit_year_file, "\n\n")

  unit_year_dt <- fread(unit_year_file)

} else {

  cat("Computing unit-level population-weighted exposure contributions...\n")

  unit_year_results <- list()

  for (yy in years) {

    cat("\n================ YEAR ", yy, " ================\n", sep = "")

    pm25_file <- file.path(pm25_dir, sprintf("grids_pm25_byunit_%d.fst", yy))

    if (!file.exists(pm25_file)) {
      warning("Missing by-unit PM2.5 file: ", pm25_file)
      next
    }

    grid_dt <- fst::read_fst(pm25_file, as.data.table = TRUE)

    unit_cols <- get_unit_cols(grid_dt)

    if (length(unit_cols) == 0) {
      warning("No unit columns found in: ", pm25_file)
      next
    }

    cat("Grid cells:", nrow(grid_dt), "\n")
    cat("Unit columns:", length(unit_cols), "\n")

    pop_y <- pop_dt[year == yy]
    nat_pop <- sum(pop_y$pop_year, na.rm = TRUE)

    if (!is.finite(nat_pop) || nat_pop <= 0) {
      warning("Invalid national population for year ", yy)
      next
    }

    county_pop_vec <- pop_y$pop_year
    names(county_pop_vec) <- pop_y$geoid

    batch_starts <- seq(1, length(unit_cols), by = batch_size)
    unit_contrib_list <- vector("list", length(batch_starts))

    for (b in seq_along(batch_starts)) {

      start_i <- batch_starts[b]
      end_i <- min(start_i + batch_size - 1, length(unit_cols))
      batch_cols <- unit_cols[start_i:end_i]

      cat(
        "  Batch ", b, "/", length(batch_starts),
        " | units ", start_i, "-", end_i, "\n",
        sep = ""
      )

      county_unit_mat <- weighted_county_means_for_batch(
        grid_dt = grid_dt,
        unit_cols_batch = batch_cols,
        counties_sf = counties_sf,
        p4s = p4s
      )

      common_counties <- intersect(rownames(county_unit_mat), names(county_pop_vec))

      county_unit_mat <- county_unit_mat[common_counties, , drop = FALSE]
      pop_vec <- county_pop_vec[common_counties]

      weighted_sum <- colSums(county_unit_mat * pop_vec, na.rm = TRUE)
      contrib_nat <- weighted_sum / nat_pop

      unit_contrib_list[[b]] <- data.table(
        year = yy,
        uID = names(contrib_nat),
        exposure_contrib = as.numeric(contrib_nat)
      )

      rm(county_unit_mat)
      gc()
    }

    unit_year <- rbindlist(unit_contrib_list, use.names = TRUE, fill = TRUE)
    unit_year[!is.finite(exposure_contrib) | is.na(exposure_contrib), exposure_contrib := 0]

    unit_year_results[[as.character(yy)]] <- unit_year

    cat("Finished year ", yy, ". Unit rows: ", nrow(unit_year), "\n", sep = "")

    rm(grid_dt)
    gc()
  }

  unit_year_dt <- rbindlist(unit_year_results, use.names = TRUE, fill = TRUE)

  if (nrow(unit_year_dt) == 0) {
    stop("No unit-level exposure contribution results were generated.")
  }

  fwrite(unit_year_dt, unit_year_file)

  cat("Saved intermediate unit-year exposure file:\n", unit_year_file, "\n\n")
}

required_unit_cols <- c("year", "uID", "exposure_contrib")
missing_unit_cols <- setdiff(required_unit_cols, names(unit_year_dt))

if (length(missing_unit_cols) > 0) {
  stop(
    "Intermediate file is missing required columns:\n",
    paste(missing_unit_cols, collapse = ", ")
  )
}

unit_year_dt[, year := as.integer(year)]
unit_year_dt[, uID := as.character(uID)]
unit_year_dt[, exposure_contrib := safe_numeric(exposure_contrib)]
unit_year_dt[!is.finite(exposure_contrib) | is.na(exposure_contrib), exposure_contrib := 0]

cat("Loaded unit-year rows:", nrow(unit_year_dt), "\n")
cat("Unique years:", paste(sort(unique(unit_year_dt$year)), collapse = ", "), "\n")
cat("Unique units:", uniqueN(unit_year_dt$uID), "\n\n")

# ===================== 4) LOAD SO2 EMISSIONS =====================

cat("Loading SO2 emissions...\n")

so2_dt <- fread(so2_file)

required_so2_cols <- c(
  "Year",
  "Plant_ID",
  "Generator_ID",
  "Plant_Name",
  "Plant_State",
  "County",
  "Longitude",
  "Latitude",
  "SO2_tons"
)

missing_so2_cols <- setdiff(required_so2_cols, names(so2_dt))

if (length(missing_so2_cols) > 0) {
  stop(
    "SO2 file is missing required columns:\n",
    paste(missing_so2_cols, collapse = ", "),
    "\nAvailable columns are:\n",
    paste(names(so2_dt), collapse = ", ")
  )
}

so2_dt[, Year := as.integer(Year)]
so2_dt[, Plant_ID := as.character(Plant_ID)]
so2_dt[, Generator_ID := as.character(Generator_ID)]
so2_dt[, Plant_Name := as.character(Plant_Name)]
so2_dt[, Plant_State := as.character(Plant_State)]
so2_dt[, County := as.character(County)]
so2_dt[, Longitude := safe_numeric(Longitude)]
so2_dt[, Latitude := safe_numeric(Latitude)]
so2_dt[, SO2_tons := safe_numeric(SO2_tons)]

# Construct unit IDs consistent with PM2.5 by-unit column names.
so2_dt[, uID := make_uid(Plant_ID, Generator_ID)]

so2_unit_year <- so2_dt[
  ,
  .(
    Plant_ID = first_nonmissing_chr(Plant_ID),
    Plant_Name = first_nonmissing_chr(Plant_Name),
    Plant_State = first_nonmissing_chr(Plant_State),
    County = first_nonmissing_chr(County),
    Longitude = first_nonmissing_num(Longitude),
    Latitude = first_nonmissing_num(Latitude),
    SO2_tons = sum(SO2_tons, na.rm = TRUE)
  ),
  by = .(year = Year, uID)
]

cat("SO2 unit-year rows:", nrow(so2_unit_year), "\n")
cat("SO2 unique units:", uniqueN(so2_unit_year$uID), "\n\n")

# ===================== 5) MERGE EXPOSURE CONTRIBUTION WITH SO2 =====================

cat("Merging exposure contribution with SO2 emissions...\n")

unit_year_merged <- merge(
  unit_year_dt,
  so2_unit_year,
  by = c("year", "uID"),
  all.x = TRUE
)

match_rate <- mean(!is.na(unit_year_merged$Plant_ID)) * 100
cat("Matched exposure units to SO2 metadata:", round(match_rate, 2), "%\n\n")

unit_year_merged[is.na(SO2_tons) | !is.finite(SO2_tons), SO2_tons := 0]

# Recover Plant_ID from uID if SO2 metadata is unavailable.
unit_year_merged[
  is.na(Plant_ID) | Plant_ID == "",
  Plant_ID := sub("^X([^.]+)\\..*$", "\\1", uID)
]

unit_year_merged[
  is.na(Plant_Name) | Plant_Name == "",
  Plant_Name := paste0("Plant ", Plant_ID)
]

unit_year_merged[
  is.na(Plant_State) | Plant_State == "",
  Plant_State := NA_character_
]

# ===================== 6) AGGREGATE TO FACILITY-YEAR =====================

cat("Aggregating unit-level contributions to facility-year...\n")

facility_year <- unit_year_merged[
  ,
  .(
    exposure_contrib = sum(exposure_contrib, na.rm = TRUE),
    SO2_tons = sum(SO2_tons, na.rm = TRUE),
    Plant_Name = first_nonmissing_chr(Plant_Name),
    Plant_State = first_nonmissing_chr(Plant_State),
    County = first_nonmissing_chr(County),
    Longitude = first_nonmissing_num(Longitude),
    Latitude = first_nonmissing_num(Latitude),
    n_units = uniqueN(uID)
  ),
  by = .(year, Plant_ID)
]

facility_year[!is.finite(exposure_contrib) | is.na(exposure_contrib), exposure_contrib := 0]
facility_year[!is.finite(SO2_tons) | is.na(SO2_tons), SO2_tons := 0]

fwrite(
  facility_year,
  file.path(out_dir, "facility_year_contribution_1940_1990.csv")
)

cat("Facility-year rows:", nrow(facility_year), "\n\n")

# ===================== 7) AGGREGATE TO CUMULATIVE FACILITY CONTRIBUTION =====================

cat("Aggregating facility-year contributions to cumulative facility contribution...\n")

facility_cum <- facility_year[
  ,
  .(
    cumulative_exposure_contrib = sum(exposure_contrib, na.rm = TRUE),
    cumulative_SO2_tons = sum(SO2_tons, na.rm = TRUE),
    Plant_Name = first_nonmissing_chr(Plant_Name),
    Plant_State = first_nonmissing_chr(Plant_State),
    County = first_nonmissing_chr(County),
    Longitude = first_nonmissing_num(Longitude),
    Latitude = first_nonmissing_num(Latitude),
    n_units = max(n_units, na.rm = TRUE)
  ),
  by = Plant_ID
]

facility_cum <- facility_cum[
  cumulative_exposure_contrib > 0 |
    cumulative_SO2_tons > 0
]

fwrite(
  facility_cum,
  file.path(out_dir, "facility_cumulative_contribution_1940_1990.csv")
)

cat("Facility cumulative rows:", nrow(facility_cum), "\n\n")

# ===================== 8) FIG. 4A: FACILITY EXPOSURE CONCENTRATION CURVE =====================

cat("Creating Fig. 4a...\n")

pareto_dt <- copy(facility_cum)
pareto_dt <- pareto_dt[order(-cumulative_exposure_contrib)]

total_exposure <- sum(pareto_dt$cumulative_exposure_contrib, na.rm = TRUE)

if (!is.finite(total_exposure) || total_exposure <= 0) {
  stop("Total cumulative exposure contribution is zero or invalid.")
}

pareto_dt[, facility_rank := seq_len(.N)]
pareto_dt[, facility_fraction := facility_rank / .N]
pareto_dt[, cumulative_exposure_share := cumsum(cumulative_exposure_contrib) / total_exposure]
pareto_dt[, facility_percent := facility_fraction * 100]
pareto_dt[, exposure_percent := cumulative_exposure_share * 100]

summary_cutoffs <- data.table(
  cutoff = target_cutoffs
)

summary_cutoffs[, top_n := pmax(1, ceiling(cutoff * nrow(pareto_dt)))]
summary_cutoffs[, exposure_share := pareto_dt$cumulative_exposure_share[top_n]]
summary_cutoffs[, exposure_percent := exposure_share * 100]
summary_cutoffs[, total_facilities := nrow(pareto_dt)]

fwrite(
  summary_cutoffs,
  file.path(out_dir, "FIG_4A_top_facility_exposure_shares.csv")
)

p_fig4a <- ggplot(pareto_dt, aes(x = facility_percent, y = exposure_percent)) +
  geom_line(linewidth = 1.1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.5) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    x = "Facilities ranked by exposure contribution (%)",
    y = "Cumulative exposure contribution (%)"
  )

ggsave(
  file.path(out_dir, "FIG_4A_facility_pareto_curve.pdf"),
  p_fig4a,
  width = 5.8,
  height = 4.4,
  dpi = 300
)

# ===================== 9) FIG. 4B: TOP FACILITY CONTRIBUTIONS =====================

cat("Creating Fig. 4b...\n")

top_facilities <- facility_cum[
  order(-cumulative_exposure_contrib)
][
  1:min(top_n_facilities, .N)
]

top_facilities[, label := ifelse(
  is.na(Plant_State) | Plant_State == "",
  paste0(Plant_Name),
  paste0(Plant_Name, ", ", Plant_State)
)]

top_facilities[, label := factor(label, levels = rev(label))]

p_fig4b <- ggplot(
  top_facilities,
  aes(x = label, y = cumulative_exposure_contrib)
) +
  geom_col(width = 0.75) +
  coord_flip() +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    x = NULL,
    y = expression("Cumulative contribution to national PWE"~(mu*g/m^3))
  )

ggsave(
  file.path(out_dir, "FIG_4B_top_facility_contributions.pdf"),
  p_fig4b,
  width = 7.2,
  height = 5.4,
  dpi = 300
)

fwrite(
  top_facilities,
  file.path(out_dir, "FIG_4B_top_facility_contributions.csv")
)

# ===================== 10) FIG. 4C: TARGETING EFFICIENCY =====================

cat("Creating Fig. 4c...\n")

target_base <- copy(facility_cum)
target_base <- target_base[
  cumulative_exposure_contrib > 0 |
    cumulative_SO2_tons > 0
]

total_exp <- sum(target_base$cumulative_exposure_contrib, na.rm = TRUE)

if (!is.finite(total_exp) || total_exp <= 0) {
  stop("Total exposure for targeting analysis is zero or invalid.")
}

n_fac <- nrow(target_base)

targeting_results <- list()

for (cc in target_cutoffs) {

  n_remove <- pmax(1, ceiling(cc * n_fac))

  exposure_ranked_ids <- target_base[
    order(-cumulative_exposure_contrib)
  ][
    1:n_remove,
    Plant_ID
  ]

  so2_ranked_ids <- target_base[
    order(-cumulative_SO2_tons)
  ][
    1:n_remove,
    Plant_ID
  ]

  exposure_reduction_exposure_ranked <- target_base[
    Plant_ID %in% exposure_ranked_ids,
    sum(cumulative_exposure_contrib, na.rm = TRUE)
  ] / total_exp

  exposure_reduction_so2_ranked <- target_base[
    Plant_ID %in% so2_ranked_ids,
    sum(cumulative_exposure_contrib, na.rm = TRUE)
  ] / total_exp

  targeting_results[[length(targeting_results) + 1]] <- data.table(
    cutoff_fraction = cc,
    cutoff_percent = cc * 100,
    strategy = "Exposure-based ranking",
    exposure_reduction_percent = exposure_reduction_exposure_ranked * 100,
    n_removed = n_remove
  )

  targeting_results[[length(targeting_results) + 1]] <- data.table(
    cutoff_fraction = cc,
    cutoff_percent = cc * 100,
    strategy = "SO2 emissions-based ranking",
    exposure_reduction_percent = exposure_reduction_so2_ranked * 100,
    n_removed = n_remove
  )
}

targeting_dt <- rbindlist(targeting_results)

fwrite(
  targeting_dt,
  file.path(out_dir, "FIG_4C_targeting_efficiency.csv")
)

# ===================== 10.1) MANUSCRIPT SUMMARY TABLES =====================

# This table gives the facility counts (n) associated with each percentage cutoff
# and the cumulative exposure share captured by the top-ranked facilities.
facility_count_summary <- summary_cutoffs[
  ,
  .(
    cutoff_fraction = cutoff,
    cutoff_percent = cutoff * 100,
    total_facilities,
    n_facilities = top_n,
    cumulative_exposure_percent = exposure_percent
  )
][order(cutoff_fraction)]

# Convert the targeting results to wide format so that the exposure-based and
# SO2-based strategies are shown side by side for manuscript writing.
targeting_wide <- dcast(
  targeting_dt,
  cutoff_fraction + cutoff_percent + n_removed ~ strategy,
  value.var = "exposure_reduction_percent"
)

setnames(
  targeting_wide,
  old = c("Exposure-based ranking", "SO2 emissions-based ranking"),
  new = c("exposure_based_reduction_percent", "so2_based_reduction_percent"),
  skip_absent = TRUE
)

manuscript_summary <- merge(
  facility_count_summary,
  targeting_wide[
    ,
    .(
      cutoff_fraction,
      cutoff_percent,
      n_removed,
      exposure_based_reduction_percent,
      so2_based_reduction_percent
    )
  ],
  by = c("cutoff_fraction", "cutoff_percent"),
  all.x = TRUE
)

# Rounded version for direct manuscript use.
manuscript_summary_rounded <- copy(manuscript_summary)
num_cols_round <- c(
  "cumulative_exposure_percent",
  "exposure_based_reduction_percent",
  "so2_based_reduction_percent"
)
manuscript_summary_rounded[
  ,
  (num_cols_round) := lapply(.SD, function(x) round(x, 1)),
  .SDcols = num_cols_round
]

fwrite(
  facility_count_summary,
  file.path(out_dir, "facility_count_summary_by_cutoff.csv")
)

fwrite(
  manuscript_summary_rounded,
  file.path(out_dir, "Facility_targeting_summary_for_manuscript.csv")
)

cat("\n========================================\n")
cat("Facility ranking summary (1940-1990 cumulative)\n")
cat("========================================\n")
cat("Total facilities:", unique(manuscript_summary_rounded$total_facilities), "\n\n")
print(
  manuscript_summary_rounded[
    ,
    .(
      cutoff_percent,
      n_facilities,
      cumulative_exposure_percent,
      exposure_based_reduction_percent,
      so2_based_reduction_percent
    )
  ]
)
cat("\nSaved manuscript summary table:\n")
cat(file.path(out_dir, "Facility_targeting_summary_for_manuscript.csv"), "\n\n")

p_fig4c <- ggplot(
  targeting_dt,
  aes(
    x = cutoff_percent,
    y = exposure_reduction_percent,
    linetype = strategy,
    shape = strategy
  )
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.4) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    x = "Facilities targeted (%)",
    y = "Cumulative exposure reduced (%)"
  )

ggsave(
  file.path(out_dir, "FIG_4C_targeting_efficiency.pdf"),
  p_fig4c,
  width = 6.2,
  height = 4.6,
  dpi = 300
)

# ===================== 11) PRINT KEY RESULTS =====================

cat("\nKey exposure concentration results:\n")
print(summary_cutoffs)

cat("\nTargeting results:\n")
print(targeting_dt)

cat("\nAll outputs saved in:\n", out_dir, "\n")
