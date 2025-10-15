# =======================================================================
# HyADS → PM2.5 conversion (annual, decadal years 1940–1990)
# Boss-style R script for Slurm array or single-run
# -----------------------------------------------------------------------
# srun -p normal --mem 50g -t 0-04:00 -c 4 -N 1 --pty /bin/bash
# module load gnu10/10.3.0
# module load openmpi
# module load netcdf-c
# module load r/4.1.2-dx
# source('~/repos/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R')
# load('/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/hyads_to_cmaq_models20230516.RData')
# =======================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
  library(sf)
  library(USAboundaries)
})

# ----------------------- Years to run -----------------------
years.run <- c(1940, 1950, 1960, 1970, 1980, 1990)

# Pick by Slurm array (0-based index example: --array=0-5)
arr <- Sys.getenv("SLURM_ARRAY_TASK_ID")
array_num <- suppressWarnings(as.numeric(arr))
if (!is.na(array_num)) {
  # convert 0-based to 1-based if you submit --array=0-5
  if (array_num == 0) array_num <- 1
  years.run <- years.run[array_num]
}
cat("[INFO] Years selected:", paste(years.run, collapse = ", "), "\n")

# ----------------------- CRS / Projection -------------------
p4s_hyads <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

# ----------------------- Functions & Models -----------------
FUNC_PATH  <- "~/repos/HyADS_to_pm25/RCode/hyads_to_pm25_functions.R"
MODELS_RDA <- "/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/hyads_to_cmaq_models20230516.RData"
source(FUNC_PATH)
load(MODELS_RDA)   # defines e.g., preds.ann.hyads06w05 and model objects

# ----------------------- IO: Inputs / Outputs ---------------
# INPUT: per-ton HyADS (schema aligned to 2022) produced in STEP F
hyads.dir <- "/scratch/xshan2/R_Code/disperseR/byunit_fst_perton"
#   byunit  prefix: grids_exposures_byunit_<YEAR>.fst
#   total   prefix: grids_exposures_total_<YEAR>.fst

# OUTPUT: PM2.5 grids
out.dir   <- "/scratch/xshan2/R_Code/disperseR/exp25/grids_model.lm.cv_single_poly"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------- USA mask (CONUS) -------------------
us_states.names <- state.abb[!(state.abb %in% c("HI", "AK"))]
us_states <- st_transform(USAboundaries::us_states(), p4s_hyads)
mask.usa  <- sf::as_Spatial(us_states)[us_states$state_abbr %in% us_states.names,]

# ----------------------- Model Settings ---------------------
model.dataset <- preds.ann.hyads06w05
model.name    <- "model.lm.cv_single_poly"
name.x        <- "hyads"
met.dir       <- "/projects/HAQ_LAB/lhennem/data/disperseR/HyADS_to_pm25/met"

# ----------------------- BY-UNIT → PM2.5 --------------------
cat("[INFO] Converting BY-UNIT per-ton HyADS to PM2.5 ...\n")
lapply(
  years.run,
  hyads_to_pm25_unit,
  fstart        = file.path(hyads.dir, "grids_exposures_byunit_"),
  fstart.total  = file.path(hyads.dir, "grids_exposures_total_"),  # required by function signature
  fstart_out    = file.path(out.dir,   "grids_pm25_byunit_"),
  model.dataset = model.dataset,
  model.name    = model.name,
  name.x        = name.x,
  mask.use      = mask.usa,
  p4s           = p4s_hyads,
  met.dest      = met.dir
)

# ----------------------- TOTAL → PM2.5 ----------------------
cat("[INFO] Converting TOTAL per-ton HyADS to PM2.5 ...\n")
lapply(
  years.run,
  hyads_to_pm25_unit,
  fstart.total  = file.path(hyads.dir, "grids_exposures_total_"),
  fstart_out    = file.path(out.dir,   "grids_pm25_total_"),
  model.dataset = model.dataset,
  model.name    = model.name,
  name.x        = name.x,
  mask.use      = mask.usa,
  p4s           = p4s_hyads,
  met.dest      = met.dir,
  total         = TRUE
)

cat("[DONE] Outputs written to:\n  ", out.dir, "\n")
