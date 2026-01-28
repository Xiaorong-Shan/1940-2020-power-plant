# ============================================================
# AQS (1980/1990) vs HyADS (decade grids) evaluation + PDF plots
# - PDF only (hopper-friendly)
# - No plot titles
# - Clear figure naming (Main vs SI)
# Main Figure 1:
#   - 1990 PM2.5 STP (AQS) vs HyADS modeled decade-mean PM2.5
# Supplementary Figures:
#   - S1: 1990 PM10
#   - S2: 1990 TSP
#   - S3: 1990 Sulfate(TSP)
#   - S4: 1980 TSP
#   - S5: 1980 Sulfate(TSP)
#
# Outputs (in out_dir):
#   - evaluation_summary.csv
#   - Fig_Main1_PM25STP_1990_scatter.pdf
#   - Fig_S1_PM10_1990_scatter.pdf
#   - Fig_S2_TSP_1990_scatter.pdf
#   - Fig_S3_SulfateTSP_1990_scatter.pdf
#   - Fig_S4_TSP_1980_scatter.pdf
#   - Fig_S5_SulfateTSP_1980_scatter.pdf
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
  library(sf)
  library(ggplot2)
})

# -----------------------------
# User inputs
# -----------------------------
p4s <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

eval_dir <- "/scratch/xshan2/R_Code/disperseR/main/output/pm25_decades_model.lm.cv_single_poly_proxy1999met/evaluation"
aqs_1980 <- file.path(eval_dir, "annual_conc_by_monitor_1980.csv")
aqs_1990 <- file.path(eval_dir, "annual_conc_by_monitor_1990.csv")

hyads_dir <- "/scratch/xshan2/R_Code/disperseR/main/output/pm25_decades_model.lm.cv_single_poly_proxy1999met"
hyads_1980 <- file.path(hyads_dir, "grids_pm25_total_1980.fst")
hyads_1990 <- file.path(hyads_dir, "grids_pm25_total_1990.fst")

out_dir <- file.path(eval_dir, "aqs_hyads_eval_out_pdfonly")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(aqs_1980), file.exists(aqs_1990))
stopifnot(file.exists(hyads_1980), file.exists(hyads_1990))

# -----------------------------
# Helpers
# -----------------------------
read_aqs_annual <- function(f) {
  dt <- fread(f)
  setnames(dt, gsub(" ", "_", names(dt)))
  req <- c("State_Code","County_Code","Site_Num","Latitude","Longitude",
           "Parameter_Name","Parameter_Code","Arithmetic_Mean","Year")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) stop("AQS file missing columns: ", paste(miss, collapse=", "))
  dt
}

read_hyads_grid <- function(f, p4s) {
  dt <- read_fst(f, as.data.table = TRUE)
  req <- c("x","y","vals.out")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) stop("HyADS fst missing columns: ", paste(miss, collapse=", "))
  st_as_sf(dt, coords = c("x","y"), crs = st_crs(p4s), remove = FALSE)
}

to_pm_sf <- function(dt) {
  pm <- dt[, .(
    site_id = paste(State_Code, County_Code, Site_Num, sep = "_"),
    lon = as.numeric(Longitude),
    lat = as.numeric(Latitude),
    obs = as.numeric(Arithmetic_Mean),
    year = as.integer(Year),
    Parameter_Name,
    Parameter_Code
  )]
  pm <- pm[is.finite(lon) & is.finite(lat) & is.finite(obs)]
  st_as_sf(pm, coords = c("lon","lat"), crs = 4326)
}

match_nearest <- function(obs_sf_ll, hyads_sf, p4s) {
  obs_sf <- st_transform(obs_sf_ll, st_crs(p4s))
  idx <- st_nearest_feature(obs_sf, hyads_sf)
  obs_sf$mod <- hyads_sf$vals.out[idx]
  obs_sf
}

metrics <- function(obs, mod) {
  ok <- is.finite(obs) & is.finite(mod)
  obs <- obs[ok]; mod <- mod[ok]
  if (length(obs) < 10) return(list(n=length(obs), pearson=NA_real_, spearman=NA_real_, r2=NA_real_))
  r  <- suppressWarnings(cor(obs, mod, method="pearson"))
  rs <- suppressWarnings(cor(obs, mod, method="spearman"))
  list(n=length(obs), pearson=r, spearman=rs, r2=r^2)
}

plot_scatter_pdf <- function(df, xlab, ylab, pdf_name) {
  m <- metrics(df$obs, df$mod)
  ann <- sprintf("n=%d\nPearson r=%.3f\nSpearman \u03C1=%.3f\nR\u00B2=%.3f",
                 m$n, m$pearson, m$spearman, m$r2)

  p <- ggplot(df, aes(x = mod, y = obs)) +
    geom_point(alpha = 0.35, size = 1) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = xlab, y = ylab) +  # ✅ no title
    annotate("text", x = Inf, y = -Inf,
             hjust = 1.1, vjust = -0.2,
             label = ann, size = 4) +
    theme_bw(base_size = 12)

  ggsave(file.path(out_dir, pdf_name), p, width = 6.8, height = 5.2)
  invisible(p)
}

eval_one <- function(sfobj, label, decade) {
  df <- st_drop_geometry(sfobj)
  m  <- metrics(df$obs, df$mod)
  data.table(
    pair = label,
    decade = decade,
    n_sites = uniqueN(df$site_id),
    n_obs = m$n,
    pearson_r = m$pearson,
    spearman_rho = m$spearman,
    r2 = m$r2
  )
}

# -----------------------------
# Load data
# -----------------------------
dt80 <- read_aqs_annual(aqs_1980)
dt90 <- read_aqs_annual(aqs_1990)

hy80_sf <- read_hyads_grid(hyads_1980, p4s)
hy90_sf <- read_hyads_grid(hyads_1990, p4s)

# -----------------------------
# Pollutant filters (AQS Parameter_Name)
# -----------------------------
pm25_stp_name <- "PM2.5 STP"
pm10_name     <- "PM10 Total 0-10um STP"
tsp_name      <- "Suspended particulate (TSP)"
is_sulf_tsp   <- function(x) grepl("^Sulfate \\(TSP\\)", x)

# -----------------------------
# Subset AQS by pollutant
# -----------------------------
# 1990
dt90_pm25 <- dt90[Parameter_Name == pm25_stp_name]
dt90_pm10 <- dt90[Parameter_Name == pm10_name]
dt90_tsp  <- dt90[Parameter_Name == tsp_name]
dt90_sulf <- dt90[is_sulf_tsp(Parameter_Name)]

# 1980
dt80_tsp  <- dt80[Parameter_Name == tsp_name]
dt80_sulf <- dt80[is_sulf_tsp(Parameter_Name)]

# Fail loudly if main figure has no data
if (nrow(dt90_pm25) == 0) {
  stop("No rows found for 1990 'PM2.5 STP' in AQS file. Check Parameter_Name values.")
}

# -----------------------------
# Convert to sf + match to HyADS (nearest grid cell)
# -----------------------------
pm90_pm25_sf <- match_nearest(to_pm_sf(dt90_pm25), hy90_sf, p4s)  # Main
pm90_pm10_sf <- match_nearest(to_pm_sf(dt90_pm10), hy90_sf, p4s)  # SI
pm90_tsp_sf  <- match_nearest(to_pm_sf(dt90_tsp),  hy90_sf, p4s)  # SI
pm90_sulf_sf <- match_nearest(to_pm_sf(dt90_sulf), hy90_sf, p4s)  # SI

pm80_tsp_sf  <- match_nearest(to_pm_sf(dt80_tsp),  hy80_sf, p4s)  # SI
pm80_sulf_sf <- match_nearest(to_pm_sf(dt80_sulf), hy80_sf, p4s)  # SI

# -----------------------------
# Metrics table
# -----------------------------
summary_dt <- rbindlist(list(
  eval_one(pm90_pm25_sf, "PM2.5 STP (AQS) vs HyADS PM2.5", 1990),
  eval_one(pm90_pm10_sf, "PM10 (AQS) vs HyADS PM2.5",      1990),
  eval_one(pm90_tsp_sf,  "TSP (AQS) vs HyADS PM2.5",       1990),
  eval_one(pm90_sulf_sf, "Sulfate(TSP) (AQS) vs HyADS PM2.5", 1990),
  eval_one(pm80_tsp_sf,  "TSP (AQS) vs HyADS PM2.5",       1980),
  eval_one(pm80_sulf_sf, "Sulfate(TSP) (AQS) vs HyADS PM2.5", 1980)
), fill = TRUE)

fwrite(summary_dt, file.path(out_dir, "evaluation_summary.csv"))
print(summary_dt)

# -----------------------------
# Axis labels (PM2.5 subscript)
# -----------------------------
xlab <- expression(paste("Modeled decade-mean ", PM[2.5], " (HyADS)"))
ylab <- expression(paste("Observed annual mean concentration (AQS, ", mu, "g/", m^3, ")"))

# -----------------------------
# PDF plots (NO titles)
# Naming with Main vs SI labels
# -----------------------------
# Main Figure 1
plot_scatter_pdf(
  st_drop_geometry(pm90_pm25_sf)[, c("obs","mod")],
  xlab, ylab,
  "Fig_Main1_PM25STP_1990_scatter.pdf"
)

# Supplementary
plot_scatter_pdf(
  st_drop_geometry(pm90_pm10_sf)[, c("obs","mod")],
  xlab, ylab,
  "Fig_S1_PM10_1990_scatter.pdf"
)

plot_scatter_pdf(
  st_drop_geometry(pm90_tsp_sf)[, c("obs","mod")],
  xlab, ylab,
  "Fig_S2_TSP_1990_scatter.pdf"
)

plot_scatter_pdf(
  st_drop_geometry(pm90_sulf_sf)[, c("obs","mod")],
  xlab, ylab,
  "Fig_S3_SulfateTSP_1990_scatter.pdf"
)

plot_scatter_pdf(
  st_drop_geometry(pm80_tsp_sf)[, c("obs","mod")],
  xlab, ylab,
  "Fig_S4_TSP_1980_scatter.pdf"
)

plot_scatter_pdf(
  st_drop_geometry(pm80_sulf_sf)[, c("obs","mod")],
  xlab, ylab,
  "Fig_S5_SulfateTSP_1980_scatter.pdf"
)

cat("\nDone. PDFs + summary table saved to:\n", out_dir, "\n")
