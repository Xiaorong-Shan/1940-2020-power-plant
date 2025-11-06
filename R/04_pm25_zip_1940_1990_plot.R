#!/usr/bin/env Rscript
# ============================================================
# ZIP-level Coal PM2.5 maps & boxplot (1940–1990)
# - Inputs: zips_pm25_total_YYYY.fst (preferred) or zips_pm25_byunit_YYYY.fst
# - Color scale auto-set by global P99
# - Headless-safe: PNG via ragg (if installed) + PDF fallback
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(data.table)
  library(fst)
  library(ggplot2)
  library(USAboundaries)
  library(grid)    # unit()
})

# ---------------- PATHS (UPDATED) ----------------
# Where your ZIP-level files were written in the previous step:
OUT_DIR <- "/scratch/xshan2/R_Code/disperseR/main/output/pm25_decades_model.lm.cv_single_poly_proxy1999met/zips_exposures"

# Where to save figures:
FIG_DIR <- "/scratch/xshan2/R_Code/disperseR/main/output/pm25_decades_model.lm.cv_single_poly_proxy1999met/figures"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ZCTA geometry + crosswalk (unchanged)
ZCTA_DIR <- "/projects/HAQ_LAB/lhennem/disperseR/main/input/zcta_500k"
ZCTA_SHP <- file.path(ZCTA_DIR, "cb_2017_us_zcta510_500k.shp")

YEARS <- c(1940,1950,1960,1970,1980,1990)

# ---------------- ZIP GEOMETRY ----------------
p4s <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
cw <- disperseR::crosswalk

zips_sf <- st_read(ZCTA_SHP, quiet = TRUE)
if ("ZCTA5CE10" %in% names(zips_sf)) names(zips_sf)[names(zips_sf) == "ZCTA5CE10"] <- "ZCTA"
zips_sf <- merge(zips_sf, cw, by = "ZCTA", all = FALSE, allow.cartesian = TRUE)
zips_sf$ZIP <- formatC(zips_sf$ZIP, width = 5, format = "d", flag = "0")
zips_sf <- st_transform(zips_sf, p4s)

states <- USAboundaries::us_states()
states48 <- states[states$name %in% c(setdiff(state.name, c("Alaska","Hawaii")), "District of Columbia"), ]
states48 <- st_transform(states48, st_crs(zips_sf))

# ---------------- READ ONE YEAR (robust) ----------------
read_zip_value <- function(fpath) {
  dt <- read.fst(fpath, as.data.table = TRUE)

  if (!"ZIP" %in% names(dt)) {
    if ("zip" %in% names(dt)) setnames(dt, "zip", "ZIP")
    if ("ZCTA" %in% names(dt)) setnames(dt, "ZCTA", "ZIP")
  }
  stopifnot("ZIP" %in% names(dt))

  yr_chr <- sub(".*_(\\d{4})\\.fst$", "\\1", basename(fpath))
  yr <- suppressWarnings(as.integer(yr_chr))
  if (is.na(yr)) stop("Cannot parse year from filename: ", fpath)

  val_cols <- setdiff(names(dt), "ZIP")
  if ("total" %in% val_cols) {
    dt[, value := get("total")]
  } else {
    dt[, value := rowSums(.SD, na.rm = TRUE), .SDcols = val_cols]
  }
  dt[, year := yr]
  dt[, .(ZIP, year, value)]
}

# ---------------- COLLECT INPUTS ----------------
files_total <- file.path(OUT_DIR, sprintf("zips_pm25_total_%d.fst", YEARS))
if (all(file.exists(files_total))) {
  files <- files_total
  message("Using zips_pm25_total_YYYY.fst")
} else {
  files_byu <- file.path(OUT_DIR, sprintf("zips_pm25_byunit_%d.fst", YEARS))
  if (!all(file.exists(files_byu))) stop("Missing inputs in: ", OUT_DIR)
  files <- files_byu
  message("Using zips_pm25_byunit_YYYY.fst (will row-sum uIDs)")
}

# ---------------- LOAD & PREPARE ----------------
zips_dt <- rbindlist(lapply(files, read_zip_value))
zips_dt[, year := factor(year, levels = YEARS, labels = YEARS)]

# Global P99 for top of the scale
p99 <- as.numeric(quantile(zips_dt$value, 0.99, na.rm = TRUE))
if (!is.finite(p99) || p99 <= 0) p99 <- max(zips_dt$value, na.rm = TRUE)
if (!is.finite(p99) || p99 <= 0) p99 <- 1e-4
LIMS <- c(0, p99)

# Neat colorbar ticks
brks <- pretty(LIMS, n = 3)
brks <- brks[brks >= 0 & brks <= p99]
if (length(brks) < 3) brks <- c(0, p99/2, p99)
labs <- formatC(brks, format = "fg", digits = 2)
labs[length(labs)] <- paste0("≥", labs[length(labs)])

# Merge geometry
zips_dt_sf <- merge(zips_sf, zips_dt, by = "ZIP", all.y = TRUE)

# ---------------- PLOTS ----------------
spat.gg <- ggplot(zips_dt_sf) +
  geom_sf(aes(geometry = geometry, fill = value), color = NA) +
  scale_fill_gradient(
    low = "cornsilk", high = "black",
    limits = LIMS, breaks = brks, labels = labs,
    na.value = NA, oob = scales::squish,
    guide = guide_colourbar(
      barheight = unit(1.6, "cm"),
      barwidth  = unit(0.30, "cm"),
      ticks = TRUE
    )
  ) +
  geom_sf(data = states48, inherit.aes = FALSE,
          fill = NA, color = "grey90", linewidth = 0.05) +
  facet_wrap(~ year, nrow = 1) +
  labs(fill = expression(paste("Coal PM"[2.5], " (", mu, "g ", m^{-3}, ")"))) +
  theme_bw() +
  theme(
    axis.text  = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.text = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    legend.position = "right"
  )

hist.gg <- ggplot(as.data.table(zips_dt), aes(x = year, y = value, group = year)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.4) +
  labs(x = NULL, y = expression(paste("Coal PM"[2.5], " (", mu, "g ", m^{-3}, ")"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11))

# ---------------- SAVE (PNG + PDF) ----------------
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png(file.path(FIG_DIR, "coal_pm25_trends_1940_1990.png"),
                width = 14, height = 4.5, units = "in", res = 300)
  print(spat.gg); dev.off()

  ragg::agg_png(file.path(FIG_DIR, "coal_pm25_box_1940_1990.png"),
                width = 6, height = 4, units = "in", res = 300)
  print(hist.gg); dev.off()
} else {
  message("ragg not installed; PNGs skipped (PDFs will still be written).")
}

pdf(file.path(FIG_DIR, "coal_pm25_trends_1940_1990.pdf"), width = 14, height = 4.5)
print(spat.gg); dev.off()

pdf(file.path(FIG_DIR, "coal_pm25_box_1940_1990.pdf"), width = 6, height = 4)
print(hist.gg); dev.off()

cat("Saved figures to:", FIG_DIR, "\n")
