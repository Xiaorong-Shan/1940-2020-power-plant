# ===============================================================
# 0. Load required libraries
# ===============================================================
library(data.table)
library(raster)
library(sf)
library(magrittr)
library(ggplot2)
library(viridis)
library(ncdf4)
library(dplyr)
library(Metrics)
library(randomForest)
library(xgboost)
library(tidyverse)
library(readxl)
library(fst)

# ===============================================================
# 1. Load CAMPD facility attribute and emission datasets
# ===============================================================
facility_info <- read.csv("/Users/xshan2/OneDrive - George Mason University - O365 Production/GMU_PhD/01_Research/02_2019fall_ADRD/prelim_analysis/PowerPlant/EF_pp_paper_code/data/facility-attributes-1995-2020-CAMPD.csv")
facility_emission <- read.csv("/Users/xshan2/OneDrive - George Mason University - O365 Production/GMU_PhD/01_Research/02_2019fall_ADRD/prelim_analysis/PowerPlant/EF_pp_paper_code/data/annual-emissions-facility-aggregation-1995-2020-CAMPD.csv")

# Find shared columns for merging
common_columns <- intersect(names(facility_info), names(facility_emission))

# Merge the two datasets
facility_info_emiss.all <- merge(facility_info, facility_emission, by = common_columns)

# ===============================================================
# 2. Data cleaning and preprocessing
# ===============================================================
facility_info_emiss.all$Commercial.Operation.Date <- as.Date(facility_info_emiss.all$Commercial.Operation.Date)

# Keep only facilities that started operating after 1940
subset_after_1940 <- facility_info_emiss.all %>%
  filter(Commercial.Operation.Date > as.Date("1939-12-31"))

# ===============================================================
# 3. Compute SO2 emission factors and additional variables
# ===============================================================
EF_SO2 <- subset_after_1940 %>%
  mutate(
    convert_fuel = case_when(
      grepl("coal", Primary.Fuel.Type, ignore.case = TRUE) ~ "Coal",
      grepl("gas", Primary.Fuel.Type, ignore.case = TRUE) ~ "Gas",
      grepl("oil|diesel|petroleum", Primary.Fuel.Type, ignore.case = TRUE) ~ "Oil",
      grepl("wood|biomass", Primary.Fuel.Type, ignore.case = TRUE) ~ "Biomass",
      TRUE ~ "Other"
    ),
    Unit.Type = gsub("\\s*\\(.*?\\)", "", Unit.Type),
    scrubber = ifelse(grepl("scrubber|FGD", SO2.Controls, ignore.case = TRUE), 1, 0),
    log_Gen_kWh = log10(Gross.Load..MWh. * 1000 + 1),
    SO2_grams = SO2.Mass..short.tons. * 907185,
    EF = SO2_grams / (Gross.Load..MWh. * 1000),
    log_EF = log10(EF + 1)
  ) %>%
  filter(EF < 10)

# ===============================================================
# 4. Split data into training and testing sets
# ===============================================================
df <- EF_SO2
df$convert_fuel <- as.factor(df$convert_fuel)
df$Unit.Type <- as.factor(df$Unit.Type)

set.seed(123)
train_index <- sample(1:nrow(df), 0.8 * nrow(df))
train_data <- df[train_index, ]
test_data  <- df[-train_index, ]

test_data$convert_fuel <- factor(test_data$convert_fuel, levels = levels(train_data$convert_fuel))
test_data$Unit.Type <- factor(test_data$Unit.Type, levels = levels(train_data$Unit.Type))

# ===============================================================
# 5. Train XGBoost model to predict log(EF)
# ===============================================================
X_train <- model.matrix(log_EF ~ convert_fuel + Unit.Type + Year + scrubber + log_Gen_kWh - 1, data = train_data)
X_test  <- model.matrix(log_EF ~ convert_fuel + Unit.Type + Year + scrubber + log_Gen_kWh - 1, data = test_data)
y_train <- train_data$log_EF
y_test  <- test_data$log_EF

dtrain <- xgb.DMatrix(data = X_train, label = y_train)
dtest  <- xgb.DMatrix(data = X_test, label = y_test)

params <- list(
  objective = "reg:squarederror",
  max_depth = 6,
  eta = 0.05,
  subsample = 0.7,
  colsample_bytree = 0.7,
  eval_metric = "rmse"
)

xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 300,
  early_stopping_rounds = 10,
  watchlist = list(train = dtrain, test = dtest),
  print_every_n = 50
)

# ===============================================================
# 6. Evaluate model performance
# ===============================================================
pred_log <- predict(xgb_model, newdata = dtest)
pred_EF  <- exp(pred_log)

mae  <- mean(abs(exp(y_test) - pred_EF))
rmse <- sqrt(mean((exp(y_test) - pred_EF)^2))
r2   <- 1 - sum((exp(y_test) - pred_EF)^2) / sum((exp(y_test) - mean(exp(y_test)))^2)

cat("📈 Model Performance (Exp[Log-EF]):\n")
cat("  MAE  =", round(mae, 5), "\n")
cat("  RMSE =", round(rmse, 5), "\n")
cat("  R²   =", round(r2, 5), "\n")

# ===============================================================
# 7. Apply model to historical generator data (EIA 1940–1990)
# ===============================================================
generator2010ret <- read_excel('/Users/xshan2/Library/CloudStorage/OneDrive-GeorgeMasonUniversity-O365Production/GMU_PhD/01_Research/02_2019fall_ADRD/prelim_analysis/PowerPlant/december_generator2022.xlsx',
                               sheet = 3, skip = 2) %>% data.table
generator2010act <- read_excel('/Users/xshan2/Library/CloudStorage/OneDrive-GeorgeMasonUniversity-O365Production/GMU_PhD/01_Research/02_2019fall_ADRD/prelim_analysis/PowerPlant/december_generator2022.xlsx',
                               sheet = 1, skip = 2) %>% data.table

if (!("Retirement Year" %in% names(generator2010act))) {
  generator2010act[, `Retirement Year` := NA_integer_]
}

cols_keep <- c("Entity ID", "Plant ID", "Plant Name", "Plant State", "County",
               "Nameplate Capacity (MW)", "Operating Year", "Retirement Year",
               "Generator ID", "Energy Source Code", "Latitude", "Longitude")

keep_act <- intersect(cols_keep, names(generator2010act))
keep_ret <- intersect(cols_keep, names(generator2010ret))

gen_all <- rbindlist(
  list(generator2010act[, ..keep_act], generator2010ret[, ..keep_ret]),
  use.names = TRUE, fill = TRUE
)[, ..cols_keep] |> unique()

# ===============================================================
# 8. Fuel code mapping to main categories
# ===============================================================
map_fuel <- function(code){
  fcase(
    code %in% c("BIT","SUB","LIG","RC","WC"), "Coal",
    code %in% c("DFO","RFO","JF","KER"), "Oil",
    code %in% c("NG","OG","PG","BFG","SGC","SGP","LFG"), "Gas",
    code %in% c("BLQ","WDS","WDL","WO","OBS","OBL","OBG","AB","MSW","TDF","SLW","PUR"), "Biomass",
    code %in% c("WAT","WND","SUN","GEO","NUC","WH"), "Renewable",
    default = "Unknown"
  )
}
gen_all[, convert_fuel := map_fuel(`Energy Source Code`)]

# ===============================================================
# 9. Function to select generators operating in a given year
# ===============================================================
select_generators_for_year <- function(dt, yr){
  dt[`Operating Year` <= yr & (is.na(`Retirement Year`) | `Retirement Year` > yr)]
}

# ===============================================================
# 10. Prepare for decade loop
# ===============================================================
target_years <- c(1940, 1950, 1960, 1970, 1980, 1990)
OUT_DIR <- "/Users/xshan2/Library/CloudStorage/OneDrive-GeorgeMasonUniversity-O365Production/GMU_PhD/01_Research/03_powerplant/data/pp_emissions_1940_1990"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

mm_formula <- ~ convert_fuel + `Unit.Type` + Year + scrubber + log_Gen_kWh - 1
train_cols <- colnames(model.matrix(mm_formula, data = train_data))

align_mm <- function(M, cols) {
  miss <- setdiff(cols, colnames(M))
  if (length(miss)) {
    M <- cbind(M, matrix(0, nrow(M), length(miss), dimnames = list(NULL, miss)))
  }
  extra <- setdiff(colnames(M), cols)
  if (length(extra)) M <- M[, setdiff(colnames(M), extra), drop = FALSE]
  M[, cols, drop = FALSE]
}

results_facility <- vector("list", length(target_years))
names(results_facility) <- as.character(target_years)

# ===============================================================
# 11. Main loop: estimate SO2 for each decade (1940–1990)
# ===============================================================
for (yr in target_years) {
  gen_y <- select_generators_for_year(gen_all, yr)
  if (nrow(gen_y) == 0L) {
    warning("No generators in operation for year ", yr)
    next
  }
  
  data_year <- gen_y[
    , .(
      Plant_ID = `Plant ID`, Generator_ID = `Generator ID`,
      Plant_Name = `Plant Name`, Plant_State = `Plant State`, County = County,
      Longitude = Longitude, Latitude = Latitude,
      convert_fuel = convert_fuel, max_capacity = `Nameplate Capacity (MW)`
    )
  ]
  
  data.table::setDT(data_year)
  data_year <- data_year[!is.na(max_capacity) & max_capacity > 0 & !is.na(convert_fuel)]
  
  data_year[, Year := as.integer(yr)]
  data_year[, convert_fuel := factor(convert_fuel, levels = levels(train_data$convert_fuel))]
  data_year[, `Unit.Type` := factor("Combustion turbine", levels = levels(train_data$Unit.Type))]
  data_year[, scrubber := 0L]
  data_year[, annual_gen_kWh := max_capacity * 1000 * 0.95 * 8760]
  data_year[, log_Gen_kWh := log10(annual_gen_kWh + 1)]
  
  data_year <- data_year[!is.na(convert_fuel)]
  
  X_tmp <- model.matrix(mm_formula, data = data_year)
  X_new <- align_mm(X_tmp, train_cols)
  
  pred_log_EF <- predict(xgb_model, newdata = xgb.DMatrix(X_new))
  data_year[, Predicted_EF_g_per_kWh := exp(pred_log_EF)]
  data_year[, SO2_grams := Predicted_EF_g_per_kWh * annual_gen_kWh]
  data_year[, SO2_tons := SO2_grams / 907185]
  
  keep_cols <- c("Year","Plant_ID","Generator_ID","Plant_Name","Plant_State","County",
                 "Longitude","Latitude","convert_fuel","max_capacity",
                 "annual_gen_kWh","Predicted_EF_g_per_kWh","SO2_tons")
  results_facility[[as.character(yr)]] <- data_year[, ..keep_cols]
  
  message("Done year ", yr, " (rows: ", nrow(data_year), ")")
}

# ===============================================================
# 12. Combine all decades and save results
# ===============================================================
combined_facility <- data.table::rbindlist(results_facility, use.names = TRUE, fill = TRUE)
out_csv <- file.path(OUT_DIR, "so2_facility_emissions_1940_1990.csv")
data.table::fwrite(combined_facility, out_csv)
message("Wrote combined facility CSV: ", out_csv, " (rows: ", nrow(combined_facility), ")")

# ===============================================================
# 13. Summaries: by Year × Fuel, and total by Year
# ===============================================================
combined_facility[, Year := as.integer(Year)]
combined_facility[, SO2_tons := as.numeric(SO2_tons)]

# 13.1 Year × Fuel summary
by_year_fuel <- combined_facility[
  , .(SO2_tons = sum(SO2_tons, na.rm = TRUE)),
  by = .(Year, convert_fuel)
][order(Year, convert_fuel)]
print(by_year_fuel)

# 13.2 Wide format
by_year_fuel_wide <- dcast(by_year_fuel, Year ~ convert_fuel, value.var = "SO2_tons", fill = 0)
print(by_year_fuel_wide)

# 13.3 Total per Year
by_year_total <- combined_facility[
  , .(SO2_tons_total = sum(SO2_tons, na.rm = TRUE)),
  by = Year
][order(Year)]
print(by_year_total)
