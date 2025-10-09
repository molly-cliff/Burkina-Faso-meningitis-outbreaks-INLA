# Libraries
library(INLA)
library(ggplot2)
library(mgcv)
library(dplyr)
library(readxl)
library(writexl)
library(sf)
library(spdep)   
library(sp)    
library(viridis)  
library(patchwork)


setwd("~/INLA project")  # change as needed

# -------------------------
# 1. Load data
# -------------------------
# If you have both CSV and Excel, read the correct one. Here we use the Excel file:
monthly_raw <- read_excel("BF_outbreak.xlsx")

# Keep only relevant variables (update names to match your sheet)
keep_vars <- c("week", "month", "year", "outbreak", "district_country", "ADMN2")
missing_vars <- setdiff(keep_vars, names(monthly_raw))
if (length(missing_vars) > 0) {
  stop("Missing variables in data: ", paste(missing_vars, collapse = ", "))
}
monthly_raw <- monthly_raw[ , keep_vars]

# Convert numeric-ish columns
monthly_raw <- monthly_raw %>%
  mutate(
    week = as.numeric(week),
    month = as.numeric(month),
    year = as.numeric(year),
    outbreak = as.numeric(outbreak),
    district_country = as.character(district_country),
    ADMN2 = as.character(ADMN2)
  )

# -------------------------
# 2. Prepare monthly dataset
# -------------------------
# Collapse to one row per district-year-month indicating whether any outbreak occurred
monthly_data <- monthly_raw %>%
  group_by(year, month, district_country, ADMN2) %>%
  summarise(
    outbreak_occur = as.integer(any(outbreak > 0, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(year, month)

# Create continuous month index (starting at 1)
monthly_data <- monthly_data %>%
  mutate(
    continuous_month = (year - min(year, na.rm = TRUE)) * 12 + month
  )

# Time index for INLA AR1
monthly_data <- monthly_data %>%
  mutate(
    time_id = as.integer(factor(continuous_month, 
                                levels = sort(unique(continuous_month))))
  )

# -------------------------
# 3. Read shapefile and build adjacency
# -------------------------
shape <- st_read("Burkinafaso.shp")  # ensure it is the adm2 shapefile

# Check expected columns exist
if (!all(c("NAME_2", "COUNTRY") %in% names(shape))) {
  stop("Expected NAME_2 and COUNTRY fields in shapefile not found.")
}

# Create a district-country key in the shapefile to match monthly_data
shape$district_country <- paste(shape$NAME_2, shape$COUNTRY, sep = " ")

# Check for mismatches between data and shape
missing_in_shape <- setdiff(unique(monthly_data$district_country), shape$district_country)
if (length(missing_in_shape) > 0) {
  stop("These districts exist in the data but not in the shapefile: ",
       paste(missing_in_shape, collapse = ", "))
}

missing_in_data <- setdiff(shape$district_country, unique(monthly_data$district_country))
if (length(missing_in_data) > 0) {
  message("Dropping ", length(missing_in_data), " districts from shapefile that are not in data.")
  shape <- shape %>% filter(district_country %in% unique(monthly_data$district_country))
}

# Ensure both are sorted in the same order for consistent ID mapping
shape <- shape %>% arrange(district_country)
monthly_data <- monthly_data %>% arrange(district_country)

# Assign numeric area ID in monthly_data that matches the order in shape
monthly_data$area <- match(monthly_data$district_country, shape$district_country)
if (any(is.na(monthly_data$area))) stop("NA area IDs after matching; check district keys.")

# Build adjacency (nb) from shapefile polygons
# Convert to Spatial for poly2nb if necessary
shape_sp <- as(shape, "Spatial")
nb <- poly2nb(shape_sp, queen = TRUE)

# Write adjacency in INLA format and read as graph
nb2INLA("map.adj", nb)
graph_adj <- inla.read.graph("map.adj")  # graph object for INLA

# -------------------------
# 4. Hyperparameters for spatial models
# -------------------------
# Note: BYM2 typically uses a specific parameterization; here we set a PC prior on precision.
shyper <- list(prec = list(prior = "pc.prec", param = c(1, 0.01))) 
# adjust the param values to suit prior beliefs (scale, tail prob)

# -------------------------
# 5. Model formulas and fits
# -------------------------
# Use an index for the AR1 time effect (time_id)
# Use family = "binomial" with outbreak_occur (0/1)
# Random effect for area uses the graph 'graph_adj'

# Helper function to compute -mean(log(cpo)) safely
compute_logcpo <- function(model) {
  cpo <- model$cpo$cpo
  cpo <- cpo[!is.na(cpo) & cpo > 0]
  if (length(cpo) == 0) return(NA_real_)
  return(-mean(log(cpo)))
}

# Model 1: BYM2
formula_bym2 <- outbreak_occur ~
  f(time_id, model = "ar1") +
  f(area, model = "bym2", graph = graph_adj, hyper = shyper, scale.model = TRUE)

fit_bym2 <- inla(
  formula_bym2,
  data = monthly_data,
  family = "binomial",
  control.predictor = list(compute = TRUE, link =  "logit"),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# Model 2: BYM (classic)
formula_bym <- outbreak_occur ~
  f(time_id, model = "ar1") +
  f(area, model = "bym", graph = graph_adj, scale.model = TRUE)

fit_bym <- inla(
  formula_bym,
  data = monthly_data,
  family = "binomial",
  control.predictor = list(compute = TRUE, link = "logit"),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# Model 3: Besag (ICAR)
formula_besag <- outbreak_occur ~
  f(time_id, model = "ar1") +
  f(area, model = "besag", graph = graph_adj, scale.model = TRUE)

fit_besag <- inla(
  formula_besag,
  data = monthly_data,
  family = "binomial",
  control.predictor = list(compute = TRUE, link = "logit"),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# Model 4: Besag proper
formula_besagproper <- outbreak_occur ~
  f(time_id, model = "ar1") +
  f(area, model = "besagproper", graph = graph_adj)

fit_besagproper <- inla(
  formula_besagproper,
  data = monthly_data,
  family = "binomial",
  control.predictor = list(compute = TRUE, link = "logit"),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# -------------------------
# 6. Model diagnostics and comparison
# -------------------------
results <- data.frame(
  model = c("bym", "bym2", "besag", "besagproper"),
  DIC = c(fit_bym$dic$dic, fit_bym2$dic$dic, fit_besag$dic$dic, fit_besagproper$dic$dic),
  WAIC = c(fit_bym$waic$waic, fit_bym2$waic$waic, fit_besag$waic$waic, fit_besagproper$waic$waic),
  LogCPO = c(
    compute_logcpo(fit_bym),
    compute_logcpo(fit_bym2),
    compute_logcpo(fit_besag),
    compute_logcpo(fit_besagproper)
  ),
  stringsAsFactors = FALSE
)

print(results)
write_xlsx(results, path = "model_comparison_spatial.xlsx")

# Print summaries (fixed and random (area) summaries)
print("BYM2 fixed effects:")
print(fit_bym2$summary.fixed)
print("BYM2 area random summary (first rows):")
print(head(fit_bym2$summary.random$area))

print("BYM fixed effects:")
print(fit_bym$summary.fixed)
print("BESAG fixed effects:")
print(fit_besag$summary.fixed)
print("BESAGPROPER fixed effects:")
print(fit_besagproper$summary.fixed)

# -------------------------
# 7. Visualization: spatial effects per model
# -------------------------
plot_spatial_effects <- function(model, shape_sf, monthly_data, title) {
  # Extract the random effect table for 'area'
  if (!("area" %in% names(model$summary.random))) {
    stop("Model does not contain a random effect called 'area'")
  }
  df_area <- model$summary.random$area %>%
    mutate(ID = as.integer(ID)) %>%
    rename(mean_effect = mean)
  
  # Create lookup table from area ID -> district_country
  lookup <- monthly_data %>%
    distinct(area, district_country) %>%
    arrange(area)
  
  df_area <- df_area %>%
    left_join(lookup, by = c("ID" = "area"))
  
  # Merge summary into spatial object
  shape_plot <- shape_sf %>%
    left_join(df_area, by = "district_country")
  
  p <- ggplot(shape_plot) +
    geom_sf(aes(fill = mean_effect), color = "gray40", size = 0.3) +
    scale_fill_viridis_c(option = "plasma", name = "Posterior mean") +
    labs(title = title, subtitle = "Spatial random effect (area)") +
    theme_minimal()
  
  return(p)
}

p_bym  <- plot_spatial_effects(fit_bym, shape, monthly_data, "BYM Model — Spatial Effects")
p_bym2 <- plot_spatial_effects(fit_bym2, shape, monthly_data, "BYM2 Model — Spatial Effects")
p_besag <- plot_spatial_effects(fit_besag, shape, monthly_data, "Besag Model — Spatial Effects")
p_besagproper <- plot_spatial_effects(fit_besagproper, shape, monthly_data, "Besag Proper Model — Spatial Effects")

# Show plots in a 2x2 grid
(p_bym | p_bym2) / (p_besag | p_besagproper)

# -------------------------
# 8. Save model objects for later use
# -------------------------
save(fit_bym, fit_bym2, fit_besag, fit_besagproper, file = "inla_fits_spatial.RData")
