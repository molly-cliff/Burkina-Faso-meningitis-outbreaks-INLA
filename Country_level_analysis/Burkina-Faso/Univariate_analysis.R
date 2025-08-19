library(INLA)
library(ggplot2)
library(corrplot)
library(reshape2)
library(car)
library(readxl)
library(sf)
library(spdep)
library(writexl)

# ==== Set Working Directory & Load Data ====
setwd("~/INLA project/landcovertest3")

data <- read_excel("testing.xlsx")


# View the first few rows of the data
head(data)
summary(data)

data$year_month2 <- as.numeric(paste0(data$year, sprintf("%02d", data$month)))
data$year_month2 <-as.numeric(data$year_month2 )

# Extract country from the last 12 characters of district_country for vaccine introduction variable
n_last <- 12


# Extract country (last 12 characters)
data$country <- substr(data$district_country, 
                       nchar(data$district_country) - n_last + 1, 
                       nchar(data$district_country))

# Define vaccine introduction condition
data$vaccine <- ifelse(
  (data$district_country == "Sanmatenga Burkina Faso" & data$year_month2 >=201009) |
    (data$country == "Burkina Faso" & data$year_month2 >= 201052),
  1, 0
)


data$total_cropland<- data$c3ann + data$c3per + data$c4ann +data$c4per
# ==== Correlation Analysis ====
plot_correlation <- function(vars, title = "Correlation Plot") {
  corr_matrix <- cor(data[, vars], use = "complete.obs")
  
  ggplot(melt(corr_matrix), aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limit = c(-1, 1), name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_fixed() +
    labs(title = title)
  
  corrplot(corr_matrix, method = "color", type = "upper", 
           addCoef.col = "black", number.cex = 1.2, tl.cex = 1.2,
           tl.col = "black", col = colorRampPalette(c("#ff0000", "white", "#195696"))(10000))
}

plot_correlation(c("rainfall", "windspeed", "eastward_wind", "north_wind", "humidity", "aod", "secmb", "total_cropland" , "vaccine","temp"))
plot_correlation(c("windspeed", "eastward_wind", "humidity", "aod", "secmb", "total_cropland" , "vaccine"))

# ==== VIF Analysis ====

test_model<-glm(outbreak_occur ~ aod + eastward_wind +  humidity + windspeed + total_cropland + secmb + vaccine + temp, data = data)
 
vif(test_model)

data$area <- as.numeric(as.factor(data$district_country))
data$ID.area <- as.numeric(data$area)

####NEIGHBOURHOOD MATRICES FOR BYM MODEL ====

setwd("~/INLA project")
shape2 <- st_read("Burkinafaso.shp")
valid_shapefile<-st_is_valid(shape2)
table(valid_shapefile)
### Make spatial neighborhood structure
shapefile_spatial <- as_Spatial(shape2)
nb <- poly2nb(shapefile_spatial, queen = TRUE)

Coords <- coordinates(shapefile_spatial)
plot(shapefile_spatial, col = adjustcolor("lightblue", alpha.f = 0.5), border = "darkgray", lwd = 0.7)
plot(nb, coords = Coords, add = TRUE, lwd = 1, col = "blue")

nb2INLA("map.adj", nb)
ken.adj <- paste0(getwd(), "/map.adj")


### UNIVARIATE INLA MODELS W/ SPATIAL & TEMPORAL COMPONENTS ====



# Add scaled to variable list
vars_inla <- c("windspeed", "eastward_wind", "aod", "secmb", "total_cropland" , "vaccine", "temp")

# Initialize result container
posterior_table_m1 <- data.frame()

for (var in vars_inla) {
  formula <- as.formula(paste(
    "outbreak_occur ~", var,
    "+ f(month, model = 'ar1', cyclic = TRUE)",
    "+ f(area, model = 'bym2', graph = ken.adj, hyper = 'shyper', scale.model = TRUE, constr = TRUE)"
  ))
  
  M5 <- inla(
    formula,
    data = data,
    family = "binomial",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
  )
  
  effect_row <- M5$summary.fixed[grepl(var, rownames(M5$summary.fixed)), ]
  lcpocpo_value <- -mean(log(M5$cpo$cpo), na.rm = TRUE)
  
  posterior_table_m1 <- rbind(posterior_table_m1, data.frame(
    covariate = var,
    estimate = effect_row$mean,
    CI_1 = effect_row$`0.025quant`,
    CI_2 = effect_row$`0.975quant`,
    waic = M5$waic$waic,
    dic = M5$dic$dic,
    cpo = lcpocpo_value
  ))
  
  cat("\n--- INLA Model Summary for:", var, "---\n")
  cat("Mean:", effect_row$mean, "\n")
  cat("95% CI:", effect_row$`0.025quant`, "-", effect_row$`0.975quant`, "\n")
  cat("Significant?:", ifelse(effect_row$`0.025quant` > 0 | effect_row$`0.975quant` < 0, "Yes", "No"), "\n")
  cat("DIC:", M5$dic$dic, ", WAIC:", M5$waic$waic, "\n")
}

# Final cleanup and export
rownames(posterior_table_m1) <- 1:nrow(posterior_table_m1)
print(posterior_table_m1)
write_xlsx(posterior_table_m1, "inla_covariate_results_with_space_time.xlsx")

# Plot the results using ggplot2
ggplot(data = posterior_table_m1, aes(y = covariate, x = estimate)) + 
  geom_pointrange(aes(xmin = CI_1, xmax = CI_2), color = "blue", fill = "white", shape = 22) +
  theme_minimal() +
  labs(title = "INLA Model Estimates with 95% Credible Intervals", x = "Estimate", y = "Covariate")

# Compute the odds ratios (OR) by exponentiating the estimates and the CI bounds
posterior_table_m1_or <- posterior_table_m1
posterior_table_m1_or[, c("estimate", "CI_1", "CI_2")] <- exp(posterior_table_m1_or[, c("estimate", "CI_1", "CI_2")])



### UNIVARIATE INLA MODELS W/ SPATIAL & TEMPORAL COMPONENTS, VARS STANDARDISED ====

test<-table(data$outbreak_occur)
test<-as.data.frame(test)
outbreak<-test[2,2]
no_outbreak<-test[1,2]
outbreak<-as.numeric(outbreak)
no_outbreak<-as.numeric(no_outbreak)
total<-no_outbreak+outbreak

non_outbreak_weight <- outbreak/total
outbreak_weight<-1-non_outbreak_weight

weights <- ifelse(data$outbreak_occur == 1, outbreak_weight, non_outbreak_weight) 


data$humidity_std <- scale(data$humidity)
data$rainfall_std <- scale(data$rainfall)
data$windspeed_std <- scale(data$windspeed)
data$north_wind_std <- scale(data$north_wind)
data$eastward_wind_std <- scale(data$eastward_wind)
data$aod_std <- scale(data$aod)
data$temp_std <- scale(data$temp)


# Add scaled to variable list
vars_inla <- c("windspeed_std","eastward_wind_std" ,"humidity_std" ,"eastward_wind_std", "secmb", "aod_std", "total_cropland", "vaccine", "temp_std")

# Initialize result container
posterior_table_m1 <- data.frame()

for (var in vars_inla) {
  formula <- as.formula(paste(
    "outbreak_occur ~", var,
    "+ f(month, model = 'ar1', cyclic = TRUE)",
    "+ f(area, model = 'bym2', graph = ken.adj, hyper = 'shyper', scale.model = TRUE, constr = TRUE)"
  ))
  
  M5 <- inla(
    formula,
    data = data,
    family = "binomial",
    weights=weights,
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
  )
  
  effect_row <- M5$summary.fixed[grepl(var, rownames(M5$summary.fixed)), ]
  lcpocpo_value <- -mean(log(M5$cpo$cpo), na.rm = TRUE)
  
  posterior_table_m1 <- rbind(posterior_table_m1, data.frame(
    covariate = var,
    estimate = effect_row$mean,
    CI_1 = effect_row$`0.025quant`,
    CI_2 = effect_row$`0.975quant`,
    waic = M5$waic$waic,
    dic = M5$dic$dic,
    cpo = lcpocpo_value
  ))
  
  cat("\n--- INLA Model Summary for:", var, "---\n")
  cat("Mean:", effect_row$mean, "\n")
  cat("95% CI:", effect_row$`0.025quant`, "-", effect_row$`0.975quant`, "\n")
  cat("Significant?:", ifelse(effect_row$`0.025quant` > 0 | effect_row$`0.975quant` < 0, "Yes", "No"), "\n")
  cat("DIC:", M5$dic$dic, ", WAIC:", M5$waic$waic, "\n")
}

# Final cleanup and export
rownames(posterior_table_m1) <- 1:nrow(posterior_table_m1)
print(posterior_table_m1)
write_xlsx(posterior_table_m1, "inla_covariate_results_with_space_time.xlsx")

# Plot the results using ggplot2
ggplot(data = posterior_table_m1, aes(y = covariate, x = estimate)) + 
  geom_pointrange(aes(xmin = CI_1, xmax = CI_2), color = "blue", fill = "white", shape = 22) +
  theme_minimal() +
  labs(title = "INLA Model Estimates with 95% Credible Intervals", x = "Estimate", y = "Covariate")

# Compute the odds ratios (OR) by exponentiating the estimates and the CI bounds
posterior_table_m1_or <- posterior_table_m1
posterior_table_m1_or[, c("estimate", "CI_1", "CI_2")] <- exp(posterior_table_m1_or[, c("estimate", "CI_1", "CI_2")])
rownames(posterior_table_m1) <- 1:nrow(posterior_table_m1)
print(posterior_table_m1)
write_xlsx(posterior_table_m1, "inla_covariate_results_with_space_time_weighted.xlsx")

final_columns <- c("code", "year", "month", "district_country", "outbreak_occur", "year_month",
                   "rainfall", "windspeed", "eastward_wind", "north_wind", "humidity", "aod",
                   "secmb", "total_cropland", "temp")

final_data <- data[, final_columns]
write_xlsx(final_data, "full_BF_data.xlsx")

