library(INLA)
library(ggplot2)
library(corrplot)
library(reshape2)
library(car)
library(readxl)
library(sf)
library(spdep)
library(writexl)
library(dplyr)

# ==== Set Working Directory & Load Data ====
setwd("~/INLA project")

data <- read_excel("expanded_bf_data.xlsx")

head(data)
summary(data)

# ================================================
# Data Preparation
# ================================================
# Extract country from the last 12 characters of district_country for vaccine introduction variable
n_last <- 12


# Extract country (last 12 characters)
data$country <- substr(data$district_country, 
                       nchar(data$district_country) - n_last + 1, 
                       nchar(data$district_country))

# Define vaccine introduction condition
data$vaccine <- ifelse(
  (data$district_country == "Sanmatenga Burkina Faso" & data$year_month >=201009) |
    (data$country == "Burkina Faso" & data$year_month >= 201052),
  1, 0
)



data$area <- as.numeric(as.factor(data$district_country))
data$ID.area <- as.numeric(data$area)

# cont month variable for time component
data<-data %>%
  arrange(year,month)%>%
  mutate(continuous_month = (year-min(year))*12+month)

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

###SCALE VARIABLES FOR MODEL

data$humidity_square <- (scale(data$humidity)^2)
data$rainfall_square <- (scale(data$rainfall)^2)
data$windspeed_square <- (scale(data$windspeed)^2)
data$north_wind_square <- (scale(data$north_wind)^2)
data$eastward_wind_square <- (scale(data$eastward_wind)^2)
data$aod_square <- (scale(data$aod)^2)
data$temp_square <- (scale(data$temp)^2)
data$pop_density_square <- (scale(data$pop_density)^2)
data$secmb_square <- (scale(data$secmb)^2)
data$total_cropland_square <- (scale(data$total_cropland)^2)



# ================================================
# Multivariate modelling
# ================================================

# first full model with all sig variables

# -----------------------------------------------
# Function to fit models with consistent settings
# -----------------------------------------------
fit_inla <- function(formula, data) {
  inla(
    formula,
    family = "binomial",
    data = data,
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
}

formula_full <- outbreak_occur ~ aod_square + eastward_wind_lag1 + temp_std + vaccine +
  humidity_lag2 + 
  f(continuous_month, model = 'ar1') +
  f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE)

model_full <- fit_inla(formula_full, data)
summary(model_full)

# ---- Extract Fixed Effects & Compute Odds Ratios ----
fixed <- model_full$summary.fixed
fixed <- fixed %>%
  mutate(
    OR = exp(mean),
    OR_low = exp(`0.025quant`),
    OR_high = exp(`0.975quant`),
    var = rownames(fixed)
  )

# ---- Plot Odds Ratios ----
ggplot(fixed, aes(x = var, y = OR, ymin = OR_low, ymax = OR_high)) +
  geom_pointrange(color = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  labs(
    x = "Predictor",
    y = "Odds ratio (95% CI)",
    title = "Odds ratios for fixed effects (INLA model)"
  ) +
  theme_minimal()




# ================================================
# 2. Backwards regression, these models  are all taking out different variables in turn
# ================================================
# Define formulas for model comparison
model_formulas <- list(
  model1 = outbreak_occur ~ aod_square + eastward_wind_lag1 + temp_std + vaccine +
    humidity_lag2 +  f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE),
  
  model2 = outbreak_occur ~ eastward_wind_lag1 + temp_std + vaccine +
    humidity_lag2 +
    f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE),
  
  model3 = outbreak_occur ~ aod_square + temp_std + vaccine +
    humidity_lag2 + 
    f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE),
  
  model4 = outbreak_occur ~ aod_square + eastward_wind_lag1 + vaccine +
    humidity_lag2 +
    f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE),
  
  model5 = outbreak_occur ~ aod_square + eastward_wind_lag1 + temp_std +
    humidity_lag2 + 
    f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE),
  
  
  model6 = outbreak_occur ~ aod_square + eastward_wind_lag1 + temp_std + vaccine +
    f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE)
)

# Fit all models
model_list <- lapply(model_formulas, fit_inla, data = data)

# ================================================
# Model Comparison Table (DIC & WAIC)
# ================================================
waic_dic_table <- data.frame(
  Model = names(model_list),
  DIC = sapply(model_list, function(x) x$dic$dic),
  WAIC = sapply(model_list, function(x) x$waic$waic)
)
print(waic_dic_table)

# Because models here ran relatively similar in terms of WAIC and DIC we decided to take out vaccination
# Vaccination (1/0) would potentially dominate over the other environmental effects and quash any outbreak risk post 2010


# ================================================
# 3.  Backwards regression part 2, these models  are all taking out different variables in turn
# ================================================

model_formulas <- list(
  model1 = outbreak_occur ~ temp_std + eastward_wind_lag1 + humidity_lag2 +
   f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE),
  
  model2 = outbreak_occur ~ aod_square + temp_std + eastward_wind_lag1 +
    f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE),
  
  model3 = outbreak_occur ~ eastward_wind_lag1 + aod_square +  humidity_lag2 + 
    f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE),
  
  model4 = outbreak_occur ~ aod_square +  humidity_lag2 +  temp_std +
    f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE)
)


# Fit all models
model_list <- lapply(model_formulas, fit_inla, data = data)

# ================================================
# Model Comparison Table (DIC & WAIC)
# ================================================
waic_dic_table <- data.frame(
  Model = names(model_list),
  DIC = sapply(model_list, function(x) x$dic$dic),
  WAIC = sapply(model_list, function(x) x$waic$waic)
)
print(waic_dic_table)
# Could potentially remove dust based on vars



# ================================================
# 4. Looking at impact of logit link vs clog log for rare events
# ================================================

formula <-
  outbreak_occur ~    humidity_lag2 + eastward_wind_lag1 + temp_std + aod_square +
  f(continuous_month, model = 'ar1') +
  f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE)

# default (logit)
res_logit <- inla(formula, family = "binomial", data = data, 
                  control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE)
)

# alternative: complementary log-log (better for rare positive events)
res_cloglog <- inla(formula, family = "binomial", data = data,
                    control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE),
                    control.family = list(link = "cloglog"))


p_hat <- res_logit$summary.fitted.values$mean

# AUPRC (PRROC expects scores for positives and negatives separately)
pr_res <- pr.curve(scores.class0 = p_hat[data$outbreak_occur == 1],
                   scores.class1 = p_hat[data$outbreak_occurk == 0],
                   curve = TRUE)
pr_res$auc.integral  # AUPRC

# AUROC
roc_obj <- roc(data$outbreak_occur, p_hat)
auc(roc_obj)
plot(roc_obj)

df_cal <- data %>%
  mutate(pred = p_hat) %>%
  mutate(bin = cut(pred, breaks = seq(0, 1, by = 0.05), include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(
    mean_pred = mean(pred, na.rm = TRUE),
    mean_obs  = mean(outbreak_occur, na.rm = TRUE),
    n = n()
  ) %>% ungroup()

ggplot(df_cal, aes(x = mean_pred, y = mean_obs, size = n)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Predicted probability", y = "Observed outbreak proportion",
       title = "Binned calibration plot") +
  xlim(0,1) + ylim(0,1) + theme_minimal()




p_hat2 <- res_cloglog$summary.fitted.values$mean



# AUPRC (PRROC expects scores for positives and negatives separately)
pr_res <- pr.curve(scores.class0 = p_hat[data$outbreak_occur == 1],
                   scores.class1 = p_hat[data$outbreak_occurk == 0],
                   curve = TRUE)
pr_res$auc.integral  # AUPRC

# AUROC
roc_obj <- roc(data$outbreak_occur, p_hat)
auc(roc_obj)
plot(roc_obj)

df_cal <- data %>%
  mutate(pred = p_hat) %>%
  mutate(bin = cut(pred, breaks = seq(0, 1, by = 0.05), include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(
    mean_pred = mean(pred, na.rm = TRUE),
    mean_obs  = mean(outbreak_occur, na.rm = TRUE),
    n = n()
  ) %>% ungroup()

ggplot(df_cal, aes(x = mean_pred, y = mean_obs, size = n)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Predicted probability", y = "Observed outbreak proportion",
       title = "Binned calibration plot") +
  xlim(0,1) + ylim(0,1) + theme_minimal()

# ================================================
# 5. Non linear effects of variables
# ================================================


# ================================================
# 5a AR1
# ================================================
model_formulas <- list(
  model1 = outbreak_occur ~ aod_square + eastward_wind_lag1 + temp_std +
    f(humidity, model = 'ar1') +
    f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE),
  
  model2 = outbreak_occur ~ humidity_lag2 + eastward_wind_lag1 + temp_std +
    f(aod, model = 'ar1') +
    f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE),
  
  model3 = outbreak_occur ~ humidity_lag2 + eastward_wind_lag1 + aod_square +
    f(temp, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE))


# Fit all models
model_list <- lapply(model_formulas, fit_inla, data = data)

# ================================================
# Model Comparison Table (DIC & WAIC)
# ================================================
waic_dic_table <- data.frame(
  Model = names(model_list),
  DIC = sapply(model_list, function(x) x$dic$dic),
  WAIC = sapply(model_list, function(x) x$waic$waic)
)
print(waic_dic_table)

# ================================================
# 5b rw1
# ================================================
 

data$hum_grp <- inla.group(data$humidity, n = 20)
data$aod_grp <- inla.group(data$humidity, n = 20)
model_formulas <- list(
  model1 = outbreak_occur ~ aod_square + eastward_wind_lag1 + temp_std +
    f(hum_grp, model = 'rw1') +
    f(continuous_month, model = 'rw1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE),
  
  model2 = outbreak_occur ~ humidity_lag2 + eastward_wind_lag1 + temp_std +
    f(aod_grp, model = 'rw1') +
    f(continuous_month, model = 'rw1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE),

  
  model3 = outbreak_occur ~ humidity_lag2 + aod_square + 
    eastward_wind_lag1 + 
    f(temp, model = "rw1") +   # non-linear temp
    f(area, model = 'besagproper', graph = ken.adj, constr = TRUE) +
    f(continuous_month, model = 'ar1'))


# Fit all models
model_list <- lapply(model_formulas, fit_inla, data = data)

# ================================================
# Model Comparison Table (DIC & WAIC)
# ================================================
waic_dic_table <- data.frame(
  Model = names(model_list),
  DIC = sapply(model_list, function(x) x$dic$dic),
  WAIC = sapply(model_list, function(x) x$waic$waic)
)

print(waic_dic_table)
# ================================================
# 5c rw2
# ================================================


model_formulas <- list(
  model1 = outbreak_occur ~ aod_square + eastward_wind_lag1 + temp_std +
    f(humidity, model = 'rw2') +
    f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE),
  
  model2 = outbreak_occur ~ humidity_lag2 + eastward_wind_lag1 + temp_std +
    f(aod, model = 'rw2') +
    f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE),
  
  model3 = outbreak_occur ~ humidity_lag2 + eastward_wind_lag1 + aod_square +
    f(temp, model = 'rw2') +
    f(continuous_month, model = 'ar1') +
    f(area, model = 'besagproper', graph = ken.adj, hyper = 'shyper', constr = TRUE))


# Fit all models
model_list <- lapply(model_formulas, fit_inla, data = data)

# ================================================
# Model Comparison Table (DIC & WAIC)
# ================================================
waic_dic_table <- data.frame(
  Model = names(model_list),
  DIC = sapply(model_list, function(x) x$dic$dic),
  WAIC = sapply(model_list, function(x) x$waic$waic)
)
print(waic_dic_table)



# ================================================
# FINAL WORKING MODEL
# ================================================


formula_nl <- outbreak_occur ~ humidity_lag2 + aod_square + 
  eastward_wind_lag1 + 
  f(temp, model = "rw1") +   # non-linear temp
  f(area, model = 'besagproper', graph = ken.adj, constr = TRUE) +
  f(continuous_month, model = 'ar1',
    hyper = list(theta = list(prior = "pc.cor1", param = c(0.7, 0.7))))

final_mod  <- inla(formula_nl,
                   family = "binomial",
                   data = data,
                   control.predictor = list(compute = TRUE),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))



pred_df <- data.frame(
  data,  # original data for district/month references
  fitted = final_mod$summary.fitted.values$mean,
  lower = final_mod$summary.fitted.values$`0.025quant`,
  upper = final_mod$summary.fitted.values$`0.975quant`
)





pred_df<-pred_df[ , c("district_country","month","year", "outbreak_occur" ,"humidity" ,"aod" , "eastward_wind", "temp","fitted", "lower", "upper"  )]


write_xlsx(pred_df,"district_month_predictions.xlsx",
)

# ================================================
# PLOTS FOR MODEL
# ================================================

# odds ratio
fixed <- final_mod$summary.fixed
fixed$OR <- exp(fixed$mean)
fixed$OR_low <- exp(fixed$`0.025quant`)
fixed$OR_high <- exp(fixed$`0.975quant`)
fixed$var <- rownames(fixed)

ggplot(fixed, aes(x = var, y = OR, ymin = OR_low, ymax = OR_high)) +
  geom_pointrange(color = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  labs(x = "Predictor", y = "Odds ratio",
       title = "Fixed effects on meningitis outbreak risk") +
  theme_minimal()


# roc curve
p_hat <- final_mod$summary.fitted.values$mean
obs   <- data$outbreak_occur  # observed outcomes (0/1)
--
library(pROC)

roc_obj <- roc(obs, p_hat)

# AUC value
auc_val <- auc(roc_obj)
cat("AUC:", round(auc_val, 3), "\n")
plot(roc_obj, col = "blue", lwd = 2, main = "ROC curve for meningitis outbreak model")
abline(a = 0, b = 1, lty = 2, col = "red")  # diagonal line = random classifier

# binned calibration plots
df_cal <- data %>%
  mutate(pred = p_hat) %>%
  mutate(bin = cut(pred, breaks = seq(0, 1, by = 0.05), include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(
    mean_pred = mean(pred, na.rm = TRUE),
    mean_obs  = mean(outbreak_occur, na.rm = TRUE),
    n = n()
  ) %>% ungroup()

ggplot(df_cal, aes(x = mean_pred, y = mean_obs, size = n)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Predicted probability", y = "Observed outbreak proportion",
       title = "Binned calibration plot") +
  xlim(0,1) + ylim(0,1) + theme_minimal()




# non linear effects of temperature plot
temp_eff <- final_mod$summary.random$temp
x_vals <- sort(unique(data$temp))

temp_eff$val <- x_vals[temp_eff$ID]


ggplot(temp_eff, aes(x = ID, y = mean)) +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`),
              fill = "skyblue", alpha = 0.4) +
  geom_line(color = "blue", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Standardised temperature", y = "Log-odds of outbreak",
       title = "RW1 effect of temperature") +
  theme_minimal()

