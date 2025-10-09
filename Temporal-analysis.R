############################################################
# SPATIOTEMPORAL OUTBREAK ANALYSIS USING INLA
# Example: Burkina Faso outbreak data
# Author: [Your Name]
# Date: [Insert Date]
############################################################


##############################
# 1. LOAD LIBRARIES
##############################
library(INLA)
library(ggplot2)
library(mgcv)
library(dplyr)
library(DAAG)
library(reshape2)
library(readxl)
library(writexl)
library(lubridate)
library(scales)  # for plogis


##############################
# 2. LOAD AND PREPARE DATA
##############################
setwd("~/INLA project")

# Option 1: From CSV
spatiotemporaloutbreaks <- read.csv("full-outbreak-matched.csv")
newdata <- subset(spatiotemporaloutbreaks, country == "Burkina Faso")

# Option 2: From Excel (preferred)
newdata <- read_excel("BF_outbreak.xlsx")

# Keep only relevant variables
myvars <- c("week", "month", "year", "outbreak", "district_country")
newdata <- newdata[myvars]


##############################
# 3. DATA CLEANING & TYPE CONVERSION
##############################
newdata$week <- as.numeric(newdata$week)
newdata$outbreak <- as.numeric(newdata$outbreak)
newdata$month <- as.numeric(newdata$month)
newdata$year <- as.numeric(newdata$year)


##############################
# 4. AGGREGATE BY MONTH & CREATE CONTINUOUS MONTH INDEX
##############################
monthly_data <- newdata %>%
  group_by(year, month, district_country) %>%
  summarise(outbreak_occur = as.integer(any(outbreak > 0, na.rm = TRUE)),
            .groups = 'drop')

monthly_data <- monthly_data %>%
  arrange(year, month) %>%
  mutate(continuous_month = (year - min(year)) * 12 + month)

monthly_data$continuous_month <- as.numeric(monthly_data$continuous_month)

cat("Monthly aggregated data:\n")
print(head(monthly_data))


############################################################
# 5. TEMPORAL EFFECTS FOR MONTH – CYCLIC MODELS
############################################################

### MODEL 1A: RW2 on MONTH (Cyclic)
formula <- outbreak_occur ~ f(month, model = "rw2", cyclic = TRUE)
rw2_month_cyclic <- inla(
  formula,
  data = monthly_data,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)
dic_rw2_month_cyclic <- rw2_month_cyclic$dic$dic
waic_rw2_month_cyclic <- rw2_month_cyclic$waic$waic
cpo_valuesrw2_month_cyclic <- -mean(log(rw2_month_cyclic$cpo$cpo))

summary(rw2_month_cyclic)


### MODEL 1B: RW1 on MONTH (Cyclic)
formula <- outbreak_occur ~ f(month, model = "rw1", cyclic = TRUE)
rw1_month_cyclic <- inla(
  formula,
  data = monthly_data,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)
dic_rw1_month_cyclic <- rw1_month_cyclic$dic$dic
waic_rw1_month_cyclic <- rw1_month_cyclic$waic$waic
cpo_valuesrw1_month_cyclic <- -mean(log(rw1_month_cyclic$cpo$cpo))

summary(rw1_month_cyclic)


### MODEL 1C: AR1 on MONTH (Cyclic)
formula <- outbreak_occur ~ f(month, model = "ar1", cyclic = TRUE)
ar1_month_cyclic <- inla(
  formula,
  data = monthly_data,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)
dic_ar1_month_cyclic <- ar1_month_cyclic$dic$dic
waic_ar1_month_cyclic <- ar1_month_cyclic$waic$waic
cpo_valuesar1_month_cyclic <- -mean(log(ar1_month_cyclic$cpo$cpo))

summary(ar1_month_cyclic)


############################################################
# 6. TEMPORAL EFFECTS – CONTINUOUS MONTH MODELS
############################################################

### MODEL 2A: RW2 on CONTINUOUS_MONTH
formula <- outbreak_occur ~ f(continuous_month, model = "rw2")
rw2_month_cont <- inla(
  formula,
  data = monthly_data,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)
dic_rw2_month_cont <- rw2_month_cont$dic$dic
waic_rw2_month_cont <- rw2_month_cont$waic$waic
cpo_valuesrw2_month_cont <- -mean(log(rw2_month_cont$cpo$cpo))

summary(rw2_month_cont)


### MODEL 2B: RW1 on CONTINUOUS_MONTH
formula <- outbreak_occur ~ f(continuous_month, model = "rw1")
rw1_month_cont <- inla(
  formula,
  data = monthly_data,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)
dic_rw1_month_cont <- rw1_month_cont$dic$dic
waic_rw1_month_cont <- rw1_month_cont$waic$waic
cpo_valuesrw1_month_cont <- -mean(log(rw1_month_cont$cpo$cpo))

summary(rw1_month_cont)


### MODEL 2C: AR1 on CONTINUOUS_MONTH
formula <- outbreak_occur ~ f(continuous_month, model = "ar1")
ar1_month_cont <- inla(
  formula,
  data = monthly_data,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)
dic_ar1_month_cont <- ar1_month_cont$dic$dic
waic_ar1_month_cont <- ar1_month_cont$waic$waic
cpo_valuesar1_month_cont <- -mean(log(ar1_month_cont$cpo$cpo))

summary(ar1_month_cont)


############################################################
# 7. MODEL COMPARISON TABLE
############################################################
dic <- data.frame(
  criteria = c("DIC", "WAIC", "LogCPO"),
  rw2_month_cyclic = c(dic_rw2_month_cyclic, waic_rw2_month_cyclic, cpo_valuesrw2_month_cyclic),
  rw1_month_cyclic = c(dic_rw1_month_cyclic, waic_rw1_month_cyclic, cpo_valuesrw1_month_cyclic),
  ar1_month_cyclic = c(dic_ar1_month_cyclic, waic_ar1_month_cyclic, cpo_valuesar1_month_cyclic),
  rw2_month_cont = c(dic_rw2_month_cont, waic_rw2_month_cont, cpo_valuesrw2_month_cont),
  rw1_month_cont = c(dic_rw1_month_cont, waic_rw1_month_cont, cpo_valuesrw1_month_cont),
  ar1_month_cont = c(dic_ar1_month_cont, waic_ar1_month_cont, cpo_valuesar1_month_cont)
)
print(dic)
write_xlsx(dic, path = "model_comparison_temporal.xlsx")


############################################################
# 8. VISUALIZE TEMPORAL EFFECTS (AR1 Continuous Month)
############################################################

# Extract temporal random effects
temporal_effects <- ar1_month_cont$summary.random$continuous_month %>%
  dplyr::select(ID, mean, `0.025quant`, `0.975quant`) %>%
  rename(month_index = ID,
         effect_mean = mean,
         effect_lci = `0.025quant`,
         effect_uci = `0.975quant`)

# Convert month index to dates (assuming start = Jan 2003)
temporal_effects$date <- ymd("2003-01-01") + months(temporal_effects$month_index - 1)

# Plot temporal effects
ggplot(temporal_effects, aes(x = date, y = effect_mean)) +
  geom_ribbon(aes(ymin = effect_lci, ymax = effect_uci),
              fill = "lightblue", alpha = 0.4) +
  geom_line(color = "blue", size = 1) +
  theme_minimal() +
  labs(
    title = "Temporal Random Effects from INLA Model",
    subtitle = "Posterior mean and 95% credible intervals",
    x = "Time (Year)",
    y = "Temporal effect (log-odds of outbreak)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
