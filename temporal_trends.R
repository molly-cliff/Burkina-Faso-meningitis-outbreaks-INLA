
library(INLA)
library(ggplot2)
library(mgcv)
library(dplyr)
library(DAAG)
library(reshape2)

# Load and prepare data
setwd("~/INLA project")


spatiotemporaloutbreaks <- read.csv("full-outbreak-matched.csv")
newdata <- subset(spatiotemporaloutbreaks, country == "Burkina Faso")
myvars <- c("week", "month", "year", "outbreak", "district_country")
newdata <- newdata[myvars]
newdata$week_ar<-newdata$week

# Continuous week and month for later
newdata<-newdata %>%
  arrange(year,week)%>%
  mutate(continuous_week = (year-min(year))*52+week)


newdata<-newdata %>%
  arrange(year,month)%>%
  mutate(continuous_month = (year-min(year))*12+month)
head(newdata)

# Convert variables to numeric
newdata$week <- as.numeric(newdata$week)
newdata$outbreak <- as.numeric(newdata$outbreak)
newdata$month <- as.numeric(newdata$month)
newdata$year <- as.numeric(newdata$year)
newdata$continuous_week <- as.numeric(newdata$continuous_week)
newdata$continuous_month <- as.numeric(newdata$continuous_month)


# Aggregate by month
monthly_data <- newdata %>%
  group_by(year, month, district_country) %>%
  summarise(outbreak_occur = as.integer(any(outbreak > 0, na.rm =TRUE)),.groups = 'drop')

cat("Monthly aggregated data:\n")
print(head(monthly_data))

# Create a quarter variable
newdata <- newdata %>%
  mutate(quarter = ceiling(month / 3))

# Aggregate by quarter
quarterly_data <- newdata %>%
  group_by(year, quarter, district_country) %>%
  summarise(outbreak_occur = as.integer(any(outbreak > 0, na.rm =TRUE)),.groups = 'drop')

cat("Quarterly aggregated data:\n")
print(head(quarterly_data))


quarterly_data <-quarterly_data  %>%
  arrange(year,quarter)%>%
  mutate(continuous_quarter = (year-min(year))*4+quarter)
head(quarterly_data)

######################################
## TEMPORAL EFFECTS FOR WEEK ######
######################################

# This script fits three different INLA models to analyze temporal effects
# on weekly outbreak data using different random walk models (rw1, rw2) and an AR(1) model.
# The outcome variable is 'outbreak' and the temporal covariate is 'week'.

#########################################
# MODEL 1: WEEK EFFECT WITH RW2 (CYCLIC)
#########################################

# Define the model formula using a second-order random walk (RW2) with cyclic constraint
formula <- outbreak ~ f(week, model = "rw2", cyclic = TRUE)

# Fit the INLA model with a binomial likelihood
zinb_model_week_rw2 <- inla(
  formula,
  data = newdata,
  family = "binomial",
  control.predictor = list(compute = TRUE),  # Request computation of fitted values
  control.compute = list(dic = TRUE, waic = TRUE)  # Request model fit statistics
)

# Print model performance metrics and summaries
cat("Summary of the INLA Model for Week (rw2):\n")
cat("\nDIC:", zinb_model_week_rw2$dic$dic, "\n")   # Deviance Information Criterion
cat("WAIC:", zinb_model_week_rw2$waic$waic, "\n") # Watanabe-Akaike Information Criterion

# Display fixed and random effects
cat("\nFixed Effects Summary:\n")
print(zinb_model_week_rw2$summary.fixed)

cat("\nRandom Effects (Week) Summary:\n")
print(zinb_model_week_rw2$summary.random$week)

cat("\nModel Summary:\n")
summary(zinb_model_week_rw2)


#########################################
# MODEL 2: WEEK EFFECT WITH RW1 (CYCLIC)
#########################################

# Define the model using a first-order random walk (RW1) with cyclic constraint
formula <- outbreak ~ f(week, model = "rw1", cyclic = TRUE)

# Fit the INLA model
zinb_model_week_rw1 <- inla(
  formula,
  data = newdata,
  family = "binomial",  # Assuming zero-inflated binomial, but only 'binomial' is used
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

# Print model performance metrics and summaries
cat("Summary of the INLA Model for Week (rw1):\n")
cat("\nDIC:", zinb_model_week_rw1$dic$dic, "\n")
cat("WAIC:", zinb_model_week_rw1$waic$waic, "\n")

cat("\nFixed Effects Summary:\n")
print(zinb_model_week_rw1$summary.fixed)

cat("\nRandom Effects (Week) Summary:\n")
print(zinb_model_week_rw1$summary.random$week)

cat("\nModel Summary:\n")
summary(zinb_model_week_rw1)


#########################################
# MODEL 3: WEEK EFFECT WITH AR(1) (CYCLIC)
#########################################

# Define the model using an autoregressive model of order 1 (AR1) with cyclic constraint
formula <- outbreak ~ f(week, model = "ar1", cyclic = TRUE)

# Fit the INLA model
zinb_model_week_ar1 <- inla(
  formula,
  data = newdata,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

# Print model performance metrics and summaries
cat("Summary of the INLA Model for Week (ar1):\n")
cat("\nDIC:", zinb_model_week_ar1$dic$dic, "\n")
cat("WAIC:", zinb_model_week_ar1$waic$waic, "\n")

cat("\nFixed Effects Summary:\n")
print(zinb_model_week_ar1$summary.fixed)

cat("\nRandom Effects (Week) Summary:\n")
print(zinb_model_week_ar1$summary.random$week)

cat("\nModel Summary:\n")
summary(zinb_model_week_ar1)




######################################
## CONTINUOUS WEEK ANALYSIS ######
######################################

# This section fits an INLA model using 'continuous_week' as a continuous temporal variable
# to capture smooth trends in outbreak probability over time.

#########################################
# MODEL: WEEK EFFECT WITH RW1 (CONTINUOUS)
#########################################

# Define the model formula using a first-order random walk (RW1) 
# applied to a continuous version of 'week'. RW1 allows for flexible modeling
# of temporal dependence where changes are modeled as small deviations from previous values.
formula <- outbreak ~ f(continuous_week, model = "rw1")

# Fit the INLA model
zinb_model_week_rw1 <- inla(
  formula,
  data = newdata,
  family = "binomial",  # Binomial likelihood (could be zero-inflated if noted)
  control.predictor = list(compute = TRUE),  # Request fitted values computation
  control.compute = list(dic = TRUE, waic = TRUE)  # Request model fit criteria
)

# Display model performance metrics
cat("Summary of the INLA Model for Week (rw1):\n")
cat("\nDIC:", zinb_model_week_rw1$dic$dic, "\n")   # Deviance Information Criterion
cat("WAIC:", zinb_model_week_rw1$waic$waic, "\n") # W_




######################################
## TEMPORAL EFFECTS FOR WEEK ######
######################################

# This section explores temporal effects of 'week' while adjusting for
# fixed effects such as 'month' and 'year' using the INLA framework.
# It compares different random walk models (rw1, rw2) and includes fixed covariates.

#########################################
# MODEL 1: WEEK (RW1, CYCLIC) + MONTH
#########################################

# Define the model formula with:
# - 'month' as a fixed effect to control for monthly seasonality
# - 'week' modeled as a first-order cyclic random walk (RW1)
formula <- outbreak ~ month + f(week, model = "rw1", cyclic = TRUE)

# Fit the model using INLA
zinb_model_week_rw1 <- inla(
  formula,
  data = newdata,
  family = "binomial",  # Binomial likelihood for binary outcome
  control.predictor = list(compute = TRUE),  # Compute fitted values
  control.compute = list(dic = TRUE, waic = TRUE)  # Compute model fit criteria
)

# Output model fit statistics and summaries
cat("Summary of the INLA Model for Week (rw1) + Month:\n")
cat("\nDIC:", zinb_model_week_rw1$dic$dic, "\n")
cat("WAIC:", zinb_model_week_rw1$waic$waic, "\n")

cat("\nFixed Effects Summary:\n")
print(zinb_model_week_rw1$summary.fixed)

cat("\nRandom Effects (Week) Summary:\n")
print(zinb_model_week_rw1$summary.random$week)

cat("\nModel Summary:\n")
summary(zinb_model_week_rw1)

#########################################
# MODEL 2: WEEK (RW1, CYCLIC) + MONTH + YEAR
#########################################

# Define a similar model, now also adjusting for 'year'
formula <- outbreak ~ year + month + f(week, model = "rw1", cyclic = TRUE)

# Fit the updated model
zinb_model_week_rw1 <- inla(
  formula,
  data = newdata,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

# Output results
cat("Summary of the INLA Model for Week (rw1) + Month + Year:\n")
cat("\nDIC:", zinb_model_week_rw1$dic$dic, "\n")
cat("WAIC:", zinb_model_week_rw1$waic$waic, "\n")

cat("\nFixed Effects Summary:\n")
print(zinb_model_week_rw1$summary.fixed)

cat("\nRandom Effects (Week) Summary:\n")
print(zinb_model_week_rw1$summary.random$week)

cat("\nModel Summary:\n")
summary(zinb_model_week_rw1)

#########################################
# MODEL 3: WEEK (RW2, CYCLIC) ONLY + PRIOR
#########################################

# Define a model with only 'week' as a temporal effect using RW2
# Includes a PC prior for the precision (controls smoothness)
formula_rw2 <- outbreak ~  
  f(week, model = "rw2", cyclic = TRUE, 
    hyper = list(prec = list(prior = "pc.prec", param = c(0.01, 0.5))))
# PC prior: (0.01, 0.5) implies 1% probability that the standard deviation > sqrt(1/0.01) ≈ 10

# Fit the model
zinb_model_week_rw2 <- inla(
  formula,
  data = newdata,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

# Output results
cat("Summary of the INLA Model for Week (rw2) Only:\n")
cat("\nDIC:", zinb_model_week_rw2$dic$dic, "\n")
cat("WAIC:", zinb_model_week_rw2$waic$waic, "\n")

cat("\nFixed Effects Summary:\n")
print(zinb_model_week_rw2$summary.fixed)

cat("\nRandom Effects (Week) Summary:\n")
print(zinb_model_week_rw2$summary.random$week)

cat("\nModel Summary:\n")
summary(zinb_model_week_rw2)


#############################################
## TEMPORAL EFFECTS FOR MONTH – MODELS 4–5 ##
#############################################

### MODEL 4A: MONTHLY OUTBREAKS ~ RW2 on MONTH (Cyclic) ###
# Uses RW2 to capture smoother non-linear seasonal effects over months

formula <- outbreak_occur ~ f(month, model = "rw2", cyclic = TRUE)

zinb_model_week_rw2 <- inla(
  formula,
  data = monthly_data,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

cat("Model 4A - RW2 on Month:\n")
cat("\nDIC:", zinb_model_week_rw2$dic$dic, "\n")
cat("WAIC:", zinb_model_week_rw2$waic$waic, "\n")
print(zinb_model_week_rw2$summary.fixed)
print(zinb_model_week_rw2$summary.random$week)
summary(zinb_model_week_rw2)


### MODEL 4B: MONTHLY OUTBREAKS ~ RW1 on MONTH (Cyclic) ###

formula <- outbreak_occur ~ f(month, model = "rw1", cyclic = TRUE)

zinb_model_week_rw1 <- inla(
  formula,
  data = monthly_data,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

cat("Model 4B - RW1 on Month:\n")
cat("\nDIC:", zinb_model_week_rw1$dic$dic, "\n")
cat("WAIC:", zinb_model_week_rw1$waic$waic, "\n")
print(zinb_model_week_rw1$summary.fixed)
print(zinb_model_week_rw1$summary.random$week)
summary(zinb_model_week_rw1)


### MODEL 4C: MONTHLY OUTBREAKS ~ AR1 on MONTH (Cyclic) ###

formula <- outbreak_occur ~ f(month, model = "ar1", cyclic = TRUE)

zinb_model_week_ar1 <- inla(
  formula,
  data = monthly_data,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

cat("Model 4C - AR1 on Month:\n")
cat("\nDIC:", zinb_model_week_ar1$dic$dic, "\n")
cat("WAIC:", zinb_model_week_ar1$waic$waic, "\n")
print(zinb_model_week_ar1$summary.fixed)
print(zinb_model_week_ar1$summary.random$week)
summary(zinb_model_week_ar1)

#############################################
### MODEL 5: CONTINUOUS MONTH EFFECTS ###
#############################################

### MODEL 5A: OUTBREAK ~ RW2 on CONTINUOUS_MONTH ###

formula <- outbreak ~ f(continuous_month, model = "rw2")

zinb_model_week_rw2 <- inla(
  formula,
  data = newdata,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

cat("Model 5A - RW2 on Continuous Month:\n")
cat("\nDIC:", zinb_model_week_rw2$dic$dic, "\n")
cat("WAIC:", zinb_model_week_rw2$waic$waic, "\n")
print(zinb_model_week_rw2$summary.fixed)
print(zinb_model_week_rw2$summary.random$week)
summary(zinb_model_week_rw2)


### MODEL 5B: OUTBREAK ~ RW1 on CONTINUOUS_MONTH ###

formula <- outbreak ~ f(continuous_month, model = "rw1")

zinb_model_week_rw1 <- inla(
  formula,
  data = newdata,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

cat("Model 5B - RW1 on Continuous Month:\n")
cat("\nDIC:", zinb_model_week_rw1$dic$dic, "\n")
cat("WAIC:", zinb_model_week_rw1$waic$waic, "\n")
print(zinb_model_week_rw1$summary.fixed)
print(zinb_model_week_rw1$summary.random$week)
summary(zinb_model_week_rw1)


### MODEL 5C: OUTBREAK ~ AR1 on CONTINUOUS_MONTH ###

formula <- outbreak ~ f(continuous_month, model = "ar1")

zinb_model_week_ar1 <- inla(
  formula,
  data = newdata,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

cat("Model 5C - AR1 on Continuous Month:\n")
cat("\nDIC:", zinb_model_week_ar1$dic$dic, "\n")
cat("WAIC:", zinb_model_week_ar1$waic$waic, "\n")
print(zinb_model_week_ar1$summary.fixed)
print(zinb_model_week_ar1$summary.random$week)
summary(zinb_model_week_ar1)

#############################################
### MODEL COMPARISON + METRICS ###
#############################################

### MODEL 6A: MONTHLY OUTBREAKS ~ AR1 + ZIB ###
# Zero-inflated binomial model with AR1 effect on month

formula <- outbreak_occur ~ f(month, model = "ar1", cyclic = TRUE)

zinb_model_week_ar1 <- inla(
  formula,
  data = monthly_data,
  family = "zeroinflatedbinomial1",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

cat("Model 6A - AR1 + ZIB:\n")
cat("\nDIC:", zinb_model_week_ar1$dic$dic, "\n")
cat("WAIC:", zinb_model_week_ar1$waic$waic, "\n")
print(zinb_model_week_ar1$summary.fixed)
print(zinb_model_week_ar1$summary.random$week)
summary(zinb_model_week_ar1)


### MODEL 6B: ZIB + AR1 + YEAR ###

formula <- outbreak_occur ~ year + f(month, model = "ar1", cyclic = TRUE)

month_year_ar1 <- inla(
  formula,
  data = monthly_data,
  family = "zeroinflatedbinomial1",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

cat("Model 6B - AR1 + ZIB + Year:\n")
cat("\nDIC:", month_year_ar1$dic$dic, "\n")
cat("WAIC:", month_year_ar1$waic$waic, "\n")
print(month_year_ar1$summary.fixed)
print(month_year_ar1$summary.random$week)
summary(month_year_ar1)

# Store fitted values for later evaluation
monthly_data$month_year_ar1 <- month_year_ar1$summary.fitted.values[,1]
monthly_data$month_ar1 <- zinb_model_week_ar1$summary.fitted.values[,1]

#############################################
### EVALUATION: SENSITIVITY & SPECIFICITY ###
#############################################

sensitivity_specificity <- function(x, y) {
  predicted_classes <- ifelse(x > 0.1, 1, 0)
  conf_matrix <- table(Actual = y, Predicted = predicted_classes)
  
  print("Confusion Matrix:")
  print(conf_matrix)
  
  predicted_classes <- as.factor(predicted_classes)
  monthly_data$outbreak <- as.factor(y)
  
  # Calculate sensitivity and specificity
  sensitivity <- sensitivity(y, predicted_classes)
  specificity <- specificity(y, predicted_classes)
  
  cat("Sensitivity:", sensitivity)
  cat(" Specificity:", specificity)
}

# Set outcome as factor
monthly_data$outbreak <- as.factor(monthly_data$outbreak_occur)

# Evaluate both models
sensitivity_specificity(monthly_data$month_year_ar1, monthly_data$outbreak)
sensitivity_specificity(monthly_data$month_ar1, monthly_data$outbreak)

