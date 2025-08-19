library(INLA)
library(tidyverse) 
library(readxl)
library(car)
library(sf)
library(spdep)
library(broom)
library(writexl)
library(mgcv)
library(pROC)
library(DAAG)
library(caret)
library(gridExtra)
library(haven)
library(tmap)
library(MLmetrics)
library(DescTools)
library(ggpubr)
library(reshape2)

sensitivity_specificity <- function(predicted_probs, actual) {
  predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
  conf_matrix <- table(Actual = actual, Predicted = predicted_classes)
  cat("Confusion Matrix:\n")
  print(conf_matrix)
  
  predicted_classes <- factor(predicted_classes, levels = c(0, 1))
  actual <- factor(actual, levels = c(0, 1))
  
  sensitivity_value <- caret::sensitivity(predicted_classes, actual, positive = "1")
  specificity_value <- caret::specificity(predicted_classes, actual, negative = "0")
  
  cat("Sensitivity:", sensitivity_value, "\n")
  cat("Specificity:", specificity_value, "\n")
  
  return(data.frame(
    Sensitivity = sensitivity_value,
    Specificity = specificity_value
  ))
}



setwd("~/INLA project")
data <- read_excel("full_BF_data.xlsx")
head(data)
summary(data)


# Create year-month composite variable
data$year_month2 <- as.numeric(paste0(data$year, sprintf("%02d", data$month)))
data$year_month2 <-as.numeric(data$year_month2 )
n_last <- 12 
data$country<-substr(data$district_country, nchar(data$district_country,) - n_last + 1, nchar(data$district_country,))

# Define vaccine introduction condition from Trotter et al
data$vaccine <- ifelse(
  (data$district_country == "Sanmatenga Burkina Faso" & data$year_month2 >=201009) |
    (data$country == "Burkina Faso" & data$year_month2 >= 201052),
  1, 0
)

summary(data)
colSums(is.na(data))

# Save cleaned data, this is the same enviromental data but with the newly created vaccine column
write_xlsx(data, "full_BF_data_vaccine.xlsx")

setwd("~/INLA project")

#Reading in outbreak data, no enviromental covariates in this data set
#This is to first show the best monthly space time model
monthly_data <- read_excel("BF_outbreak.xlsx")

table(monthly_data$outbreak)
setwd("~/INLA project")
shape2 <- st_read("Burkinafaso.shp")
windows(record = T)
plot(st_geometry(shape2))


bf1<- monthly_data %>%
  group_by(year, month, ADMN2) %>%
  summarise(outbreak_occur = as.integer(any(outbreak > 0, na.rm =TRUE)),.groups = 'drop')


#Organises/sorts out area/time covariates
bf1$area <- as.numeric(as.factor(bf1$ADMN2))
bf1$ID.area <- as.numeric(bf1$area)
bf1$ID.area1 <- as.numeric(bf1$area)
bf1$ID.year <- as.numeric(as.factor(bf1$year))
bf1$ID.year1 <- as.numeric(as.factor(bf1$year))
bf1$ID.area.time <- seq(1,length(bf1$ADMN2))
bf1$ID.year.int <- bf1$ID.year
bf1$ID.area.int <- bf1$ID.area
bf1$ID.area.int<-as.numeric(bf1$ID.area.int)


####NEIGHBOURHOOD MATRICES FOR BYM MODEL ====
valid_shapefile<-st_is_valid(shape2)
table(valid_shapefile)
shapefile_spatial <- as_Spatial(shape2)
nb <- poly2nb(shapefile_spatial, queen = TRUE)

Coords <- coordinates(shapefile_spatial)
plot(shapefile_spatial, col = adjustcolor("lightblue", alpha.f = 0.5), border = "darkgray", lwd = 0.7)
plot(nb, coords = Coords, add = TRUE, lwd = 1, col = "blue")

nb2INLA("map.adj", nb)
ken.adj <- paste0(getwd(), "/map.adj")

####CURRENT BEST WORKING SPACE TIME MODEL, WITH NO ENVIROMENTAL COVARIATES ====

form5 <- outbreak_occur ~  f(month, model = 'ar1', cyclic = TRUE) +
  f(ID.area, model = "bym2", graph = ken.adj,  hyper = "shyper",scale.model = TRUE, constr = TRUE)

M5 <- inla(
  form5,
  data = bf1, 
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
)

cat("Summary of the INLA Model for Week (rw1):\n")
cat("\nDIC:", M5$dic$dic, "\n")
cat("WAIC:", M5$waic$waic, "\n")


cat("\nModel Summary:\n")
summary(M5)

bf1$M5 <- M5$summary.fitted.values[,1]
cpo_valuesM5<-M5$cpo$cpo
cpo_valuesM5<-mean(log(cpo_valuesM5))
cpo_valuesM5



###JOIN OUTBREAK AND ENVIRONMENTAL DATA ====

bf1$country <- "Burkina Faso"
bf1$district_country <- paste(bf1$ADMN2, bf1$country, sep = " ")
bf1$year_month <- paste0(bf1$year, " ", sprintf("%02d", bf1$month))
bf1$code <- paste(bf1$year_month, bf1$district_country, sep = " ")

data$month <- as.numeric(data$month)
data$code <- paste(data$year_month, data$district_country, sep = " ")

bf3 <- merge(bf1, data, by = "code")

###WEIGHTING STRATEGY 1 ====
test<-table(bf3$outbreak_occur.y)
test<-as.data.frame(test)
outbreak<-test[2,2]
no_outbreak<-test[1,2]
outbreak<-as.numeric(outbreak)
no_outbreak<-as.numeric(no_outbreak)
total<-no_outbreak+outbreak

non_outbreak_weight <- outbreak/total
outbreak_weight<-1-non_outbreak_weight

weights <- ifelse(bf3$outbreak_occur.y == 1, outbreak_weight, non_outbreak_weight) 

bf3$humidity_std <- scale(bf3$humidity)
bf3$rainfall_std <- scale(bf3$rainfall)
bf3$windspeed_std <- scale(bf3$windspeed)
bf3$north_wind_std <- scale(bf3$north_wind)
bf3$eastward_wind_std <- scale(bf3$eastward_wind)
bf3$aod_std <- scale(bf3$aod)
bf3$temp_std <- scale(bf3$temp)
###FULL MODEL, BEFORE TAKING OUT VARIABLES ====

form5 <- outbreak_occur.y ~ aod_std + vaccine + temp_std + humidity_std + total_cropland +
  f(month.x, model = "ar1", cyclic = TRUE)
  f(area, model = "bym2", graph = ken.adj,  hyper = "shyper",scale.model = TRUE, constr = TRUE)

M5 <- inla(
  form5,
  data = bf3,  # Use your aggregated monthly dataset
  family = "binomial",  # Binary outbreak occurrence
  weights=weights,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
)

cat("Summary of the INLA Model for Week (rw1):\n")
cat("\nDIC:", M5$dic$dic, "\n")
cat("WAIC:", M5$waic$waic, "\n")


cat("\nModel Summary:\n")
summary(M5)

bf3$M5 <- M5$summary.fitted.values[,1]
cpo_valuesM5<-M5$cpo$cpo
cpo_valuesM5<-mean(log(cpo_valuesM5))
cpo_valuesM5


predicted_probs<- M5$summary.fitted.values[,1]
test <- as.data.frame(
  sensitivity_specificity(
    predicted_probs = predicted_probs,
    actual = bf3$outbreak_occur.y))


specificity<-test[1,2]
sensitivity<-test[1,1]
# Create a data frame with ORs and 95% credible intervals
fixed_effects <- M5$summary.fixed
ORs <- data.frame(
  Variable = rownames(fixed_effects),
  OR = exp(fixed_effects$mean),
  Lower = exp(fixed_effects$`0.025quant`),
  Upper = exp(fixed_effects$`0.975quant`)
)



library(ggplot2)

# Reorder factors for plotting
ORs$Variable <- factor(ORs$Variable, levels = ORs$Variable[order(ORs$OR, decreasing = TRUE)])

# Create the plot
test<-ggplot(ORs, aes(x = Variable, y = OR)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +  # Flip coordinates for better readability
  labs(
    title = "Odds Ratios with 95% Credible Intervals",
    y = "Odds Ratio (logit scale exp(beta))",
    x = "Variable"
  ) +
  theme_minimal(base_size = 14)

plot(test)



form5 <- outbreak_occur.y ~  aod_std + vaccine + temp_std + humidity_std  +
  f(month.x, model = "ar1", cyclic = TRUE) + 
  f(area, model = "bym2", graph = ken.adj,  hyper = "shyper",scale.model = TRUE, constr = TRUE)

M5 <- inla(
  form5,
  data = bf3,  # Use your aggregated monthly dataset
  family = "binomial",  # Binary outbreak occurrence
  weights=weights,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
)

cat("Summary of the INLA Model for Week (rw1):\n")
cat("\nDIC:", M5$dic$dic, "\n")
cat("WAIC:", M5$waic$waic, "\n")


cat("\nModel Summary:\n")
summary(M5)

bf3$M5 <- M5$summary.fitted.values[,1]
cpo_valuesM5<-M5$cpo$cpo
cpo_valuesM5<-mean(log(cpo_valuesM5))
cpo_valuesM5


predicted_probs<- M5$summary.fitted.values[,1]
test <- as.data.frame(
  sensitivity_specificity(
    predicted_probs = predicted_probs,
    actual = bf3$outbreak_occur.y))


specificity<-test[1,2]
sensitivity<-test[1,1]
# Create a data frame with ORs and 95% credible intervals
fixed_effects <- M5$summary.fixed
ORs <- data.frame(
  Variable = rownames(fixed_effects),
  OR = exp(fixed_effects$mean),
  Lower = exp(fixed_effects$`0.025quant`),
  Upper = exp(fixed_effects$`0.975quant`)
)



library(ggplot2)

# Reorder factors for plotting
ORs$Variable <- factor(ORs$Variable, levels = ORs$Variable[order(ORs$OR, decreasing = TRUE)])

# Create the plot
test<-ggplot(ORs, aes(x = Variable, y = OR)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +  # Flip coordinates for better readability
  labs(
    title = "Odds Ratios with 95% Credible Intervals",
    y = "Odds Ratio (logit scale exp(beta))",
    x = "Variable"
  ) +
  theme_minimal(base_size = 14)

plot(test)


library(INLA)
library(dplyr)
library(INLA)
library(dplyr)
library(caret)

  
  # Unique districts to loop over
  districts <- unique(bf3$ADMN2)
  
  # Standardize covariates
  bf3 <- bf3 %>%
    mutate(
      humidity_std = scale(humidity),
      rainfall_std = scale(rainfall),
      windspeed_std = scale(windspeed),
      north_wind_std = scale(north_wind),
      eastward_wind_std = scale(eastward_wind),
      aod_std = scale(aod),
      temp_std = scale(temp)
    )
  
  # Weighting strategy
  test <- table(bf3$outbreak_occur.y)
  outbreak <- as.numeric(test["1"])
  no_outbreak <- as.numeric(test["0"])
  total <- outbreak + no_outbreak
  non_outbreak_weight <- outbreak / total
  outbreak_weight <- 1 - non_outbreak_weight
  bf3$weights <- ifelse(bf3$outbreak_occur.y == 1, outbreak_weight, non_outbreak_weight)
  # Store predictions and confusion matrices
  all_predictions <- data.frame()
  confusion_matrices <- list()

  

  # Vector of districts you want to iterate over
  districts_to_test <- c("Boulgou", "Boulkiemdé", "Houet","Kompienga",
                         "Loroum","Sanmatenga", "Yatenga"
                         )
  
  # Initialize storage for predictions and confusion matrices
  all_predictions <- data.frame()
  confusion_matrices <- list()
  
  # Loop over each district
  for (leave_out_ADMN2 in districts_to_test) {
    message("Leaving out: ", leave_out_ADMN2)
    
    data_copy <- bf3
    
    # Identify test rows
    test_indices <- which(data_copy$ADMN2 == leave_out_ADMN2)
    
    # Mask test district outcomes
    data_copy$outbreak_occur.y[test_indices] <- NA
    
    # Fit INLA model
    M_loocv <- inla(
      outbreak_occur.y ~ aod_std + temp_std + vaccine + humidity_std  +
        f(month.x, model = "ar1", cyclic = TRUE) +
        f(area, model = "bym2", graph = ken.adj,
          
          scale.model = TRUE,
          constr = TRUE
        ),
      data = data_copy,
      family = "binomial",
      weights= weights,
      control.fixed = list(mean = 0, prec = 0.001),
      control.predictor = list(compute = TRUE, link = 1),
      control.compute = list(dic = FALSE, waic = FALSE)
    )
    
    # Predict on left-out district only
    predicted_probs <- M_loocv$summary.fitted.values$mean[test_indices]
    actual_vals <- bf3$outbreak_occur.y[test_indices]
    predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
    
    # Store predictions
    district_preds <- data.frame(
      ADMN2 = leave_out_ADMN2,
      actual = actual_vals,
      predicted_prob = predicted_probs,
      predicted_class = predicted_classes
    )
    
    all_predictions <- bind_rows(all_predictions, district_preds)
    
    # Store confusion matrix
    cm_table <- table(Actual = actual_vals, Predicted = predicted_classes)
    confusion_matrices[[leave_out_ADMN2]] <- cm_table
  }
  
  all_predictions
  
  
  overall_conf_mat <- table(
    Actual = all_predictions$actual,
    Predicted = all_predictions$predicted_class
  )
  print(overall_conf_mat)
  
  result = cor.test( all_predictions$predicted_prob,  all_predictions$actual, method = "pearson")
  
  result

  

  
  # Replace with your actual data
  actual <- all_predictions$actual
  predicted_prob <- all_predictions$predicted_prob
  predicted_class <- all_predictions$predicted_classes
  # 1. Brier Score
  brier <- BrierScore(predicted_classes, actual)
  cat("Brier Score:", brier, "\n")
  
  # 2. ROC and AUC
  roc_curve <- roc(actual,predicted_prob, auc=TRUE)
  # Plot ROC curve
  plot(roc_curve, col = "blue", main = "ROC Curve", print.auc = TRUE)
  
  # Add additional information to the plot
  abline(a = 0, b = 1, lty = 2, col = "red")  # Diagonal line for random guess

  model_formulas <- list(
    model1 = outbreak_occur.y ~ aod_std + temp_std + vaccine + humidity_std  +
      f(month.x, model = "rw1", cyclic = TRUE) +
      f(area, model = "bym2", graph = ken.adj, scale.model = TRUE, constr = TRUE),
    
model2 = outbreak_occur.y ~ aod_std + temp_std + vaccine + humidity_std  +
      f(area, model = "bym2", graph = ken.adj, scale.model = TRUE, constr = TRUE) +
      f(month.x, model = "rw2", cyclic = TRUE) ,
    model4 = outbreak_occur.y ~ aod_std + temp_std + vaccine + humidity_std  +
      f(area, model = "bym2", graph = ken.adj, scale.model = TRUE, constr = TRUE) +
      f(month.x, model = "ar1", cyclic = TRUE),
model5 = outbreak_occur.y ~ aod_std * temp_std + vaccine + humidity_std  +
  f(area, model = "bym2", graph = ken.adj, scale.model = TRUE, constr = TRUE) +
  f(month.x, model = "ar1", cyclic = TRUE)
  )
  library(dplyr)
  
  run_loocv_models <- function(data, districts, formulas, weights, adjacency_graph) {
    all_results <- list()
    
    for (model_name in names(formulas)) {
      message("Running model: ", model_name)
      
      model_preds <- data.frame()
      waic_values <- numeric(length(districts))
      
      for (i in seq_along(districts)) {
        district <- districts[i]
        data_copy <- data
        test_indices <- which(data_copy$ADMN2 == district)
        data_copy$outbreak_occur.y[test_indices] <- NA
        
        model_fit <- inla(
          formulas[[model_name]],
          data = data_copy,
          family = "binomial",
          weights = weights,
          control.fixed = list(mean = 0, prec = 0.01),
          control.predictor = list(compute = TRUE, link = 1),
          control.compute = list(waic = TRUE, dic = FALSE)
        )
        
        predicted_probs <- model_fit$summary.fitted.values$mean[test_indices]
        actual_vals <- data$outbreak_occur.y[test_indices]
        predicted_class <- ifelse(predicted_probs > 0.5, 1, 0)
        
        model_preds <- bind_rows(model_preds,
                                 data.frame(
                                   ADMN2 = district,
                                   actual = actual_vals,
                                   predicted_prob = predicted_probs,
                                   predicted_class = predicted_class
                                 ))
        
        waic_values[i] <- model_fit$waic$waic
      }
      
      overall_conf_mat <- table(
        Actual = model_preds$actual,
        Predicted = model_preds$predicted_class
      )
      
      pearson_corr <- cor.test(model_preds$predicted_prob, model_preds$actual, method = "pearson")
      
      all_results[[model_name]] <- list(
        predictions = model_preds,
        confusion_matrix = overall_conf_mat,
        waic = mean(waic_values),  # average WAIC over all left-out districts
        pearson_correlation = pearson_corr
      )
    }
    
    return(all_results)
  }
  results <- run_loocv_models(bf3, districts_to_test, model_formulas, weights, ken.adj)
  
  # Compare WAIC for each model
  waic_summary <- sapply(results, function(x) x$waic)
  print(waic_summary)
  
  # Show confusion matrices and correlation for each model
  for (model_name in names(results)) {
    cat("\nModel:", model_name, "\n")
    print(results[[model_name]]$confusion_matrix)
    cat("Pearson correlation (predicted_prob vs actual):\n")
    print(results[[model_name]]$pearson_correlation)
  }
  
  
  
  
  
  # Vector of districts to iterate over
  districts_to_test <- c("Boulgou", "Boulkiemdé", "Houet", "Kompienga",
                         "Loroum", "Sanmatenga", "Yatenga")
  
  # Initialize storage for predictions and confusion matrices
  all_predictions <- data.frame()
  confusion_matrices <- list()
  
  # Loop over each district
  for (leave_out_ADMN2 in districts_to_test) {
    message("Leaving out: ", leave_out_ADMN2)
    
    data_copy <- bf3
    
    # Identify test rows
    test_indices <- which(data_copy$ADMN2 == leave_out_ADMN2)
    
    # Mask test district outcomes
    data_copy$outbreak_occur.y[test_indices] <- NA
    
    # Fit INLA model with informative priors
    M_loocv <- inla(
      outbreak_occur.y ~ aod_std * temp_std + vaccine + humidity_std +
        f(area, model = "bym2", graph = ken.adj, scale.model = TRUE, constr = TRUE,
          hyper = list(
            prec = list(prior = "pc.prec", param = c(0.5, 0.01)),  # SD ~ 0.5 with P(SD > 0.5) = 0.01
            phi = list(prior = "pc", param = c(0.5, 0.5))           # Moderate shrinkage toward no spatial structure
          )
        ) +
        f(month.x, model = "ar1", cyclic = TRUE,
          hyper = list(
            theta = list(prior = "pc.cor1", param = c(0.7, 0.7))  # Shrink AR1 correlation toward 0
          )
        ),
      data = data_copy,
      family = "binomial",
      weights = weights,
      control.fixed = list(
        mean = 0,
        prec = 1   # SD = 1 for all fixed effects
      ),
      control.predictor = list(compute = TRUE, link = 1),
      control.compute = list(dic = FALSE, waic = FALSE)
    )
    
    # Predict on left-out district only
    predicted_probs <- M_loocv$summary.fitted.values$mean[test_indices]
    actual_vals <- bf3$outbreak_occur.y[test_indices]
    predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
    
    # Store predictions
    district_preds <- data.frame(
      ADMN2 = leave_out_ADMN2,
      actual = actual_vals,
      predicted_prob = predicted_probs,
      predicted_class = predicted_classes
    )
    
    all_predictions <- bind_rows(all_predictions, district_preds)
    
    # Store confusion matrix
    cm_table <- table(Actual = actual_vals, Predicted = predicted_classes)
    confusion_matrices[[leave_out_ADMN2]] <- cm_table
  }
  
  # Final outputs
  print("All predictions:")
  print(all_predictions)
  
  # Overall confusion matrix
  overall_conf_mat <- table(
    Actual = all_predictions$actual,
    Predicted = all_predictions$predicted_class
  )
  print("Overall confusion matrix:")
  print(overall_conf_mat)
  
  # Correlation between predicted probabilities and actual values
  result <- cor.test(all_predictions$predicted_prob, all_predictions$actual, method = "pearson")
  print(result)
  library(ggplot2)
  ggplot(all_predictions, aes(x = predicted_prob, y = actual)) +
    geom_jitter(height = 0.1, width = 0) +
    geom_smooth(method = "loess") +
    labs(title = "Calibration Plot", x = "Predicted Probability", y = "Actual Outcome")
  
  
  # ---- RMSE calculation ----
  rmse <- sqrt(mean((all_predictions$predicted_prob - all_predictions$actual)^2, na.rm = TRUE))
  cat("RMSE:", round(rmse, 4), "\n")
  
  
