
library(INLA)
library(ggplot2)
library(mgcv)
library(dplyr)
library(DAAG)
library(reshape2)
library(caret)
library(sf)
library(spdep)
library(readxl)
library(gridExtra)
library(haven)
library(tidyverse)
library(tmap)
library(pROC)


setwd("~/INLA project")
monthly_data <- read_excel("BF_outbreak.xlsx")
view(monthly_data)
table(monthly_data$outbreak)
setwd("C:/Users/mvc32/OneDrive - University of Cambridge/Documents/Climate_meningitis_belt")
shape2 <- st_read("Burkinafaso.shp")
windows(record = T)
plot(st_geometry(shape2))


#Monthly model
bf1<- monthly_data %>%
  group_by(year, month, ADMN2) %>%
  summarise(outbreak_occur = as.integer(any(outbreak > 0, na.rm =TRUE)),.groups = 'drop')
print(head(monthly_data))

bf1$area <- as.numeric(as.factor(bf1$ADMN2))
bf1$ID.area <- as.numeric(bf1$area)
bf1$ID.area1 <- as.numeric(bf1$area)
bf1$ID.year <- as.numeric(as.factor(bf1$year))
bf1$ID.year1 <- as.numeric(as.factor(bf1$year))
bf1$ID.area.time <- seq(1,length(bf1$ADMN2))
bf1$ID.year.int <- bf1$ID.year
bf1$ID.area.int <- bf1$ID.area

#### SIMPLE TEMPORAL MODEL ====
# no prior assumptions on data
form1 <- outbreak_occur ~ year+ f(month, model = "ar1", cyclic = TRUE)  

m1 <- inla(
  form1,
  data = bf1,
  family = "binomial",  # no zero inflation here, look at and add in
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)


####NEIGHBOURHOOD MATRICES ====
valid_shapefile<-st_is_valid(shape2)
table(valid_shapefile)
### Make spatial neighborhood structure
shapefile_spatial <- as_Spatial(shape2)
nb <- poly2nb(shapefile_spatial, queen = TRUE)


Coords <- coordinates(shapefile_spatial)
plot(shapefile_spatial, col = adjustcolor("lightblue", alpha.f = 0.5), border = "darkgray", lwd = 0.7)
plot(nb, coords = Coords, add = TRUE, lwd = 1, col = "blue")


# Save adjacency matrix for INLA
nb2INLA("map.adj", nb)
ken.adj <- paste0(getwd(), "/map.adj")
# nb2 <- inla.read.graph(filename = "map.adj")
# inla.debug.graph("map.adj")
# 
# W.boston <- nb2mat(nb, style = "B", zero.policy = TRUE)
# W.boston.rs <- nb2mat(nb, style = "W", zero.policy = TRUE)

shyper = list(prec = list(prior = "pcprec", param=c(1, 0.0001)))
####SPATIOTEMPORAL MODEL FIXED YEAR ====

form2 <- outbreak_occur ~ year +
  f(month, model = "ar1", cyclic = TRUE) +  # Cyclic AR1 for month
  f(ID.area, model = "bym2", graph = ken.adj, hyper = "shyper", scale.model = TRUE) 

M2 <- inla(
  form2,
  data = bf1,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

cat("Summary of the INLA Model for Week (rw1):\n")
cat("\nDIC:", M2$dic$dic, "\n")
cat("WAIC:", M2$waic$waic, "\n")

cat("\nModel Summary:\n")
summary(M2)

bf1$M2 <- M2$summary.fitted.values[,1]
cpo_valuescyclic_ar1<-M2$cpo$cpo
lcpocpo_valueM2<--mean(log(cpo_valuescyclic_ar1))
cat("CPO value:")
print(lcpocpo_valueM2)

sensitivity_specificity <- function(predicted_probs, actual) {
  predicted_classes <- ifelse(predicted_probs > 0.4, 1, 0)
  
  conf_matrix <- table(Actual = actual, Predicted = predicted_classes)
  print("Confusion Matrix:")
  print(conf_matrix)
  
  # Convert to factors for `sensitivity()` and `specificity()`
  predicted_classes <- as.factor(predicted_classes)
  actual <- as.factor(actual)
  
  # Calculate sensitivity and specificity
  sensitivity_value <- sensitivity(predicted_classes, actual, positive = "1")
  specificity_value <- specificity(predicted_classes, actual, negative = "0")
  
  # Print results
  cat("Sensitivity:", sensitivity_value, "\n")
  cat("Specificity:", specificity_value, "\n")
}


bf1$outbreak_occur <- as.numeric(as.character(bf1$outbreak_occur))
bf1$M2 <- as.numeric(as.character(bf1$M2))
sensitivity_specificity(bf1$M2, bf1$outbreak_occur)


sensitivity_specificity_dynamic <- function(predicted_probs, actual) {
  # Ensure input is numeric
  predicted_probs <- as.numeric(predicted_probs)
  actual <- as.numeric(actual)
  
  # Compute ROC curve
  roc_curve <- roc(actual, predicted_probs)
  
  # Find the optimal threshold using the Youden Index
  optimal_index <- which.max(roc_curve$sensitivities + roc_curve$specificities - 1)
  optimal_threshold <- roc_curve$thresholds[optimal_index]
  
  # Convert predictions to binary using the optimal threshold
  predicted_classes <- ifelse(predicted_probs > optimal_threshold, 1, 0)
  
  # Create confusion matrix
  conf_matrix <- table(Actual = actual, Predicted = predicted_classes)
  print("Confusion Matrix:")
  print(conf_matrix)
  
  # Convert to factors for caret's sensitivity and specificity functions
  predicted_classes <- as.factor(predicted_classes)
  actual <- as.factor(actual)
  
  # Calculate sensitivity and specificity
  sensitivity_value <- sensitivity(predicted_classes, actual, positive = "1")
  specificity_value <- specificity(predicted_classes, actual, negative = "0")
  
  # Print results
  cat("Optimal Threshold:", optimal_threshold, "\n")
  cat("Sensitivity:", sensitivity_value, "\n")
  cat("Specificity:", specificity_value, "\n")
}

# Run the function
sensitivity_specificity_dynamic(bf1$M2, bf1$outbreak_occur)


####SPATIOTEMPORAL MODEL PIECEWISE REGRESSION ====


bf1$Year_Centered <- bf1$year - 2010
bf1$After2009 <- ifelse(bf1$year > 2010, bf1$Year_Centered, 0)

form3 <- outbreak_occur ~ Year_Centered + After2009 +
  f(month, model = "ar1", cyclic = TRUE) +
  f(ID.area, model = "bym2", graph = ken.adj,  hyper = "shyper",scale.model = TRUE, constr = TRUE)

  M3 <- inla(
    form3,
    data = bf1,  # Use your aggregated monthly dataset
    family = "binomial",  # Binary outbreak occurrence
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
  )
  
  cat("Summary of the INLA Model for Week (rw1):\n")
  cat("\nDIC:", M3$dic$dic, "\n")
  cat("WAIC:", M3$waic$waic, "\n")
  
  cat("\nModel Summary:\n")
  summary(M3)
  
  bf1$M3 <- M3$summary.fitted.values[,1]
  cpo_valuesM3<-M3$cpo$cpo
  cpo_valuesM3<--mean(log(cpo_valuesM3))
  
  cpo_valuesM3
bf1$outbreak_occur <- as.numeric(as.character(bf1$outbreak_occur))
  bf1$M3 <- as.numeric(as.character(bf1$M3))
  
  sensitivity_specificity(bf1$M3, bf1$outbreak_occur)
  
  sensitivity_specificity_dynamic(bf1$M3, bf1$outbreak_occur)
  
  
  
  
  ####SPATIOTEMPORAL MODEL ADDED RANDOM YEAR EFFECTS ====

  form4 <- outbreak_occur ~ f(month, model = "ar1", cyclic = TRUE) +
    f(ID.year, model = "rw2") + 
    f(ID.year1, model = "iid") +
    f(ID.area, model = "bym2", graph = ken.adj,  hyper = "shyper",scale.model = TRUE, constr = TRUE)
  
  M4 <- inla(
    form4,
    data = bf1,  # Use your aggregated monthly dataset
    family = "binomial",  # Binary outbreak occurrence
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
  )
  
  cat("Summary of the INLA Model for Week (rw1):\n")
  cat("\nDIC:", M4$dic$dic, "\n")
  cat("WAIC:", M4$waic$waic, "\n")

  
  cat("\nModel Summary:\n")
  summary(M4)
  
  bf1$M4 <- M4$summary.fitted.values[,1]
  cpo_valuesM4<-M4$cpo$cpo
  cpo_valuesM4<--mean(log(cpo_valuesM4))
  cpo_valuesM4
  
  bf1$outbreak_occur <- as.numeric(as.character(bf1$outbreak_occur))
  bf1$M4 <- as.numeric(as.character(bf1$M4))
  
  sensitivity_specificity(bf1$M4, bf1$outbreak_occur)
  
  sensitivity_specificity_dynamic(bf1$M4, bf1$outbreak_occur)
  
  
  
  
  
  form5 <- outbreak_occur ~ f(month, model = "ar1", cyclic = TRUE) +
    f(ID.year1, model = "iid") +
    f(ID.area, model = "bym2", graph = ken.adj,  hyper = "shyper",scale.model = TRUE, constr = TRUE)
  
  M5 <- inla(
    form5,
    data = bf1,  # Use your aggregated monthly dataset
    family = "binomial",  # Binary outbreak occurrence
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
  )
  
  cat("Summary of the INLA Model for Week (rw1):\n")
  cat("\nDIC:", M5$dic$dic, "\n")
  cat("WAIC:", M5$waic$waic, "\n")
  
  
  cat("\nModel Summary:\n")
  summary(M4)
  
  bf1$M5 <- M5$summary.fitted.values[,1]
  cpo_valuesM5<-M5$cpo$cpo
  cpo_valuesM5<--mean(log(cpo_valuesM5))
  cpo_valuesM5
  
  bf1$outbreak_occur <- as.numeric(as.character(bf1$outbreak_occur))
  bf1$M5 <- as.numeric(as.character(bf1$M5))
  
  sensitivity_specificity(bf1$M5, bf1$outbreak_occur)
  
  sensitivity_specificity_dynamic(bf1$M5, bf1$outbreak_occur)
  
  ######## DIC, WAIC, LOGCPO VALUES
  

  dic <- data.frame(criteria = rep(c("DIC", "WAIC")), M1= rep(c(m1$dic$dic,m1$waic$waic)),
                    M2= rep(c(M2$dic$dic,M2$waic$waic)),
                    M3= rep(c(M3$dic$dic,M3$waic$waic)),
                    M4= rep(c(M4$dic$dic,M4$waic$waic)), 
                    M5= rep(c(M5$dic$dic,M5$waic$waic)))
  
  
  dic

   ####PLOTTING MONTHLY OUTBREAK PROBABILITY OF BEST MODEL ====
  
  merged_data$year_month <- paste(merged_data$year, merged_data$month, sep = "_")

  split_df <- split(merged_data, merged_data$year_month)
  
  for (year_month in names(split_df)) {
    shapefile_sf <- st_as_sf(split_df[[year_month]])
    assign(paste0("shapefile_", year_month), shapefile_sf, envir = .GlobalEnv)
  }
  shapefiles_list <- mget(ls(pattern = "shapefile_\\d{4}_\\d{1,2}"))

  plot_list <- list()
  
  
  for (year_month in names(shapefiles_list)) {
    current_shapefile <- shapefiles_list[[year_month]]
    current_shapefile$M5 <- as.numeric(current_shapefile$M4)
    
    plot <- ggplot(data = current_shapefile) +
      geom_sf(aes(fill = M5), color = "black", size = 0.1) +
      scale_fill_gradient(low = "white", high = "red", na.value = "gray") + 

      labs(title = paste("Outbreak Probability in", year_month),
           fill = "Probability") +
      
      theme_minimal() +
      theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
    
    plot_list[[year_month]] <- plot
  }
  
  
  plot_list[["shapefile_2007_3"]] 
