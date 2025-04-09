
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
weekly_data <- read_excel("BF_outbreak.xlsx")
view(weekly_data)
table(weekly_data$outbreak)
setwd("C:/Users/mvc32/OneDrive - University of Cambridge/Documents/Climate_meningitis_belt")
shape2 <- st_read("Burkinafaso.shp")
windows(record = T)
plot(st_geometry(shape2))




weekly_data$area <- as.numeric(as.factor(weekly_data$ADMN2))
weekly_data$ID.area <- as.numeric(weekly_data$area)
weekly_data$ID.area1 <- as.numeric(weekly_data$area)

weekly_data$ID.year <- as.numeric(as.factor(weekly_data$year))
weekly_data$ID.year1 <- as.numeric(as.factor(weekly_data$year))
weekly_data$ID.area.time <- seq(1,length(weekly_data$ADMN2))
weekly_data$ID.year.int <- weekly_data$ID.year
weekly_data$ID.area.int <- weekly_data$ID.area

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


shyper = list(prec = list(prior = "pcprec", param=c(1, 0.0001)))
shyper = list(prec = list(prior = "pcprec", param=c(1, 0.0001)))
#### MODEL 1, FIXED MONTH AND AR1 WEEK ====
form1 <- outbreak ~ month + f(week, model = "ar1", cyclic = TRUE) +
  f(ID.area, model = "bym2", graph = ken.adj,  hyper = "shyper",scale.model = TRUE, constr = TRUE)  

M1 <- inla(
  form1,
  data = weekly_data,
  family = "binomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)


cat("\nModel Summary:\n")
summary(M1)

weekly_data$M1 <- M1$summary.fitted.values[,1]
cpo_valuescyclic_ar1<-M1$cpo$cpo
lcpocpo_valueM1<--mean(log(cpo_valuescyclic_ar1))
cat("CPO value:")
print(lcpocpo_valueM1)

sensitivity_specificity <- function(predicted_probs, actual) {
  predicted_classes <- ifelse(predicted_probs > 0.4, 1, 0)
  
  conf_matrix <- table(Actual = actual, Predicted = predicted_classes)
  print("Confusion Matrix:")
  print(conf_matrix)
  predicted_classes <- as.factor(predicted_classes)
  actual <- as.factor(actual)
  
  sensitivity_value <- sensitivity(predicted_classes, actual, positive = "1")
  specificity_value <- specificity(predicted_classes, actual, negative = "0")
  
  cat("Sensitivity:", sensitivity_value, "\n")
  cat("Specificity:", specificity_value, "\n")
}


weekly_data$outbreak <- as.numeric(as.character(weekly_data$outbreak))
weekly_data$M1 <- as.numeric(as.character(weekly_data$M1))
sensitivity_specificity(weekly_data$M1, weekly_data$outbreak)



####MODEL 2, FUXED YEAR AND MONTH, RANDOM WEEK ====

form2 <- outbreak ~ year + month + f(week, model = "ar1", cyclic = TRUE) +
  f(ID.area, model = "bym2", graph = ken.adj,  hyper = "shyper",scale.model = TRUE, constr = TRUE)  

M2 <- inla(
  form2,
  data = weekly_data,
  family = "binomial",  # Binary outbreak occurrence
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
)



cat("\nModel Summary:\n")
summary(M2)

weekly_data$M2 <- M2$summary.fitted.values[,1]
cpo_valuesM2<-M2$cpo$cpo
cpo_valuesM2<--mean(log(cpo_valuesM2))
cpo_valuesM2

weekly_data$outbreaK <- as.numeric(as.character(weekly_data$outbreak))
weekly_data$M2 <- as.numeric(as.character(weekly_data$M2))

sensitivity_specificity(weekly_data$M2, weekly_data$outbreak)



####MODEL 3, FIXED MONTH, RANDOM WEEK AND YEAR ====

form3 <- outbreak ~ month + f(week, model = "ar1", cyclic = TRUE) +
  f(ID.year, model = "rw2") + 
  f(ID.year1, model = "iid") +
  f(ID.area, model = "bym2", graph = ken.adj,  hyper = "shyper",scale.model = TRUE, constr = TRUE)  

M3 <- inla(
  form3,
  data = weekly_data,
  family = "binomial",  # Binary outbreak occurrence
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
)




cat("\nModel Summary:\n")
summary(M3)

weekly_data$M3 <- M3$summary.fitted.values[,1]
cpo_valuesM3<-M3$cpo$cpo
cpo_valuesM3<--mean(log(cpo_valuesM3))
cpo_valuesM3

weekly_data$outbreaK <- as.numeric(as.character(weekly_data$outbreak))
weekly_data$M3 <- as.numeric(as.character(weekly_data$M3))

sensitivity_specificity(weekly_data$M3, weekly_data$outbreak)




####MODEL 4, RANDOM WEEK AND YEAR ====

form4 <- outbreak ~ year +f(week, model = "ar1", cyclic = TRUE) +
  f(ID.year1, model = "iid") +
  f(ID.area, model = "bym2", graph = ken.adj,  hyper = "shyper",scale.model = TRUE, constr = TRUE)  

M4 <- inla(
  form4,
  data = weekly_data,
  family = "binomial",  # Binary outbreak occurrence
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
)


summary(M4)



form4 <- outbreak ~  f(week, model = "ar1", cyclic = TRUE) +
  f(ID.year1, model = "iid") +
  f(ID.area, model = "bym2", graph = ken.adj,  hyper = "shyper",scale.model = TRUE, constr = TRUE)  

M4 <- inla(
  form4,
  data = weekly_data,
  family = "binomial",  # Binary outbreak occurrence
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
)


summary(M4)

weekly_data$M4 <- M4$summary.fitted.values[,1]
cpo_valuesM4<-M4$cpo$cpo
cpo_valuesM4<--mean(log(cpo_valuesM4))
cpo_valuesM4

weekly_data$outbreaK <- as.numeric(as.character(weekly_data$outbreak))
weekly_data$M4 <- as.numeric(as.character(weekly_data$M4))

sensitivity_specificity(weekly_data$M4, weekly_data$outbreak)




####MODEL 5, RANDOM WEEK AND YEAR AND MONTH ====

form5 <- outbreak ~ year+ month +f(week, model = "ar1", cyclic = TRUE) +
  f(ID.year, model = "rw2") + 
  f(ID.year1, model = "iid") +
  f(ID.area, model = "bym2", graph = ken.adj,  hyper = "shyper",scale.model = TRUE, constr = TRUE)  

M5 <- inla(
  form5,
  data = weekly_data,
  family = "binomial",  # Binary outbreak occurrence
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
)


summary(M5)

weekly_data$M5 <- M5$summary.fitted.values[,1]
cpo_valuesM5<-M5$cpo$cpo
cpo_valuesM5<--mean(log(cpo_valuesM5))
cpo_valuesM5

weekly_data$outbreaK <- as.numeric(as.character(weekly_data$outbreak))
weekly_data$M5 <- as.numeric(as.character(weekly_data$M5))

sensitivity_specificity(weekly_data$M5, weekly_data$outbreak)




######## DIC, WAIC VALUES


dic <- data.frame(criteria = rep(c("DIC", "WAIC")), M1= rep(c(M1$dic$dic,M1$waic$waic)),
                  M2= rep(c(M2$dic$dic,M2$waic$waic)),
                  M3= rep(c(M3$dic$dic,M3$waic$waic)),
                  M4= rep(c(M4$dic$dic,M4$waic$waic)),
                  M5= rep(c(M5$dic$dic,M5$waic$waic)) )



dic
