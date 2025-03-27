
library(INLA)
library(ggplot2)
library(mgcv)
library(dplyr)
library(DAAG)
library(reshape2)
library(sf)
library(spdep)
library(caret)
library(PBSmapping)

# Load and prepare data
setwd("~/INLA project")
spatiotemporaloutbreaks <- read.csv("full-outbreak-matched.csv")


newdata<-spatiotemporaloutbreaks %>%
  group_by(country, name_1, name_2,district_country)%>%
  summarise(outbreak=ifelse(sum(outbreak)>0,1,0), .groups="drop")

newdata$outbreak <- as.numeric(newdata$outbreak)
zero_proportion <- sum(newdata$outbreak == 0) / nrow(newdata)
zero_proportion 


setwd("C:/Users/mvc32/OneDrive - University of Cambridge/Documents/Climate_meningitis_belt")
shape2 <- st_read("Shapefile_improved.shp")
countries_to_keep<-"Burkina Faso"
shape2<-shape2%>%
  filter(COUNTRY %in% countries_to_keep)
newdata<-newdata%>%
  filter(country %in% countries_to_keep)

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

# Read adjacency matrix
nb2 <- inla.read.graph(filename = "map.adj")
inla.debug.graph("map.adj")

W.boston <- nb2mat(nb, style = "B", zero.policy = TRUE)
W.boston.rs <- nb2mat(nb, style = "W", zero.policy = TRUE)


newdata$ID <- seq_len(nrow(newdata))
newdata$ID <- as.numeric(factor(newdata$ID, levels = 1:nrow(W.boston)))
newdata <- newdata[order(newdata$ID), ]


### **Define Function to Fit Different Models**
fit_inla_model <- function(formula, model_name) {
  model <- inla(
    formula,
    data = newdata,
    family = "binomial",  # Changed to binomial for binary outcomes
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
  cat("\nModel:", model_name)
  print(summary(model))
  return(model)
}



### **Run Different Models**
# IID (Independent Effects)
formula_iid <- outbreak ~ f(ID, model = "iid")
model_iid <- fit_inla_model(formula_iid, "IID")
cpo_valuesiid<-model_iid$cpo$cpo
lcpoiid<--mean(log(cpo_valuesiid))

# BYM2 (Spatial + Unstructured Random Effects)
formula_bym2 <- outbreak ~ f(ID, model = "bym2", graph = ken.adj, scale.model = TRUE, constr = TRUE)
model_bym2 <- fit_inla_model(formula_bym2, "BYM2")

cpo_valuesbym2<-model_bym2$cpo$cpo
lcpobym2<--mean(log(cpo_valuesbym2))

# BYM (Spatial + Unstructured Random Effects)
formula_bym <- outbreak ~ f(ID, model = "bym", graph = ken.adj, scale.model = TRUE, constr = TRUE)
model_bym <- fit_inla_model(formula_bym2, "BYM")

cpo_valuesbym<-model_bym$cpo$cpo
lcpobym<--mean(log(cpo_valuesbym))

# BESAG (Spatial Only)
formula_besag <- outbreak ~ f(ID, model = "besag", graph = ken.adj, scale.model = TRUE)
model_besag <- fit_inla_model(formula_besag, "BESAG")

cpo_valuesbesag<-model_besag$cpo$cpo
lcpobesag<--mean(log(cpo_valuesbesag))

# BESAGPROPER (Improved Spatial Model)
formula_besagprop <- outbreak ~ f(ID, model = "besagproper", graph = ken.adj)
model_besagprop <- fit_inla_model(formula_besagprop, "BESAGPROPER")


cpo_valuesbesagproper<-model_besagprop$cpo$cpo
lcpobesagproper<--mean(log(cpo_valuesbesagproper))


### **Compare Model Performance**
model_comparison <- data.frame(
  Model = c("IID", "BYM2", "BYM",  "BESAG", "BESAGPROPER"),
  DIC = c(model_iid$dic$dic, model_bym2$dic$dic,  model_bym$dic$dic,model_besag$dic$dic, model_besagprop$dic$dic),
  WAIC = c(model_iid$waic$waic, model_bym2$waic$waic,  model_bym$dic$dic,model_besag$waic$waic, model_besagprop$waic$waic),
  LOGCPO = c(lcpoiid, lcpobym2, lcpobym, lcpobesag,lcpobesagproper )
  
  )

print(model_comparison)

### **Visualize Model Predictions**
newdata$P_IID <- model_iid$summary.fitted.values[,1]
newdata$P_BYM2 <- model_bym2$summary.fitted.values[,1]
newdata$P_BESAG <- model_besag$summary.fitted.values[,1]
newdata$P_BESAGPROP <- model_besagprop$summary.fitted.values[,1]
newdata$P_BYM <- model_bym$summary.fitted.values[,1]



### MODEL SENSITIVITY AND SPECIFICITY
sensitivity_specificity <- function(x, y) {
  predicted_classes <- ifelse(x > 0.5, 1, 0)
  conf_matrix <- table(Actual = y, Predicted = predicted_classes)
  print("Confusion Matrix:")
  print(conf_matrix)
  predicted_classes<-as.factor(predicted_classes)
  newdata$outbreak<-as.factor(y)
  # Calculate sensitivity and specificity
  sensitivity<- sensitivity(y, predicted_classes)
  #calculate specificity
  specificity<-  specificity(y, predicted_classes)
  
  # Return sensitivity and specificity
  cat("Sensitivity:", sensitivity)
  
  
  cat(" Specificity:", specificity)
  
}
newdata$outbreak<-as.factor(newdata$outbreak)

sensitivity_specificity(newdata$P_IID,newdata$outbreak)
sensitivity_specificity(newdata$P_BYM,newdata$outbreak)
sensitivity_specificity(newdata$P_BYM2,newdata$outbreak)
sensitivity_specificity(newdata$P_BESAG,newdata$outbreak)
sensitivity_specificity(newdata$P_BESAGPROP,newdata$outbreak)

