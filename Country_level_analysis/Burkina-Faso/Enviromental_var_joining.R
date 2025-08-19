library(INLA)
library(ggplot2)
library(mgcv)
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
library(raster)
library(terra)
library(stringr)
library(sf)
library(readxl)
library(lubridate)
library(dplyr)

# Set working directory and read data
setwd("~/INLA project")
weekly_data <- read_excel("BF_outbreak.xlsx")

# Check outbreak table
table(weekly_data$outbreak)

# Load shapefile
shape2 <- st_read("Burkinafaso.shp")

# Plot geometry
windows(record = TRUE)
plot(st_geometry(shape2))

# Summarize outbreaks by group
weekly_data <- weekly_data %>%
  group_by(year, month, district_country) %>%
  summarise(outbreak_occur = as.integer(any(outbreak > 0, na.rm = TRUE)), .groups = 'drop')

# Create year_month and code fields
weekly_data$year_month <- paste(weekly_data$year, weekly_data$month, sep = " ")
weekly_data$code <- paste(weekly_data$year_month, weekly_data$district_country, sep = " ")



path <- "~/INLA project/temperatuture"
files <- list.files(path, pattern = "nc4$", full.names = TRUE)

# Filter out duplicate files with (1), (2), etc. in their names
files <- files[!grepl("\\([0-9]+\\)", files)]

# Read as SpatRaster
r <- rast(files)

# Summary of raster layers
summary_df <- data.frame(
  Layer = names(r),
  Min = sapply(1:nlyr(r), function(i) min(minmax(r[[i]])[1], na.rm = TRUE)),
  Max = sapply(1:nlyr(r), function(i) max(minmax(r[[i]])[2], na.rm = TRUE))
)
print(summary_df)
print(names(r))

# Check and align CRS between raster and shapefile
if (!identical(crs(r), crs(vect(shape2)))) {
  message("CRS mismatch - Reprojecting raster to match vector")
  r <- project(r, crs(vect(shape2)))
} else {
  message("CRS aligned")
}

# Crop and mask
r <- mask(crop(r, vect(shape2)), vect(shape2))

# Rename layers with timestamps
start_date <- ymd("2003-01-01") 
date_seq <- seq(start_date, by = "1 month", length.out = nlyr(r))
layer_names <- paste0("temp_", year(date_seq), "_", sprintf("%02d", month(date_seq)))
names(r) <- layer_names
print(names(r))

# Prepare final_data structure
final_data <- weekly_data %>%
  mutate(
    year_month = paste(year, sprintf("%02d", month), sep = " "),
    code = paste(year_month, district_country, sep = " "),
    temp = NA  # Initialize temp column
  )

# Loop through each raster layer and extract temperature data
for (i in 1:nlyr(r)) {
  lyr <- r[[i]]
  lyr_name <- names(r)[i]
  cat("Processing:", lyr_name, "\n")
  
  ym_parts <- str_match(lyr_name, "temp_(\\d{4})_(\\d{2})")
  year_val <- as.integer(ym_parts[, 2])
  month_val <- as.integer(ym_parts[, 3])
  
  extracted <- terra::extract(lyr, vect(shape2), fun = mean, na.rm = TRUE)
  extracted_df <- data.frame(
    shape2,
    temp_temp = extracted[, 2]
  ) %>%
    mutate(
      year = year_val,
      month = month_val,
      year_month = paste(year, sprintf("%02d", month), sep = " "),
      district_country = paste(NAME_2, COUNTRY, sep = " "),
      code = paste(year_month, district_country, sep = " ")
    )
  
  # Merge with final_data and update temp values where missing
  final_data <- merge(final_data, extracted_df[, c("code", "temp_temp")], by = "code", all.x = TRUE)
  # Replace NA temp with extracted values
  final_data$temp[is.na(final_data$temp)] <- final_data$temp_temp[is.na(final_data$temp)]
  
  # Remove temp_temp column for next iteration
  final_data$temp_temp <- NULL
}

# Final preview
print(head(final_data))

# Load required libraries
library(terra)
library(dplyr)
library(stringr)

# Set path and read rasters
path <- "~/INLA project/CHIRPS_data_monthly"
files <- list.files(path, pattern = "tif$", full.names = TRUE)
r_stack <- terra::rast(files)  # Use terra::rast instead of raster::stack

# Summary of raster stack
summary_df <- data.frame(
  Layer = names(r_stack),
  Min = sapply(1:nlyr(r_stack), function(i) min(values(r_stack[[i]]), na.rm = TRUE)),
  Max = sapply(1:nlyr(r_stack), function(i) max(values(r_stack[[i]]), na.rm = TRUE))
)
print(summary_df)

# Convert shape2 to terra vector if not already
shape_vect <- vect(shape2)



# Check and align CRS directly using CRS strings
if (!identical(crs(r_stack), crs(shape_vect))) {
  message("CRS mismatch - Reprojecting raster to match vector")
  r_stack <- terra::project(r_stack, crs(shape_vect))
} else {
  message("CRS aligned")
}

# Crop and mask raster stack
r_stack <- terra::mask(terra::crop(r_stack, shape_vect), shape_vect)

#### RAINFALL DATA PROCESSING ####

# Assuming final_data already exists
final_data$rainfall <- NA
#### RAINFALL DATA PROCESSING ####

# Ensure shape2 and final_data exist and have 'code' column
final_data$rainfall <- NA  # Initialize rainfall column


start_date <- ymd("2003-01-01") 
date_seq <- seq(start_date, by = "1 month", length.out = nlyr(r_stack))
layer_names <- paste0("rainfall_", year(date_seq), "_", sprintf("%02d", month(date_seq)))
names(r_stack) <- layer_names
r<-r_stack
print(names(r_stack))

for (i in 1:nlyr(r)) {
  lyr <- r[[i]]
  lyr_name <- names(r)[i]
  cat("Processing:", lyr_name, "\n")
  ym_parts <- str_match(lyr_name, "rainfall_(\\d{4})_(\\d{2})")
  year <- as.integer(ym_parts[,2])
  month <- as.integer(ym_parts[,3])
  extracted <- data.frame(
    shape2,
    rainfall_temp = terra::extract(lyr, shape2, fun = mean, na.rm = TRUE)[,2]
  )
  
  extracted <- extracted %>%
    mutate(
      year = year,
      month = month,
      year_month = paste(year, sprintf("%02d", month), sep = " "),  # Ensure month is two digits
      district_country = paste(NAME_2, COUNTRY, sep = " "),
      code = paste(year_month, district_country, sep = " ")
    )
  
  final_data <- merge(final_data, extracted[, c("code", "rainfall_temp")], by = "code", all.x = TRUE)
  final_data$rainfall[is.na(final_data$rainfall)] <- final_data$rainfall_temp[is.na(final_data$rainfall)]
  final_data$rainfall_temp <- NULL
}

print(head(final_data))
#### AVERAGE WINDSPEED DATA####

path <- "~/INLA project/windspeed"
files <- list.files(path, pattern = "\\.nc$", full.names = TRUE)

r <- rast(files[1])
print(names(r))
r <- rotate(r)

if (!crs(r) == crs(shape2)) {
  message("CRS mismatch - Reprojecting")
  r <- project(r, crs(shape2))
} else {
  message("CRS aligned")
}

plot(r[[1]], main = paste("Rotated raster"))

r <- mask(crop(r, shape2), shape2)
start_date <- ymd("2003-01-01") 
date_seq <- seq(start_date, by = "1 month", length.out = nlyr(r))
layer_names <- paste0("windspeed_", year(date_seq), "_", sprintf("%02d", month(date_seq)))
names(r) <- layer_names
print(names(r))


start_filter <- ymd("2003-01-01")
end_filter <- ymd("2022-12-31")

keep_indices <- which(date_seq >= start_filter & date_seq <= end_filter)
r <- r[[keep_indices]]
final_data <-   final_data  %>%
  mutate(
    year_month = paste(year, sprintf("%02d", month), sep = " "),  # Ensure month is two digits
    code = paste(year_month, district_country, sep = " "),
    windspeed = NA
  )

for (i in 1:nlyr(r)) {
  lyr <- r[[i]]
  lyr_name <- names(r)[i]
  cat("Processing:", lyr_name, "\n")
  ym_parts <- str_match(lyr_name, "windspeed_(\\d{4})_(\\d{2})")
  year <- as.integer(ym_parts[,2])
  month <- as.integer(ym_parts[,3])
  extracted <- data.frame(
    shape2,
    windspeed_temp = terra::extract(lyr, shape2, fun = mean, na.rm = TRUE)[,2]
  )
  
  extracted <- extracted %>%
    mutate(
      year = year,
      month = month,
      year_month = paste(year, sprintf("%02d", month), sep = " "),  # Ensure month is two digits
      district_country = paste(NAME_2, COUNTRY, sep = " "),
      code = paste(year_month, district_country, sep = " ")
    )
  
  final_data <- merge(final_data, extracted[, c("code", "windspeed_temp")], by = "code", all.x = TRUE)
  final_data$windspeed[is.na(final_data$windspeed)] <- final_data$windspeed_temp[is.na(final_data$windspeed)]
  final_data$windspeed_temp <- NULL
}

print(head(final_data))

#### MERIDONAL AND ZONAL WINDSPEED DATA####
#### EASTWARD WIND DATA####

# Set file path and read in NetCDF files
path <- "~/INLA project/meriodonal_zonal"
files <- list.files(path, pattern = "\\.nc$", full.names = TRUE)

# Load raster and split into two sets
r <- rast(files[1])
r1 <- r[[1:240]]
r2 <- r[[241:480]]

# Initial check
print(names(r))
r1 <- rotate(r1)
plot(r1[[1]])

# Reproject if CRS mismatch
if (!crs(r1) == crs(shape2)) {
  message("CRS mismatch - reprojecting")
  r1 <- project(r1, crs(shape2))
} else {
  message("CRS aligned")
}

plot(r1[[1]])

# Crop and mask raster
r1 <- mask(crop(r1, shape2), shape2)

# Create time-based names for layers
start_date <- ymd("2003-01-01")
date_seq <- seq(start_date, by = "1 month", length.out = nlyr(r1))
layer1_names <- paste0("windspeed_", year(date_seq), "_", sprintf("%02d", month(date_seq)))
names(r1) <- layer1_names

print(names(r1))

# Prepare final_data structure
final_data <- final_data %>%
  mutate(
    year_month = paste(year, sprintf("%02d", month), sep = " "),
    code = paste(year_month, district_country, sep = " "),
    eastward_wind = NA  # Initialize eastward_wind column
  )

# Loop through each raster layer and extract windspeed data
for (i in 1:nlyr(r1)) {
  lyr <- r1[[i]]
  lyr_name <- names(r1)[i]
  cat("Processing:", lyr_name, "\n")
  
  ym_parts <- str_match(lyr_name, "windspeed_(\\d{4})_(\\d{2})")
  year <- as.integer(ym_parts[, 2])
  month <- as.integer(ym_parts[, 3])
  
  extracted <- data.frame(
    shape2,
    windspeed_temp = terra::extract(lyr, shape2, fun = mean, na.rm = TRUE)[, 2]
  ) %>%
    mutate(
      year = year,
      month = month,
      year_month = paste(year, sprintf("%02d", month), sep = " "),
      district_country = paste(NAME_2, COUNTRY, sep = " "),
      code = paste(year_month, district_country, sep = " ")
    )
  
  # Merge with final_data and update eastward_wind where missing
  final_data <- merge(final_data, extracted[, c("code", "windspeed_temp")], by = "code", all.x = TRUE)
  final_data$eastward_wind[is.na(final_data$eastward_wind)] <- final_data$windspeed_temp[is.na(final_data$eastward_wind)]
  final_data$windspeed_temp <- NULL
}

# Final preview
print(head(final_data))

#### NORTHWARD WIND DATA####  

# Rotate raster
r2 <- rotate(r2)
plot(r2[[1]])

# Reproject if CRS mismatch
if (!crs(r2) == crs(shape2)) {
  message("CRS mismatch - reprojecting")
  r2 <- project(r2, crs(shape2))
} else {
  message("CRS aligned")
}


# Crop and mask raster
r2 <- mask(crop(r2, shape2), shape2)

# Create time-based names for layers
start_date <- ymd("2003-01-01")
date_seq <- seq(start_date, by = "1 month", length.out = nlyr(r2))
layer2_names <- paste0("windspeed_", year(date_seq), "_", sprintf("%02d", month(date_seq)))
names(r2) <- layer2_names

print(names(r2))

# Prepare final_data structure
final_data <- final_data %>%
  mutate(
    year_month = paste(year, sprintf("%02d", month), sep = " "),
    code = paste(year_month, district_country, sep = " "),
    north_wind = NA  # Initialize north_wind column
  )

# Loop through each raster layer and extract windspeed data
for (i in 1:nlyr(r2)) {
  lyr <- r2[[i]]
  lyr_name <- names(r2)[i]
  cat("Processing:", lyr_name, "\n")
  
  ym_parts <- str_match(lyr_name, "windspeed_(\\d{4})_(\\d{2})")
  year <- as.integer(ym_parts[, 2])
  month <- as.integer(ym_parts[, 3])
  
  extracted <- data.frame(
    shape2,
    windspeed_temp = terra::extract(lyr, shape2, fun = mean, na.rm = TRUE)[, 2]
  ) %>%
    mutate(
      year = year,
      month = month,
      year_month = paste(year, sprintf("%02d", month), sep = " "),
      district_country = paste(NAME_2, COUNTRY, sep = " "),
      code = paste(year_month, district_country, sep = " ")
    )
  
  # Merge with final_data and update north_wind where missing
  final_data <- merge(final_data, extracted[, c("code", "windspeed_temp")], by = "code", all.x = TRUE)
  final_data$north_wind[is.na(final_data$north_wind)] <- final_data$windspeed_temp[is.na(final_data$north_wind)]
  final_data$windspeed_temp <- NULL
}

# Final preview
print(head(final_data))


#### HUMIDITY DATA####  

# Set file path and read in NetCDF files
path <- "~/INLA project/humidity"
files <- list.files(path, pattern = "\\.nc$", full.names = TRUE)
r <- rast(files[1])

print(names(r))
r <- rotate(r)

r2<- r

# Reproject if CRS mismatch
if (!crs(r2) == crs(shape2)) {
  message("CRS mismatch - reprojecting")
  r2 <- project(r2, crs(shape2))
} else {
  message("CRS aligned")
}

# Crop and mask raster
r2 <- mask(crop(r2, shape2), shape2)
plot(r2[[1]])
# Create time-based names for layers
start_date <- ymd("2003-01-01")
date_seq <- seq(start_date, by = "1 month", length.out = nlyr(r2))
layer2_names <- paste0("humidity_", year(date_seq), "_", sprintf("%02d", month(date_seq)))
names(r2) <- layer2_names

print(names(r2))


final_data <- final_data %>%
  mutate(
    year_month = paste(year, sprintf("%02d", month), sep = " "),
    code = paste(year_month, district_country, sep = " "),
    humidity = NA  
  )

# Loop through each raster layer and extract humidity data
for (i in 1:nlyr(r2)) {
  lyr <- r2[[i]]
  lyr_name <- names(r2)[i]
  cat("Processing:", lyr_name, "\n")
  
  ym_parts <- str_match(lyr_name, "humidity_(\\d{4})_(\\d{2})")
  year <- as.integer(ym_parts[, 2])
  month <- as.integer(ym_parts[, 3])
  
  extracted <- data.frame(
    shape2,
    humidity_temp = terra::extract(lyr, shape2, fun = mean, na.rm = TRUE)[, 2]
  ) %>%
    mutate(
      year = year,
      month = month,
      year_month = paste(year, sprintf("%02d", month), sep = " "),
      district_country = paste(NAME_2, COUNTRY, sep = " "),
      code = paste(year_month, district_country, sep = " ")
    )
  
  # Merge with final_data and update humidity_test2 where missing
  final_data <- merge(final_data, extracted[, c("code", "humidity_temp")], by = "code", all.x = TRUE)
  final_data$humidity[is.na(final_data$humidity)] <- final_data$humidity_temp[is.na(final_data$humidity)]
  final_data$humidity_temp <- NULL
}

# Final preview
print(head(final_data))

#### DUST DATA####  


setwd("~/INLA project/aod")


# List all .nc4 files
nc_files <- list.files(pattern = "\\.nc4$", full.names = TRUE)

# Initialize an empty list to hold the monthly rasters
monthly_rasters <- list()

# Loop through each file
for (file in nc_files) {
  message("Processing: ", basename(file))
  
  r <- rast(file)
  
  # Select only AODANA layers
  aod_layers <- r[[grep("^AODANA", names(r))]]
  
  # Average across the 8 time steps
  monthly_mean <- mean(aod_layers)
  
  # Optional: name it by date (e.g., 2003_01)
  date_str <- gsub(".*\\.(\\d{6})\\.nc4$", "\\1", basename(file))  # e.g., "200301"
  names(monthly_mean) <- paste0("AOD_", date_str)
  
  # Optional: Save individual GeoTIFFs
  writeRaster(monthly_mean, paste0("AOD_", date_str, ".tif"), overwrite = TRUE)
  
  # Store in list
  monthly_rasters[[date_str]] <- monthly_mean
}

# Stack all monthly rasters
aod_stack <- rast(monthly_rasters)

r_stack<-aod_stack

print(r_stack)
plot(r_stack[[1]])  # First month's AOD



r2<-r_stack
# Reproject if CRS mismatch
if (!crs(r2) == crs(shape2)) {
  message("CRS mismatch - reprojecting")
  r2 <- project(r2, crs(shape2))
} else {
  message("CRS aligned")
}

plot(r2[[1]])



# Crop and mask the raster stack
r2 <- mask(crop(r2, shape2), shape2)

# Plot first few layers for sanity check
plot(r2[[1]])
plot(r2[[2]])

# Prepare final_data for merging
final_data <- final_data %>%
  mutate(
    year_month = paste(year, sprintf("%02d", month), sep = " "),
    code = paste(year_month, district_country, sep = " "),
    aod = NA  # initialize
  )

# Loop through each raster layer
for (i in 1:nlyr(r2)) {
  lyr <- r2[[i]]
  lyr_name <- names(r2)[i]
  cat("Processing:", lyr_name, "\n")
  
  # Extract year and month from layer name like "AOD_200301"
  ym_parts <- str_match(lyr_name, "(\\d{4})(\\d{2})")
  year <- as.integer(ym_parts[, 2])
  month <- as.integer(ym_parts[, 3])
  
  # Extract mean AOD for each district in shape2
  extracted <- data.frame(
    shape2,
    aod_temp = terra::extract(lyr, shape2, fun = mean, na.rm = TRUE)[, 2]
  ) %>%
    mutate(
      year = year,
      month = month,
      year_month = paste(year, sprintf("%02d", month), sep = " "),
      district_country = paste(NAME_2, COUNTRY, sep = " "),
      code = paste(year_month, district_country, sep = " ")
    )
  
  # Merge with final_data and fill missing aod
  final_data <- merge(final_data, extracted[, c("code", "aod_temp")], by = "code", all.x = TRUE)
  final_data$aod[is.na(final_data$aod)] <- final_data$aod_temp[is.na(final_data$aod)]
  final_data$aod_temp <- NULL
}

print(head(final_data))


setwd("~/INLA project")

#### LANDCOVER DATA####  

setwd("~/INLA project/landcovertest")


nc_files <- list.files(pattern = "\\.nc$", full.names = TRUE)

# Initialize an empty list to hold the monthly rasters
monthly_rasters <- list()

# Loop through each file

for (file in nc_files) {
  message("Processing: ", basename(file))
  
  r <- rast(file)
  var_names <- names(r)
  print(var_names)
}

plot(r[[1]])
varnames(r)


r2<-subset(r,1:1166)
r2 <- r2[[ (nlyr(r2) - 12) : nlyr(r2) ]]
r3<-subset(r,1167:2332)
r3 <- r3[[ (nlyr(r3) - 12) : nlyr(r3) ]]

r4<-subset(r,2333:3498)
r4 <- r4[[ (nlyr(r4) - 12) : nlyr(r4) ]]
r5<-subset(r,3499:4664)
r5 <- r5[[ (nlyr(r5) - 12) : nlyr(r5) ]]
r6<-subset(r,4665:5830)
r6 <- r6[[ (nlyr(r6) - 12) : nlyr(r6) ]]
r7<-subset(r,5831:6996)
r7 <- r7[[ (nlyr(r7) - 12) : nlyr(r7) ]]
r8<-subset(r,6997:8162)
r8 <- r8[[ (nlyr(r8) - 12) : nlyr(r8) ]]
r9<-subset(r,8163:9328)
r9 <- r9[[ (nlyr(r9) - 12) : nlyr(r9) ]]
r10<-subset(r,9329:10494)
r10 <- r10[[ (nlyr(r10) - 12) : nlyr(r10) ]]
r11<-subset(r,10495:11660)
r11 <- r11[[ (nlyr(r11) - 12) : nlyr(r11) ]]
r12<-subset(r,11661:13992)
r12 <- r12[[ (nlyr(r12) - 12) : nlyr(r12) ]]
r13<-subset(r,11662:15158)
r13 <- r13[[ (nlyr(r13) - 12) : nlyr(r13) ]]
r14<-subset(r,15158:16324)
r14 <- r14[[ (nlyr(r14) - 12) : nlyr(r14) ]]



align_crs <- function(raster_obj, vector_obj) {
  if (!crs(raster_obj) == crs(vector_obj)) {
    message("CRS mismatch - Reprojecting")
    raster_obj <- project(raster_obj, crs(vector_obj))
  } else {
    message("CRS aligned")
  }
  return(raster_obj)
}

align_crs(r2,shape2)
align_crs(r3,shape2)
align_crs(r4,shape2)
align_crs(r5,shape2)
align_crs(r6,shape2)
align_crs(r7,shape2)
align_crs(r8,shape2)
align_crs(r9,shape2)
align_crs(r10,shape2)
align_crs(r11,shape2)
align_crs(r12,shape2)
align_crs(r13,shape2)
align_crs(r14,shape2)



r2 <- mask(crop(r2, shape2), shape2)
r3 <- mask(crop(r3, shape2), shape2)
r4 <- mask(crop(r4, shape2), shape2)
r5 <- mask(crop(r5, shape2), shape2)
r6 <- mask(crop(r6, shape2), shape2)
r7 <- mask(crop(r7, shape2), shape2)
r8 <- mask(crop(r8, shape2), shape2)
r9 <- mask(crop(r9, shape2), shape2)
r10 <- mask(crop(r10, shape2), shape2)
r11 <- mask(crop(r11, shape2), shape2)
r12 <- mask(crop(r12, shape2), shape2)
r13 <- mask(crop(r13, shape2), shape2)
r14 <- mask(crop(r14, shape2), shape2)


final_data$code <- paste(final_data$year, final_data$district_country, sep = " ")


# Prepare final_data for merging
final_data <- final_data %>%
  mutate(
    code = paste(year, district_country, sep = " "),
    landcover = NA  # initialize
  )


# Assign the correct names (2003 to 2015) to the layers
years <- 2003:2015

# Function to rename layers for each raster stack
rename_raster_layers <- function(raster_stack, years) {
  # Make sure the number of layers matches the number of years
  if (nlyr(raster_stack) == length(years)) {
    names(raster_stack) <- as.character(years)
  } else {
    stop("The number of layers does not match the number of years.")
  }
  return(raster_stack)
}

# Now, apply this to all your raster stacks
r2 <- rename_raster_layers(r2, years)
r3 <- rename_raster_layers(r3, years)
r4 <- rename_raster_layers(r4, years)
r5 <- rename_raster_layers(r5, years)
r6 <- rename_raster_layers(r6, years)
r7 <- rename_raster_layers(r7, years)
r8 <- rename_raster_layers(r8, years)
r9 <- rename_raster_layers(r9, years)
r10 <- rename_raster_layers(r10, years)
r11 <- rename_raster_layers(r11, years)
r12 <- rename_raster_layers(r12, years)
r13 <- rename_raster_layers(r13, years)
r14 <- rename_raster_layers(r14, years)




#### LANDCOVER DATA####  

setwd("~/INLA project/landcovertest3")


nc_files <- list.files(pattern = "\\.nc$", full.names = TRUE)

for (file in nc_files) {
  message("Processing: ", basename(file))
  
  r <- rast(file)
  var_names <- names(r)
  print(var_names)
}

plot(r[[1]])
varnames(r)




r2a<-subset(r,1:86)
r2a <- r2a[[2:8]]
r3a<-subset(r,87:172)
r3a <- r3a[[2:8]]

r4a<-subset(r,173:258)
r4a <- r4a[[2:8]]
r5a<-subset(r,259:344)
r5a <- r5a[[2:8]]
r6a<-subset(r,345:430)
r6a <- r6a[[2:8]]
r7a<-subset(r,431:516)
r7a <- r7a[[2:8]]
r8a<-subset(r,517:602)
r8a <- r8a[[2:8]]
r9a<-subset(r,603:774)
r9a <- r9a[[2:8]]
r10a<-subset(r,775:860)
r10a <- r10a[[2:8]]
r11a<-subset(r,861:946)
r11a <- r11a[[2:8]]
r12a<-subset(r,947:1033)
r12a <- r12a[[2:8]]
r13a<-subset(r,1033:1120)
r13a <- r13a[[2:8]]
r14a<-subset(r,1119:1205)
r14a <- r14a[[2:8]]



align_crs <- function(raster_obj, vector_obj) {
  if (!crs(raster_obj) == crs(vector_obj)) {
    message("CRS mismatch - Reprojecting")
    raster_obj <- project(raster_obj, crs(vector_obj))
  } else {
    message("CRS aligned")
  }
  return(raster_obj)
}

align_crs(r2a,shape2)
align_crs(r3a,shape2)
align_crs(r4a,shape2)
align_crs(r5a,shape2)
align_crs(r6a,shape2)
align_crs(r7a,shape2)
align_crs(r8a,shape2)
align_crs(r9a,shape2)
align_crs(r10a,shape2)
align_crs(r11a,shape2)
align_crs(r12a,shape2)
align_crs(r13a,shape2)
align_crs(r14a,shape2)



r2a <- mask(crop(r2a, shape2), shape2)
r3a <- mask(crop(r3a, shape2), shape2)
r4a <- mask(crop(r4a, shape2), shape2)
r5a <- mask(crop(r5a, shape2), shape2)
r6a <- mask(crop(r6a, shape2), shape2)
r7a <- mask(crop(r7a, shape2), shape2)
r8a <- mask(crop(r8a, shape2), shape2)
r9a <- mask(crop(r9a, shape2), shape2)
r10a <- mask(crop(r10a, shape2), shape2)
r11a <- mask(crop(r11a, shape2), shape2)
r12a <- mask(crop(r12a, shape2), shape2)
r13a <- mask(crop(r13a, shape2), shape2)
r14a <- mask(crop(r14a, shape2), shape2)


years <- 2016:2022

r2a <- rename_raster_layers(r2a, years)
r3a <- rename_raster_layers(r3a, years)
r4a <- rename_raster_layers(r4a, years)
r5a <- rename_raster_layers(r5a, years)
r6a <- rename_raster_layers(r6a, years)
r7a <- rename_raster_layers(r7a, years)
r8a <- rename_raster_layers(r8a, years)
r9a <- rename_raster_layers(r9a, years)
r10a <- rename_raster_layers(r10a, years)
r11a <- rename_raster_layers(r11a, years)
r12a <- rename_raster_layers(r12a, years)
r13a <- rename_raster_layers(r13a, years)
r14a <- rename_raster_layers(r14a, years)



r2 <- c(r2, r2a)
r3 <- c(r3, r3a)
r4 <- c(r4, r4a)
r5 <- c(r5, r5a)
r6 <- c(r6, r6a)
r7 <- c(r7, r7a)
r8 <- c(r8, r8a)
r9 <- c(r9, r9a)
r10 <- c(r10, r10a)
r11 <- c(r11, r11a)
r12 <- c(r12, r12a)
r13 <- c(r13, r13a)
r14 <- c(r14, r14a)


# Year-mapping aware processing function
process_raster_stack <- function(raster_stack, shape2, final_data, value_column = "landcover") {
  
  for (i in 1:nlyr(raster_stack)) {
    lyr <- raster_stack[[i]]
    lyr_name <- names(raster_stack)[i]
    cat("Processing:", lyr_name, "\n")
    
    # Directly use the layer name as the year
    year <- as.integer(lyr_name)
    
    # Extract mean values per district
    extracted <- data.frame(
      shape2,
      value_temp = terra::extract(lyr, shape2, fun = mean, na.rm = TRUE)[, 2]
    ) %>%
      mutate(
        year = year,
        district_country = paste(NAME_2, COUNTRY, sep = " "),
        code = paste(year, district_country, sep = " ")
      ) %>%
      dplyr::select(code, value_temp)
    
    # Merge and fill NA values only
    final_data <- merge(final_data, extracted, by = "code", all.x = TRUE)
    final_data[[value_column]][is.na(final_data[[value_column]])] <- final_data$value_temp[is.na(final_data[[value_column]])]
    
    # Remove the temporary column
    final_data$value_temp <- NULL
  }
  
  return(final_data)
}


final_data$landcover <- NULL
final_data$landcover <- NA
final_data <- process_raster_stack(r2, shape2, final_data, value_column = "landcover")
final_data$"primf"<-final_data$landcover
final_data$landcover <- NULL
final_data$landcover <- NA
final_data <- process_raster_stack(r3, shape2, final_data, value_column = "landcover")
final_data$"primn"<-final_data$"landcover"
final_data$landcover <- NULL
final_data$landcover <- NA
final_data <- process_raster_stack(r4, shape2, final_data, value_column = "landcover")
final_data$"secdf"<-final_data$"landcover"
final_data$landcover <- NULL
final_data$landcover <- NA
final_data <- process_raster_stack(r5, shape2, final_data, value_column = "landcover")
final_data$"secdn"<-final_data$"landcover"
final_data$landcover <- NULL
final_data$landcover <- NA
final_data <- process_raster_stack(r6, shape2, final_data, value_column = "landcover")
final_data$"urban"<-final_data$"landcover"
final_data$landcover <- NULL
final_data$landcover <- NA
final_data <- process_raster_stack(r7, shape2, final_data, value_column = "landcover")
final_data$"c3ann"<-final_data$"landcover"
final_data$landcover <- NULL
final_data$landcover <- NA
final_data <- process_raster_stack(r8, shape2, final_data, value_column = "landcover")
final_data$"c4ann"<-final_data$"landcover"
final_data$landcover <- NULL
final_data$landcover <- NA
final_data <- process_raster_stack(r9, shape2, final_data, value_column = "landcover")
final_data$"c3per"<-final_data$"landcover"
final_data$landcover <- NULL
final_data$landcover <- NA
final_data <- process_raster_stack(r10, shape2, final_data, value_column = "landcover")
final_data$"c4per"<-final_data$"landcover"
final_data$landcover <- NULL
final_data$landcover <- NA
final_data <- process_raster_stack(r11, shape2, final_data, value_column = "landcover")
final_data$"c3nfx"<-final_data$"landcover"
final_data$landcover <- NULL
final_data$landcover <- NA
final_data <- process_raster_stack(r12, shape2, final_data, value_column = "landcover")
final_data$"range"<-final_data$"landcover"
final_data$landcover <- NULL
final_data$landcover <- NA
final_data <- process_raster_stack(r13, shape2, final_data, value_column = "landcover")
final_data$"secmb"<-final_data$"landcover"
final_data$landcover <- NULL
final_data$landcover <- NA
final_data <- process_raster_stack(r14, shape2, final_data, value_column = "landcover")
final_data$"secma"<-final_data$"landcover"
final_data$landcover <- NULL


write_xlsx(final_data, "testing.xlsx")

