library(sf)
# Read in shapefile of Africa
# Set working directory
setwd("C:/Users/mvc32/OneDrive - University of Cambridge/Documents/Climate_meningitis_belt")
shape2 <- st_read("Shapefile_improved.shp")

####MESH CREATION ====

#MESH allows INLA to be computationally efficentnas it uses SPDE to estimate spatial autocorrelation of data. 
#This involves ussing a mesh of discrete sampling location, interpolated to
#estimate continuous rocess in space



# Make geometries valid
shape2 <- st_make_valid(shape2)

# Recompute centroids and extract coordinates, rather than using full geometry of shapefile
Locations <- shape2 |> 
  st_centroid() |> 
  st_coordinates()

# Convert to data frame for clarity
Locations_df <- data.frame(Longitude = Locations[, 1], Latitude = Locations[, 2])

Locations_jittered <- jitter(Locations)

# Create the mesh with appropriate max.edge settings
#MESH A juittering to show spacing
MeshA <- inla.mesh.2d(loc = Locations_jittered, max.edge = c(20, 40))

# View the first rows of the result
head(Locations_df)
# Convert to a data frame for clarity if needed
#MESH B exploratory analysis
MeshB <- inla.mesh.2d(Locations, max.edge = c(20, 40))
# MESH C use within paper, but will likely need some adjusting/playing with
MeshC <- inla.mesh.2d(Locations, max.edge = c(10, 20))

Mesh <- MeshB

plot(MeshA)

plot(MeshB)

plot(MeshC)

points(Locations, col = "red", pch = 2)

#Triangle size in MESH is v important
#determined using combination of max.edge and cut off
#using smaller trinagles, more precision but more computational power


library(sf)
library(spdep)

####NEIGHBOURHOOD MATRICES ====

#This function returns a neigh-
# bors list nb based on counties with contiguous boundaries. Each element of
#the list nb represents one county and contains the indices of its neighbors. For
#example, nb[[2]] contains the neighbors of county 2.
# Convert to Spatial object
shapefile_spatial <- as_Spatial(shape2)

# Create neighborhood matrix
nb <- poly2nb(shapefile_spatial)

# Plot neighborhood structure
# Extract centroids directly from the `sf` object
centroids <- st_centroid(shape2)

# Plot neighborhood structure
plot(st_geometry(shape2), border = "gray")  # Plot the shapefile
plot(nb, coordinates(as_Spatial(centroids)), add = TRUE, col = "red")  # Add neighborhood structure
