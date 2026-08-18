# Load required libraries
library(sf)
library(ggplot2)
library(rgplates)
library(readxl)
library(dplyr)
library(geojsonsf)  

#### GPS Reconstruction ####
# Read diversity data from an Excel file
diversity <- read_excel("Database.xlsx", sheet = "Diversity")

# Function to reconstruct points for a given time slice
reconstruct_points <- function(time) {
    d <- diversity %>%
        filter(Time == time) %>%
        select(N, Longitude, Latitude) %>%
        as.matrix()
    
    rownames(d) <- d[, 1]  # Set the first column as row names
    d <- d[, -1]  # Remove the first column from the matrix
    
    points <- reconstruct(d, model = "PALEOMAP", age = time) %>%
        as.data.frame()
    return(points)
}

# Generate reconstructed points for multiple time slices
time_slices <- c(280, 300, 320, 340, 360)
points_list <- lapply(time_slices, reconstruct_points)
combined_points <- do.call(rbind, points_list)

# Save the reconstructed points to a CSV file
write.csv(combined_points, "reconstructed_points.csv", row.names = TRUE)

#### Paleogeographic Map Plotting Function ####
plot_paleomap <- function(time_ma) {
    # Reconstruct coastlines for the specified time slice
    coastlines <- reconstruct(x="coastlines", model = "PALEOMAP", age = time_ma)
    reconstructed_sf <- st_as_sf(coastlines)
    
    # Define Robinson projection
    robinson_proj <- "+proj=robin +lon_0=0 +datum=WGS84"
    reconstructed_robinson <- st_transform(reconstructed_sf, crs = robinson_proj)
    
    # Load Köppen-Geiger climate data
    data_koppen <- read.csv(paste0(time_ma, "Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv"))
    
    # Clean and transform climate data
    data_cleaned <- data_koppen[, c("lat", "lon", "koppen")]
    data_cleaned <- na.omit(data_cleaned)
    colnames(data_cleaned) <- c("latitude", "longitude", "koppen")
    data_sf <- st_as_sf(data_cleaned, coords = c("longitude", "latitude"), crs = 4326)
    
    # Load GPS data
    gps_data <- read.csv(paste0(time_ma, ".csv"))
    gps_sf <- st_as_sf(gps_data, coords = c("lon", "lat"), crs = 4326)
    
    # Extract coordinates for later use
    gps_coords <- st_coordinates(gps_sf)
    gps_sf$longitude <- gps_coords[, 1]
    gps_sf$latitude <- gps_coords[, 2]
    
    # Define color palette for Köppen climate zones
    cores_koppen <- c(
        "1" = "purple", "2" = "purple",   
        "3" = "deepskyblue", "4" = "deepskyblue",   "5" = "deepskyblue",
        "6" = "blue", "7" = "blue", "8" = "blue",
        "9" = "darkorange", "10" = "darkorange",
        "11" = "darkgreen", "12" = "darkgreen", "13" = "darkgreen"
    )
    
    # Generate the map plot
    ggplot() +
        geom_sf(data = reconstructed_robinson, fill = "black", color = "black") +
        geom_sf(data = data_sf, aes(color = as.factor(koppen)), size = 2, shape = 15, alpha = 0.08) +
        geom_sf(data = gps_sf, aes(color = "red"), shape = 16, size = 2, alpha = 1) +
        coord_sf(crs = robinson_proj) +
        theme_minimal() +
        theme(legend.position = "none") +
        scale_color_manual(values = cores_koppen, name = "Koppen Climate Zones") +
        ggtitle(paste0("Paleogeographic Map - ", time_ma, " Ma"))
}

# Generate and plot maps for all time slices
lapply(time_slices, plot_paleomap)
