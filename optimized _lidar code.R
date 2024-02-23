

#pdf(NULL)
rm(list=ls())
############################
#library(rgdal)
#library(raster)
#library(tmap)
#library(tmaptools)
library(lidR)
#library(RStoolbox)
library(sf)
#library(sp)
#library(mapview)
# Load required libraries
#library(lidR)
library(future)
library(future.apply)


#Filtering the noise and eliminating it 
filter_noise = function(las, sensitivity)
{
  if (is(las, "LAS"))
  {
    p95 <- pixel_metrics(las, ~quantile(Z, probs = 0.95), 10)
    las <- merge_spatial(las, p95, "p95")
    las <- filter_poi(las, Z > 0,Z < p95*sensitivity)
    las$p95 <- NULL
    return(las)
  }
  
  if (is(las, "LAScatalog"))
  {
    res <- catalog_map(las, filter_noise, sensitivity = sensitivity)
    return(res)
  }
}





# Set the number of cores to use
num_cores <- 8

# Define the LAS catalog
ctg <- ctg<-readLAScatalog("S:\\Test")

# Set the output folder for normalized LAS files
output_folder <- "S:\\Test\\Filter"

# Define a function to normalize heights for a single LAS file
filter_noise_single <- function(las_file, output_folder) {
  # Normalize height using tin()
  las <- readLAS(las_file)
  las_norm <- filter_noise(las,sensitivity=1.2)
  
  # Construct the output file path
  output_file <- file.path(output_folder, paste0("filter_", basename(las_file)))
  
  # Write the normalized LAS to the output folder
  writeLAS(las_norm, output_file)
}

# Set up parallel processing with future
plan(multisession, workers = num_cores)

# Use future_lapply to normalize heights in parallel
future_lapply(ctg$filename, function(las_file) {
  filter_noise_single(las_file, output_folder)
})

# Stop parallel processing
plan(sequential)



############################# optional if the ground classification has not been done in lidar#####



# Set the number of cores to use
num_cores <- 8

# Define the LAS catalog
ctg <- readLAScatalog("S:\\Test\\Filter")

# Set the output folder for ground classified LAS files

output_folder <- "S:\\Tessara\\Andrew\\Classsified"

# Define a function to normalize heights for a single LAS file
classify_ground_height_single <- function(las_file, output_folder) {
  # Normalize height using tin()
  las <- readLAS(las_file)
  las_norm <- classify_ground(las, csf())
  
  # Construct the output file path
  output_file <- file.path(output_folder, paste0("classified_", basename(las_file)))
  
  # Write the normalized LAS to the output folder
  writeLAS(las_norm, output_file)
}

# Set up parallel processing with future
plan(multisession, workers = num_cores)

# Use future_lapply to normalize heights in parallel
future_lapply(ctg$filename, function(las_file) {
  classify_ground_height_single(las_file, output_folder)
})

# Stop parallel processing
plan(sequential)


################################################################### normalize the height ###########################


# Set the number of cores to use
num_cores <- 8

# Define the LAS catalog
ctg <- readLAScatalog("S:\\Test\\Filter")

# Set the output folder for normalized LAS files
output_folder <- "S:\\Test\\Normal"

# Define a function to normalize heights for a single LAS file
normalize_height_single <- function(las_file, output_folder) {
  # Normalize height using tin()
  las <- readLAS(las_file)
  las_norm <- normalize_height(las, tin())
  
  # Construct the output file path
  output_file <- file.path(output_folder, paste0("norm_", basename(las_file)))
  
  # Write the normalized LAS to the output folder
  writeLAS(las_norm, output_file)
}

# Set up parallel processing with future
plan(multisession, workers = num_cores)

# Use future_lapply to normalize heights in parallel
future_lapply(ctg$filename, function(las_file) {
  normalize_height_single(las_file, output_folder)
})

# Stop parallel processing
plan(sequential)




############################## Segment the trees and output as the shape file ###########

## Set the number of cores to use
num_cores <- 8

# Define the LAS catalog
ctg <- readLAScatalog("S:\\Test\\Normal")

# Set the output folder for normalized LAS files
output_folder <- "S:\\Test\\Chm"

# Define a function to normalize heights for a single LAS file
process_lidar_file <- function(las_file, output_folder) {
  # Read the LAS file
  las <- readLAS(las_file)
  
  # Perform processing without creating large objects in the global environment
  chm <- rasterize_canopy(las, res = 1, algorithm = p2r(subcircle = 0.2))
  ttops <- locate_trees(las, lmf(ws=15))
  
  # Skip the LAS file if no trees are detected
  if (length(ttops) == 0) {
    cat("No trees detected in:", las_file, "\n")
    return(NULL)
  }
  
  algo <- dalponte2016(chm, ttops)
  las <- segment_trees(las, algo)
  crowns <- crown_metrics(las, func = .stdmetrics, geom = "convex")
  crowns_df <- data.frame(crowns) # Convert to a data.frame
  
  # Extract the index from the LAS file name
  las_index <- sub(".las", "", basename(las_file))
  
  # Construct the output file path with the index
  output_file <- file.path(output_folder, paste0("Tree_", las_index, ".shp"))
  
  st_write(crowns, output_file, row.names=FALSE)
}

# Set up parallel processing with future
plan(multisession, workers = num_cores)

# Use future_lapply to process LAS files and save results in parallel
future_lapply(ctg$filename, function(las_file) {
  process_result <- process_lidar_file(las_file, output_folder)
  
  # Filter out NULL results (LAS files with no trees)
  if (!is.null(process_result)) {
    return(process_result)
  }
})

# Stop parallel processing
plan(sequential)
