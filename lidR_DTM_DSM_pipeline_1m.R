#!/usr/bin/env Rscript
# LiDAR processing pipeline for LOCAL MACHINE (RStudio)
# STRICT ADHERENCE to lidR for DTM/DSM and terra for Slope/Aspect
#
# ============================================================================
# LOCAL CONFIGURATION - MODIFY THESE FOR YOUR DATA
# ============================================================================
# Campaign identification and CRS settings
campaign_name <- "ARRA_CA_GOLDENGATE_2010"  # Used for logging/identification
epsg_code <- "26910"
conversion_needed <- FALSE
conversion_factor <- 1.0

# Local directory paths - MODIFY THESE TO YOUR LOCAL PATHS
input_dir <- "/Volumes/WesternDigital/ARRA_CA_GOLDENGATE_2010"  # Directory containing .laz files
output_dir <- "/Volumes/WesternDigital/ARRA_CA_GOLDENGATE_2010_dem"  # Where to save outputs

# Processing settings
num_cores <- 1  # Adjust based on your local machine (use detectCores() - 1 for auto)
# ============================================================================

cat("Processing campaign:", campaign_name, "\n")

# Load required libraries
suppressPackageStartupMessages({
  library(terra)
  library(lidR)
  library(doParallel)
  library(foreach)
  library(parallel)
})

# ============================================================================
# SETUP OUTPUT DIRECTORY
# ============================================================================
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Setup log file
log_file <- file.path(output_dir, paste0(campaign_name, "_processing_log.txt"))

# Initialize log file
writeLines(paste("LiDAR Processing Log for", campaign_name, "-", Sys.time()), log_file)

cat("\nConfiguration:\n")
cat("  Number of cores:", num_cores, "\n")
cat("\nCampaign-specific parameters:\n")
cat("  Campaign name:", campaign_name, "\n")
cat("  EPSG code:", epsg_code, "\n")
cat("  Conversion needed:", conversion_needed, "\n")
cat("  Conversion factor:", conversion_factor, "\n")
cat("\nPaths:\n")
cat("  Input directory:", input_dir, "\n")
cat("  Output directory:", output_dir, "\n")
cat("  Log file:", log_file, "\n\n")

# Find all .laz files
laz_files <- list.files(input_dir, pattern="\\.(las|laz)$", full.names=TRUE)
cat("Found", length(laz_files), "LiDAR files to process\n")

if (length(laz_files) == 0) {
  stop("No .las or .laz files found in the input directory!")
}

# ============================================================================
# PER-FILE PROCESSING FUNCTION
# ============================================================================
process_file <- function(laz_file) {
  base_name <- tools::file_path_sans_ext(basename(laz_file))
  
  # Output file paths
  dtm_file   <- file.path(output_dir, paste0(base_name, "_dtm.tif"))
  dsm_file   <- file.path(output_dir, paste0(base_name, "_dsm.tif"))
  chm_file   <- file.path(output_dir, paste0(base_name, "_chm.tif"))
  slope_file <- file.path(output_dir, paste0(base_name, "_slope.tif"))
  aspect_file <- file.path(output_dir, paste0(base_name, "_aspect.tif"))
  
  # Check if already processed (using CHM as indicator)
  if (file.exists(chm_file)) {
    return(list(
      status = "skipped",
      message = paste("✔", base_name, "– already processed"),
      error = NULL
    ))
  }
  
  # Wrap processing in tryCatch for error handling
  result <- tryCatch({
    cat("Processing:", base_name, "...\n")
    
    # Step 1: Extract bounding box information from the .laz file
    header <- lidR::readLASheader(laz_file)
    bbox <- header@PHB
    res <- 1.0
    
    xrange <- bbox$`Max X` - bbox$`Min X`
    yrange <- bbox$`Max Y` - bbox$`Min Y`
    target_ext <- terra::ext(bbox$`Min X`, bbox$`Max X`, bbox$`Min Y`, bbox$`Max Y`)
    
    # Step 2: Read LAS file with lidR
    las <- lidR::readLAS(laz_file)
    
    # Set CRS for the LAS object
    #lidR::epsg(las) <- as.numeric(epsg_code) old-version didn't work
    # New version
    lasR::set_crs(as.numeric(epsg_code))
    
    # Apply unit conversion to point cloud if needed
    if (conversion_needed) {
      cat("  Applying unit conversion to point cloud (factor:", conversion_factor, ")\n")
      las@data$Z <- las@data$Z * conversion_factor
    }
    
    # Step 3a.i: Create DTM from Class 2 (ground) points - THIS WILL BE OUR TEMPLATE
    cat("  Creating DTM (this will serve as template for all other products)...\n")
    # Filter for Class 2 points only
    ground_las <- lidR::filter_poi(las, Classification == 2)
    
    # Check if ground points exist
    if (nrow(ground_las@data) == 0) {
      stop("No ground points found. Impossible to compute a DTM.")
    }
    
    # Generate DTM using triangulation (TIN)
    dtm_r <- lidR::rasterize_terrain(ground_las, res = res, algorithm = tin(), extent = target_ext)
    
    # Critical check: If DTM generation fails, skip this file
    if (is.null(dtm_r) || all(is.na(terra::values(dtm_r)))) {
      stop("Failed to generate DTM from Class 2 points")
    }
    
    # Ensure CRS is set for DTM (our template)
    terra::crs(dtm_r) <- paste0("EPSG:", epsg_code)
    
    # Step 3a.ii: Create water mask from Class 9 (water) points and conform to DTM
    water_mask_r <- NULL
    water_las <- lidR::filter_poi(las, Classification == 9)
    
    # Check if water points exist
    if (nrow(water_las@data) > 0) {
      cat("  Found", nrow(water_las@data), "water points\n")
      # Generate water mask using lasR point-to-raster method
      water_mask_raw <- tryCatch({
        # Extract coordinates and create temporary template raster
        coords <- data.frame(X = water_las@data$X, Y = water_las@data$Y, Z = water_las@data$Z)
        
        # Create raster template matching target extent and resolution
        water_template <- terra::rast(target_ext, resolution = res, crs = paste0("EPSG:", epsg_code))
        
        # Rasterize water points using terra
        water_rast <- terra::rasterize(coords[,1:2], water_template, field = coords$Z, fun = "mean")
        water_rast
      }, error = function(e) {
        cat("  Warning: Could not create water mask -", e$message, "\n")
        NULL  # No valid water surface could be created
      })
      
      # If water mask was created, resample to match DTM exactly
      if (!is.null(water_mask_raw)) {
        cat("  Conforming water mask to DTM template...\n")
        water_mask_r <- terra::resample(water_mask_raw, dtm_r, method = "near")
        rm(water_mask_raw)
      }
    } else {
      cat("  No water points found\n")
    }
    
    # Step 3b: Create DSM and conform to DTM template
    cat("  Creating DSM and conforming to DTM template...\n")
    dsm_raw <- lidR::rasterize_canopy(las, res = res, algorithm = dsmtin(), extent = target_ext)
    
    # Resample DSM to match DTM exactly
    dsm_r <- terra::resample(dsm_raw, dtm_r, method = "bilinear")
    rm(dsm_raw)
    
    # VERIFICATION: All rasters now have identical extents to DTM
    cat("  Verifying DSM matches DTM template...\n")
    if (!terra::compareGeom(dtm_r, dsm_r, lyrs = FALSE, crs = FALSE, warncrs = FALSE, ext = TRUE, rowcol = TRUE, res = TRUE)) {
      stop("DSM doesn't match DTM template after resampling - this should not happen")
    }
    
    # Helper function for true nearest neighbor interpolation
    nearest_neighbor_fill <- function(raster_data) {
      # Get coordinates and values
      coords <- terra::xyFromCell(raster_data, 1:terra::ncell(raster_data))
      values <- terra::values(raster_data)[,1]
      
      # Find NA and valid cells
      na_indices <- which(is.na(values))
      valid_indices <- which(!is.na(values))
      
      if (length(na_indices) > 0 && length(valid_indices) > 0) {
        cat("    Filling", length(na_indices), "NA cells using", length(valid_indices), "valid neighbors...\n")
        
        # Get coordinates for NA and valid cells
        na_coords <- coords[na_indices, , drop = FALSE]
        valid_coords <- coords[valid_indices, , drop = FALSE]
        
        # For each NA cell, find its nearest neighbor
        for (i in seq_len(nrow(na_coords))) {
          # Calculate distances to all valid cells using vectorized operations
          distances <- sqrt((valid_coords[,1] - na_coords[i,1])^2 + (valid_coords[,2] - na_coords[i,2])^2)
          
          # Find nearest neighbor and assign its value
          nearest_idx <- valid_indices[which.min(distances)]
          values[na_indices[i]] <- values[nearest_idx]
        }
      }
      
      # Set values back to raster
      terra::values(raster_data) <- values
      return(raster_data)
    }
    
    # Step 3b.i: Fill data voids in DTM using true nearest neighbor interpolation
    cat("  Filling DTM voids with true nearest neighbor interpolation...\n")
    if (any(is.na(terra::values(dtm_r)))) {
      dtm_r <- nearest_neighbor_fill(dtm_r)
    }
    
    # Step 3b.ii: Fill data voids in DSM using true nearest neighbor interpolation
    cat("  Filling DSM voids with true nearest neighbor interpolation...\n")
    if (any(is.na(terra::values(dsm_r)))) {
      dsm_r <- nearest_neighbor_fill(dsm_r)
    }
    
    # Step 3c: Compute CHM and clamp negative values to 0
    # This now works perfectly since DTM and DSM have identical extents
    chm_r <- dsm_r - dtm_r
    chm_r[chm_r < 0] <- 0
    
    # Step 3d: Calculate Slope using DTM template
    cat("  Calculating slope from DTM...\n")
    slope_r <- terra::terrain(dtm_r, v = "slope", unit = "degrees")
    
    # Step 3e: Calculate Aspect using DTM template
    cat("  Calculating aspect from DTM...\n")
    aspect_r <- terra::terrain(dtm_r, v = "aspect", unit = "degrees")
    
    # Step 3f: Apply water mask to all products (all now conform to DTM template)
    if (!is.null(water_mask_r)) {
      cat("  Applying water mask to all products...\n")
      
      # Create binary water mask (1 where water exists, NA elsewhere)
      water_mask_binary <- water_mask_r
      water_mask_binary[!is.na(water_mask_binary)] <- 1
      
      # Apply mask: set water pixels to -9999
      dtm_r[!is.na(water_mask_binary)] <- -9999
      dsm_r[!is.na(water_mask_binary)] <- -9999
      chm_r[!is.na(water_mask_binary)] <- -9999
      slope_r[!is.na(water_mask_binary)] <- -9999
      aspect_r[!is.na(water_mask_binary)] <- -9999
    }
    
    # Ensure CRS is set for all derived rasters (they inherit DTM's grid)
    terra::crs(dsm_r) <- paste0("EPSG:", epsg_code)
    terra::crs(chm_r) <- paste0("EPSG:", epsg_code)
    terra::crs(slope_r) <- paste0("EPSG:", epsg_code)
    terra::crs(aspect_r) <- paste0("EPSG:", epsg_code)
    
    # Write output rasters
    cat("  Writing output files...\n")
    terra::writeRaster(dtm_r, dtm_file, overwrite = TRUE, NAflag = -9999, 
                       gdal = c("COMPRESS=DEFLATE"))
    terra::writeRaster(dsm_r, dsm_file, overwrite = TRUE, NAflag = -9999, 
                       gdal = c("COMPRESS=DEFLATE"))
    terra::writeRaster(chm_r, chm_file, overwrite = TRUE, NAflag = -9999, 
                       gdal = c("COMPRESS=DEFLATE"))
    terra::writeRaster(slope_r, slope_file, overwrite = TRUE, NAflag = -9999, 
                       gdal = c("COMPRESS=DEFLATE"))
    terra::writeRaster(aspect_r, aspect_file, overwrite = TRUE, NAflag = -9999, 
                       gdal = c("COMPRESS=DEFLATE"))
    
    # Write water mask if it exists
    if (!is.null(water_mask_r)) {
      terra::writeRaster(water_mask_r, water_mask_file, overwrite = TRUE, NAflag = -9999, 
                         gdal = c("COMPRESS=DEFLATE"))
    }
    
    # Cleanup
    rm(las, ground_las, dtm_r, dsm_r, chm_r, slope_r, aspect_r)
    if (!is.null(water_las)) rm(water_las)
    if (!is.null(water_mask_r)) rm(water_mask_r, water_mask_binary)
    gc()
    
    return(list(
      status = "success",
      message = paste("✓", base_name, "- processed successfully"),
      error = NULL
    ))
    
  }, error = function(e) {
    return(list(
      status = "failed",
      message = paste("✗", base_name, "- FAILED"),
      error = as.character(e)
    ))
  })
  
  return(result)
}

# ============================================================================
# PROCESSING OPTIONS: PARALLEL OR SEQUENTIAL
# ============================================================================
# For testing or if parallel processing causes issues, set use_parallel to FALSE
use_parallel <- TRUE

if (use_parallel && num_cores > 1) {
  # ============================================================================
  # PARALLEL PROCESSING
  # ============================================================================
  cat("\nSetting up parallel backend with", num_cores, "cores\n")
  cl <- makeCluster(num_cores)
  
  # Export required libraries to worker nodes
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(terra)
      library(lidR)
    })
  })
  
  registerDoParallel(cl)
  
  # Export variables and functions to worker nodes
  clusterExport(cl,
                varlist = c("process_file", "output_dir", "epsg_code", 
                            "conversion_needed", "conversion_factor"),
                envir = environment())
  
  cat("Starting parallel processing of", length(laz_files), "files\n\n")
  
  results <- foreach(laz = laz_files, .combine = rbind) %dopar% {
    result <- process_file(laz)
    data.frame(
      file = basename(laz),
      status = result$status,
      message = result$message,
      error = ifelse(is.null(result$error), "", result$error),
      stringsAsFactors = FALSE
    )
  }
  
  # Cleanup parallel cluster
  stopCluster(cl)
  
} else {
  # ============================================================================
  # SEQUENTIAL PROCESSING
  # ============================================================================
  cat("\nProcessing files sequentially\n")
  cat("Starting processing of", length(laz_files), "files\n\n")
  
  results <- data.frame(
    file = character(),
    status = character(),
    message = character(),
    error = character(),
    stringsAsFactors = FALSE
  )
  
  for (laz in laz_files) {
    result <- process_file(laz)
    results <- rbind(results, data.frame(
      file = basename(laz),
      status = result$status,
      message = result$message,
      error = ifelse(is.null(result$error), "", result$error),
      stringsAsFactors = FALSE
    ))
  }
}

# ============================================================================
# CLEANUP AND REPORTING
# ============================================================================
closeAllConnections()
gc()

# Write detailed results to log
cat("\n", file = log_file, append = TRUE)
cat("Processing completed at:", as.character(Sys.time()), "\n", file = log_file, append = TRUE)
cat("Total files processed:", nrow(results), "\n", file = log_file, append = TRUE)
cat("Successful:", sum(results$status == "success"), "\n", file = log_file, append = TRUE)
cat("Failed:", sum(results$status == "failed"), "\n", file = log_file, append = TRUE)
cat("Skipped:", sum(results$status == "skipped"), "\n\n", file = log_file, append = TRUE)

# Log failed files with error messages
failed_files <- results[results$status == "failed", ]
if (nrow(failed_files) > 0) {
  cat("FAILED FILES:\n", file = log_file, append = TRUE)
  for (i in 1:nrow(failed_files)) {
    cat(sprintf("  %s: %s\n", failed_files$file[i], failed_files$error[i]), 
        file = log_file, append = TRUE)
  }
}

# Print summary to console
cat("\n========================================\n")
cat("Processing summary for campaign", campaign_name, ":\n")
cat("========================================\n")
for (i in 1:nrow(results)) {
  cat(results$message[i], "\n")
}
cat("\nTotal files:", nrow(results), "\n")
cat("Successful:", sum(results$status == "success"), "\n")
cat("Failed:", sum(results$status == "failed"), "\n")
cat("Skipped:", sum(results$status == "skipped"), "\n")
cat("\nAll outputs in:", output_dir, "\n")
cat("Log file:", log_file, "\n")

# Show failed files summary if any
if (nrow(failed_files) > 0) {
  cat("\n⚠️  WARNING: Some files failed to process. Check the log file for details.\n")
  cat("Failed files:\n")
  for (i in 1:nrow(failed_files)) {
    cat("  -", failed_files$file[i], "\n")
  }
}