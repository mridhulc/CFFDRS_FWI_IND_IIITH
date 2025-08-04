# Load libraries
library(terra)
library(cffdrs)

# 1. Load original rasters (update paths as needed)
temp <- rast("D:/mridhu/exp/temperature.nc")       # Temperature in Kelvin
prec <- rast("D:/mridhu/exp/precipitation.nc")     # Precipitation in meters
u10 <- rast("D:/mridhu/exp/wind.nc", lyrs = 1)     # u-component wind (m/s)
v10 <- rast("D:/mridhu/exp/wind.nc", lyrs = 2)     # v-component wind (m/s)
rh <- rast("D:/mridhu/exp/humidity.nc")            # Relative Humidity (% or fraction)

# 2. Convert units
temp_c <- temp - 273.15       # Kelvin to Celsius
prec_mm <- prec * 1000        # meters to millimeters

# 3. Compute wind speed (m/s) and convert to km/h
ws_ms <- sqrt(u10^2 + v10^2)
ws_kmh <- ws_ms * 3.6

# 4. Ensure relative humidity is in %
if (max(values(rh), na.rm=TRUE) <= 1) {
  rh <- rh * 100
}

# 5. Align all rasters to the reference raster (temp_c)
rh_aligned      <- resample(rh, temp_c, method="bilinear")
ws_kmh_aligned  <- resample(ws_kmh, temp_c, method="bilinear")
prec_mm_aligned <- resample(prec_mm, temp_c, method="bilinear")

# 6. Check number of layers for each raster
n_temp <- nlyr(temp_c)
n_rh   <- nlyr(rh_aligned)
n_ws   <- nlyr(ws_kmh_aligned)
n_prec <- nlyr(prec_mm_aligned)

cat("Layers - temp:", n_temp, " RH:", n_rh, " WS:", n_ws, " Prec:", n_prec, "\n")

# 7. Synchronize layer counts: replicate single-layer rasters to match temp layers
if (n_rh == 1 && n_temp > 1) {
  rh_aligned <- rast(lapply(1:n_temp, function(i) rh_aligned[[1]]))
}
if (n_ws == 1 && n_temp > 1) {
  ws_kmh_aligned <- rast(lapply(1:n_temp, function(i) ws_kmh_aligned[[1]]))
}
if (n_prec == 1 && n_temp > 1) {
  prec_mm_aligned <- rast(lapply(1:n_temp, function(i) prec_mm_aligned[[1]]))
}

# 8. Clamp RH values to [0, 100] to avoid invalid values
rh_aligned[rh_aligned < 0] <- 0
rh_aligned[rh_aligned > 100] <- 100

# 9. Verify layer counts match after replication and clamping
stopifnot(
  nlyr(temp_c) == nlyr(rh_aligned),
  nlyr(rh_aligned) == nlyr(ws_kmh_aligned),
  nlyr(ws_kmh_aligned) == nlyr(prec_mm_aligned)
)

# 10. Compute FWI for each timestep individually if multi-layer
nt <- nlyr(temp_c)

if (nt == 1) {
  # Single timestep
  input_stack <- c(temp_c, rh_aligned, ws_kmh_aligned, prec_mm_aligned)
  names(input_stack) <- c("temp", "rh", "ws", "prec")
  fwi_result <- fwiRaster(input_stack)
} else {
  # Multiple timesteps
  fwi_list <- vector("list", nt)
  
  for (i in 1:nt) {
    input_stack <- c(
      temp_c[[i]],
      rh_aligned[[i]],
      ws_kmh_aligned[[i]],
      prec_mm_aligned[[i]]
    )
    names(input_stack) <- c("temp", "rh", "ws", "prec")
    fwi_list[[i]] <- fwiRaster(input_stack)
  }
  
  # Combine all time steps into one multi-layer raster
  fwi_result <- do.call(c, fwi_list)
}

# 11. Ensure unique layer names to avoid plotting and subsetting issues
nm <- names(fwi_result)
if (anyDuplicated(nm) > 0) {
  # Generate unique names by appending indices to duplicates
  unique_types <- unique(nm)
  new_names <- unlist(lapply(unique_types, function(type_name) {
    inds <- which(nm == type_name)
    if (length(inds) > 1) {
      sprintf("%s_%02d", type_name, inds)
    } else {
      type_name
    }
  }))
  names(fwi_result) <- new_names
}

# 12. Save the FWI result as a NetCDF file
writeRaster(fwi_result,
            filename = "D:/mridhu/exp/fwi_output.nc",
            filetype = "netCDF",
            overwrite = TRUE)

# 13. Plot the first timestep Fire Weather Index (FWI)
# Find the first FWI layer name (handles renamed layers)
fwi_layers <- grep("^FWI", names(fwi_result), value = TRUE)
if (length(fwi_layers) > 0) {
  plot(fwi_result[[fwi_layers[1]]], main = "CFFDRS Fire Weather Index (FWI) - First Time Step")
} else {
  # Fallback to first layer
  plot(fwi_result[[1]], main = "FWI - First Time Step")
}

