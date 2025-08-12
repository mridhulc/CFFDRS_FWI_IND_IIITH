library(terra)
library(ncdf4)
library(cffdrs)

# =========================
# 1. File paths
# =========================
meteo_file <- "D:/mridhu/exp/data_set/2023/weather.nc"
prec_file  <- "D:/mridhu/exp/data_set/2023/prec.nc"

# --- Load specific variables from meteo.nc ---
temp <- rast(paste0("NETCDF:", meteo_file, ":t2m"))   # Kelvin
d2m  <- rast(paste0("NETCDF:", meteo_file, ":d2m"))   # Kelvin
u10  <- rast(paste0("NETCDF:", meteo_file, ":u10"))   # m/s
v10  <- rast(paste0("NETCDF:", meteo_file, ":v10"))   # m/s

# --- Load precipitation from separate file ---
prec <- rast(paste0("NETCDF:", prec_file, ":tp"))     # meters

# =========================
# 2. Extract time dimension
# =========================
nc <- nc_open(meteo_file)
time_vals <- ncvar_get(nc, "valid_time")
time_units <- ncatt_get(nc, "valid_time", "units")$value
nc_close(nc)

time_origin <- as.POSIXct(strsplit(time_units, "since ")[[1]][2], tz = "UTC")
date_seq <- time_origin + time_vals * 3600
print(range(date_seq))

# =========================
# 3. Unit conversions
# =========================
temp_c   <- temp - 273.15
d2m_c    <- d2m  - 273.15
prec_mm  <- prec * 1000
ws_ms    <- sqrt(u10^2 + v10^2)
ws_kmh   <- ws_ms * 3.6

# =========================
# 4. Relative Humidity
# =========================
rh_calc <- function(t_c, td_c) {
  exp_dew  <- exp((17.625 * td_c) / (243.04 + td_c))
  exp_temp <- exp((17.625 * t_c)  / (243.04 + t_c))
  rh <- 100 * exp_dew / exp_temp
  rh[rh > 100] <- 100
  rh[rh < 0]   <- 0
  return(rh)
}
rh <- rh_calc(temp_c, d2m_c)

# =========================
# 5. Align rasters
# =========================
rh        <- resample(rh, temp_c, method="bilinear")
ws_kmh    <- resample(ws_kmh, temp_c, method="bilinear")
prec_mm   <- resample(prec_mm, temp_c, method="bilinear")

# =========================
# 6. FWI calculation
# =========================
nt <- nlyr(temp_c)
fwi_list <- vector("list", nt)

for (i in 1:nt) {
  input_stack <- c(
    temp_c[[i]],
    rh[[i]],
    ws_kmh[[i]],
    prec_mm[[i]]
  )
  names(input_stack) <- c("temp", "rh", "ws", "prec")
  fwi_list[[i]] <- fwiRaster(input_stack)
}

fwi_result <- do.call(c, fwi_list)

# =========================
# 7. Name layers by component and date
# =========================
var_types <- unique(names(fwi_result))  # e.g., FFMC, DMC, ...
names(fwi_result) <- as.vector(sapply(var_types, function(v) {
  paste0(v, "_", format(date_seq, "%Y%m%d"))
}))

# =========================
# 8. Save to NetCDF with proper variable names
# =========================
out_file <- "D:/mridhu/exp/fwi_output_2023.nc"
if (file.exists(out_file)) file.remove(out_file)

# Create a list of rasters (one per variable)
comp_list <- list()
for (comp in var_types) {
  comp_layers <- fwi_result[[grep(paste0("^", comp, "_"), names(fwi_result))]]
  names(comp_layers) <- format(date_seq, "%Y%m%d")  # only dates inside layer names
  comp_list[[comp]] <- comp_layers
}

# Convert to SpatRasterDataset (SDS)
fwi_sds <- sds(comp_list)

# Write all variables in one go
writeCDF(fwi_sds, filename = out_file, overwrite = TRUE)

# Check
nc <- nc_open(out_file)
print(names(nc$var))
nc_close(nc)


# =========================
# 10. Plot first timestep of FWI
# =========================
custom_palette <- colorRampPalette(c("blue", "green", "yellow", "red"))
plot(fwi_result[[which(grepl("^FWI_", names(fwi_result)))[1]]],
     main = paste0("FWI - ", format(date_seq[1], "%Y-%m-%d")),
     col = custom_palette(100))
