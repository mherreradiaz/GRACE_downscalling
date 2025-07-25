library(tidyverse)
library(dplyr)
library(terra)
library(glue)
library(fs)
library(rnaturalearth)

# preprocesar imagenes diarias

chl <- ne_countries(country='chile', returnclass = 'sv') |> 
  project('EPSG:4326')

gldas_files <- dir_ls('data/raw/raster/GLDAS2.2/', regexp ='.nc')
fechas <- str_extract(gldas_files, '\\d{8}') |> 
  as.Date(format = '%Y%m%d')

gldas_vars <- c('GWS_tavg','SWE_tavg','SoilMoist_P_tavg','CanopInt_tavg')

for (i in seq_along(gldas_files)) {
  
  file <- gldas_files[i]
  r_full <- rast(file)
  
  for (var in gldas_vars) {

      r <- r_full[[var]] |> 
        project('EPSG:4326') |> 
        crop(chl) |> 
        setNames(fechas[i])
      
      var <- gsub('_tavg','',var)

      output_dir <- glue('data/processed/raster/GLDAS_daily/{var}/')
      dir_create(output_dir)

      writeRaster(r, 
                  filename = glue('{output_dir}{var}_{fechas[i]}.tif'), 
                  overwrite = TRUE)
      rm(r); gc()
    }
  rm(r_full); gc()
}

# agregacion mensual

exportRast <- \(r, output.dir, band.name, names = NULL, overwrite = F) {
  require(terra)
  require(glue)
  
  if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
  
  if (!grepl('/$', output.dir)) output.dir <- glue('{output.dir}/')
  
  if (!is.null(names)) {
    if (length(names) != nlyr(r)) stop("El largo de 'names' no coincide con el número de layers.")
    names(r) <- names
  }
  
  file_names <- glue('{output.dir}{band.name}{names(r)}.tif')
  
  lapply(1:nlyr(r), \(i) {
    ly <- r[[i]]
    writeRaster(ly, file_names[i], overwrite = overwrite)
  })
  
  return(invisible(NULL))
}

gldas_vars <- c('GWS','SWE','SoilMoist_P','CanopInt')

monthly_var_r <- list()

for (var in gldas_vars) {
  
  var_files <- dir_ls(glue('data/processed/raster/GLDAS_daily/{var}/'))
  fechas <- str_extract(var_files, "\\d{4}-\\d{2}-\\d{2}")
  
  daily_var_r <- rast(var_files)
  time(daily_var_r) <- as.Date(fechas)
  
  monthly_var_r <- tapp(daily_var_r, "yearmonths", mean, cores = 10, na.rm=T)
  
  names(monthly_var_r) <- substr(as.Date(paste0(gsub('ym_','',names(monthly_var_r)),'01'),format = '%Y%m%d'),1,7)
  
  exportRast(monthly_var_r,glue('data/processed/raster/GLDAS_monthly/{var}/'),glue('{var}_'))
  
}

# vars anomalies

gldas_vars <- c('GWS','SWE','SoilMoist_P','CanopInt')

for (var in gldas_vars) {
  
  var_files <- dir_ls(glue('data/processed/raster/GLDAS_monthly/{var}/'))
  fechas <- str_extract(var_files, "\\d{4}-\\d{2}")
  
  var_r <- rast(var_files)
  anomaly_var_r <- var_r-app(var_r[[which(names(var_r)<'2010-01')]],mean,na.rm=T)
  exportRast(anomaly_var_r,glue('data/processed/raster/GLDAS_anomaly_monthly/{var}/'),glue('{var}_anomaly_'))
  
}

