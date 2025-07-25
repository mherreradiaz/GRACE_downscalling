library(dplyr)
library(terra)
library(glue)
library(fs)
library(rnaturalearth)

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

chl <- ne_countries(country='chile',returnclass = 'sv') |> 
  project('EPSG:4326')

# TWSA

grace_files <- dir_ls('data/raw/raster/GRACE')
vars <- c('lwe_thickness', 'scale_factor', 'GAD', 'land_mask')

grace_vars <- lapply(vars, function(v) {
  rast(grace_files, subds = v) |> rotate() |> project('EPSG:4326')
}) |> 
  setNames(vars)

lwe_r <- (grace_vars$lwe_thickness*grace_vars$scale_factor*grace_vars$land_mask)-grace_vars$GAD

fechas <- substr(time(lwe_r),1,7)
fechas[c(112,145)] <- c('2011-12','2015-05')

grace_chl <- lwe_r |> 
  crop(chl)

grace_r <- grace_chl |> 
  setNames(fechas)

fechas_completas <- substr(seq(as.Date(paste0(min(fechas),'-01')),as.Date(paste0(max(fechas),'-01')),by='1 month'),1,7)
fechas_faltantes <- fechas_completas[!(fechas_completas %in% fechas)]

grace_faltantes <- grace_r[[seq_along(fechas_faltantes)]] |> 
  setValues(NA) |> 
  setNames(fechas_faltantes)

grace_completo <- lapply(fechas_completas, \(fecha) c(grace_r,grace_faltantes)[[fecha]]) |> rast()

plot(grace_r[[1:12]])

exportRast(grace_r,'data/processed/raster/GRACE_monthly/','GRACE_')

# GWSA

grace_twsa <- dir_ls('data/processed/raster/GRACE_monthly/TWSA/')
fechas_grace <- str_extract(grace_twsa, "\\d{4}-\\d{2}")

gldas_vars <- c('GWS','SWE','SoilMoist_P','CanopInt')
files <- dir_ls(glue('data/processed/raster/GLDAS_anomaly_monthly/{gldas_vars[1]}/'))
fechas_gldas <- str_extract(files, "\\d{4}-\\d{2}")

fechas <- fechas_grace[fechas_grace %in% fechas_gldas]

gldas_vars_r <- lapply(gldas_vars,\(var) {
  r <- dir_ls(glue('data/processed/raster/GLDAS_anomaly_monthly/{var}/')) |> 
    rast() |> 
    subset(fechas) |> 
    project(rast(grace_twsa[1]))
  r/10
}) |> setNames(gldas_vars)

grace_twsa_r <- rast(grace_twsa) |> 
  subset(fechas)
gldas_swe_r <- gldas_vars_r$SWE
gldas_sm_r <- gldas_vars_r$SoilMoist_P
gldas_can_r <- gldas_vars_r$CanopInt

grace_gwsa_r <- grace_twsa_r -(gldas_swe_r + gldas_sm_r + gldas_can_r) 

exportRast(grace_gwsa_r,'data/processed/raster/GRACE_monthly/GWSA/','GRACE_GWSA_')


