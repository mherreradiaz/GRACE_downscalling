library(tidyverse)
library(tidyterra)
library(terra)
library(parallel)
library(rnaturalearth)
library(glue)

terraOptions(threads = detectCores()-2)
extract_qa_bits <- \(x) {
  cbind(
    MODLAND_QA = bitwAnd(x, 3),
    VI_usefulness = bitwAnd(bitwShiftR(x, 2), 15),
    Aerosol_quantity = bitwAnd(bitwShiftR(x, 6), 3),
    Adjacent_cloud = bitwAnd(bitwShiftR(x, 8), 1),
    BRDF_corr = bitwAnd(bitwShiftR(x, 9), 1),
    Mixed_clouds = bitwAnd(bitwShiftR(x, 10), 1),
    Land_water_flag = bitwAnd(bitwShiftR(x, 11), 7),
    Snow_ice = bitwAnd(bitwShiftR(x, 14), 1),
    Shadow = bitwAnd(bitwShiftR(x, 15), 1)
  )
}
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
  
  # preview <- head(file_names, 10)
  # cat('Las capas se guardarán de la siguiente forma (10 primeras):\n')
  # cat(paste0(preview, collapse = '\n'), '\n')
  # 
  # respuesta <- readline(prompt = '¿Desea proceder con la exportación? [s/n]: ')
  # if (tolower(respuesta) != 's') {
  #   cat('Exportación cancelada.\n')
  #   return(invisible(NULL))
  # }
  
  lapply(1:nlyr(r), \(i) {
    ly <- r[[i]]
    writeRaster(ly, file_names[i], overwrite = overwrite)
  })
  
  cat('Exportación finalizada.\n')
}

chl <- ne_countries(country='chile',returnclass = 'sv') |> 
  project('EPSG:4326')

# preprocesar NDVI

qa_files <- list.files('data/raw/raster/MOD13A3.061', full.names=T, pattern='.tif') |>
  grep(pattern = 'VI_Q', value=T) |>
  grep(pattern = 'aux', value=T, invert=T)

fechas <- str_extract(qa_files, "(?<=doy)\\d{7}") |> 
  strptime(format='%Y%j') |> 
  as.Date()

r_qa_stack <- rast(qa_files) |> 
  crop(chl)

valores <- c(2112, 2116, 2181, 2372, 4160, 4164, 6208, 6212, 6277, 6468, 18500, 20548, 22596)
qa_mask <- classify(r_qa_stack, rcl=cbind(valores, 1), others=NA)

names(qa_mask) <- fechas

# qa_mask <- rast('data/raw/raster/MOD13A3.061/QA_mask/QA_mask.tif')

ndvi_files <- list.files('data/raw/raster/MOD13A3.061',full.names=T,pattern = '.tif') |> 
  grep(pattern = 'y_NDVI_d',value=T) |> 
  grep(pattern = 'aux',value = T, invert = T)

r_ndvi <- rast(ndvi_files) |> 
  crop(chl)

ndvi <- r_ndvi |>
  setNames(fechas)

ndvi_mask <- r_ndvi |>
  mask(qa_mask,maskvalues=F) |>
  setNames(fechas)

exportRast(ndvi,'data/processed/raster/MOD13A3/RAW/NDVI','monthly_NDVI_')
exportRast(ndvi_mask,'data/processed/raster/MOD13A3/FILTERED/NDVI','monthly_NDVI_')

evi_files <- list.files('data/raw/raster/MOD13A3.061',full.names=T,pattern = '.tif') |> 
  grep(pattern = 'y_EVI_d',value=T) |> 
  grep(pattern = 'aux',value = T, invert = T)

r_evi <- rast(evi_files) |> 
  crop(chl)

evi <- r_evi |>
  setNames(fechas)

evi_mask <- r_evi |>
  mask(qa_mask,maskvalues=F) |>
  setNames(fechas)

exportRast(evi,'data/processed/raster/MOD13A3/RAW/EVI','monthly_EVI_')
exportRast(evi_mask,'data/processed/raster/MOD13A3/FILTERED/EVI','monthly_EVI_')
