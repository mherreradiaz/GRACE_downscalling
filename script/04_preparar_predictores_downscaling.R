library(fs)
library(terra)
library(rnaturalearth)
library(ncdf4)
library(sf)
library(tidyverse)
library(glue)

chl <- ne_countries(country='chile',returnclass = 'sv') |> 
  project('EPSG:4326')

# Preparar indicadores TerraClimate

dir <- 'data/raw/raster/TerraClimate/'
var_names <- c('AET','PPT','Q','SOIL')

tc_preds <- lapply(var_names, \(var) {
  files <- dir_ls(dir,recurse = T,regexp = glue('{var}_'), type = 'file')
  fechas <- str_extract(files,glue('(?<={var}_).*?(?=\\.tif)'))
  files |>
    rast() |> 
    setNames(fechas) |>
    crop(chl)
}) |> 
  setNames(var_names)

# Preparar predictores indicadores de sequía

dir <- 'data/raw/raster/DroughtIndex/'
index_names <- c('SPI','SPEI','EDDI','SSI','SETI')

drought_pred <- lapply(index_names, \(index) {
  files <- dir_ls(dir,recurse = T,regexp = index, type = 'file')
  fechas <- str_extract(files,'\\d{4}-\\d{2}')
  files |>
    rast() |> 
    setNames(fechas) |>
    resample(tc_preds[[1]],method = 'cubicspline') |> 
    crop(chl)
}) |> 
  setNames(index_names)

# Preparar DEM

dem_pred <- rast('data/raw/raster/DEM/CHL_elv_msk.tif') |> 
  setNames('DEM') |> 
  resample(tc_preds[[1]],method = 'cubicspline')

# Preparar predictores indicadores de vegetacion

dir <- 'data/processed/raster/MOD13A3/FILTERED/'
index_names <- c('NDVI','EVI')
  
vi_pred <- lapply(index_names, \(index) {
  files <- dir_ls(dir,recurse = T,regexp = index, type = 'file')
  fechas <- str_extract(files,'\\d{4}-\\d{2}')
  files |>
    rast() |> 
    setNames(fechas) |>
    resample(tc_preds[[1]],method = 'cubicspline') |> 
    crop(chl)
}) |> 
  setNames(index_names)

# Guardar predictores

lapply(c('AET','PPT','Q','SOIL'), \(var) 
       writeRaster(tc_preds[[var]],glue('data/processed/raster/downscaling_preds/{var}.tif')))
lapply(c('SPI','SPEI','EDDI','SSI','SETI'), \(var) 
       writeRaster(drought_pred[[var]],glue('data/processed/raster/downscaling_preds/{var}.tif')))
writeRaster(dem_pred,'data/processed/raster/downscaling_preds/DEM.tif',
            overwrite = T)
lapply(c('NDVI','EVI'), \(var) 
       writeRaster(vi_pred[[var]],glue('data/processed/raster/downscaling_preds/{var}.tif')))
