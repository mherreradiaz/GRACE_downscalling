library(tidyverse)
library(readr)
library(purrr)
library(fs)
library(glue)
options(timeout = 3000)

download_list <- read_csv('C:/Users/mauricio.herrera/Documents/Projectos/GRACE_downscalling/data/raw/txt/CHELSA_download.txt') |> 
  pull(1)

vars <- str_extract(download_list,'(?<=monthly/)[^/]+(?=/CHELSA)') |>  unique()

lapply(c('rsds'),\(var) {
  urls <- unique(grep(var,download_list,value=T))
  urls <- urls[which(str_extract(urls,'\\d{4}')< 1998)]
  dir.out <- glue('data/raw/{var}/')
  dir_create(dir.out)
  walk(urls, \(url) {
    out_file <- path(dir.out, path_file(url))
    download.file(url, out_file, mode = 'wb', quiet = F)
  })
})


library(terra)
library(rnaturalearth)

chl <- ne_countries(country='chile', returnclass = 'sv') |> 
  project('EPSG:4326')

long_map <- c(
  clt = "Total cloud cover",
  cmi = "Climate Moisture Index",
  hurs = 'Nearsurface relative humidity',
  pet = "Potential evapotranspiration (Penman-Monteith)",
  pr  = "Precipitation amount",
  rsds = 'Surface downwelling shortwave flux in air',
  tas = "Mean daily air temperature"
)
unit_map <- c(
  clt = "percent",
  cmi = "kg/m2/month",
  hurs = 'percent',
  pet = "mm/month",
  pr  = "mm/month",
  rsds = 'MJ/m2/day',
  tas = "degC"
)

dir <- 'C:/Users/mauricio.herrera/Documents/Projectos/raster_dataset/procesamiento/CHELSA/data'

var_names <- list.files(glue('{dir}/raw/raster'))[-1]

lapply(var_names,\(var) {
  
  var_files <- list.files(glue('{dir}/raw/raster/{var}'),full.names=T)
  
  fechas <- str_extract(var_files, "\\d{2}_\\d{4}") |> 
    str_replace("(\\d{2})_(\\d{4})", "\\2-\\1-01") |> 
    as.Date()
  
  var_files_sort <- var_files[order(fechas)]
  fechas_sort <- sort(fechas)
  
  var_r <- rast(var_files_sort) |> 
    crop(chl)
  names(var_r) <- glue('CHELSA_{var}_{fechas_sort}')
  time(var_r) <- fechas_sort
  varnames(var_r) <- long_map[var]
  units(var_r) <- unit_map[var]
  
  # var_nc <- rast(grep(var,list.files(glue('{dir}/processed/'),full.names=T),value = T))
  
  # var_total <- c(var_r,var_nc)
  # fechas_total <- time(var_total)
  
  out_file <- glue('{dir}/processed/CHELSA_{var}_{format(min(fechas), "%Y%m")}_{format(max(fechas), "%Y%m")}.nc')
  # out_file <- glue('{dir}/processed/CHELSA_{var}_{format(min(fechas_total), "%Y%m")}_{format(max(fechas_total), "%Y%m")}.nc')
  
  writeCDF(var_r, out_file,
  # writeCDF(var_total, out_file,
           varname    = var,
           longname   = long_map[var],
           unit       = unit_map[var],
           zname      = "time",
           prec       = "float",
           compression= 4,
           missval    = -9999,
           overwrite=T)
  
  
})

