library(tidyverse)
library(tidyterra)
library(terra)
library(parallel)
library(rnaturalearth)
library(glue)
library(fs)

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

dir <- 'data/raw/raster/CHELSA'

chelsa_files <- dir_ls(dir,recurse = T, glob = '*.tif') 
var_names <- list.files(dir)

lapply('tas',\(var) {
  
  var_files <- grep(var,chelsa_files,value = T)
  
  mes_año <- str_extract(var_files, "\\d{2}_\\d{4}")
  año_mes <- gsub("(\\d{2})_(\\d{4})", "\\2-\\1", mes_año)
  
  var_r <- rast(var_files) |> 
    crop(chl) |> 
    setNames(año_mes)
  
  out.dir <- glue('data/processed/raster/CHELSA/{var}/')
  
  exportRast(var_r,out.dir,glue('CHELSA_{var}_'))
  
})
