library(terra)
library(sf)
library(tidyverse)
library(glue)

pe <- function(x){
  if (!all(is.na(x))){
    r <- rank(x, na.last = 'keep')
    out <- (r-0.33)/(sum(!is.na(x))+0.33)
  } else out <- rep(NA,length(x))
  out
}
eddi <- function(x){
  if (!all(is.na(x))){
    C0 = 2.515517;C1 = 0.802853;C2 = 0.010328
    d1 = 1.432788;d2 = 0.189269;d3 = 0.001308
    pe <- pe(x)
    
    data <- data.frame(x=x,pe=pe)
    
    out <- data |> 
      dplyr::mutate(W = case_when(pe <= .5 ~ sqrt(-2*log(1-pe)),
                                  pe > .5 ~ sqrt(-2*log(pe))),
                    EDDI = W - (C0+C1*W+C2*W^2)/(1+d1*W+d2*W^2+d3*W^3),
                    EDDI = dplyr::case_when(pe <= .5 ~ EDDI,
                                            pe > .5 ~ -EDDI)
      ) 
  } else out <- data.frame(EDDI = rep(NA,length(x)))
  
  return(out)
}
windowStd <- \(x, scale = 36, minPred = 1) {
  require(dplyr)
  require(zoo)
  
  if (length(x) < scale | all(is.na(x))) {return(rep(NA,length(x)))}
  
  win_accum <- lapply(length(x):36, \(i) {
    idx <- i:(i-scale+1)
    win_i <- x[idx]
    na_idx <- which(is.na(win_i))
    if(scale-length(na_idx) < scale*minPred) return(NA)
    win_i_filled <- zoo::na.approx(win_i, na.rm = F) |> 
      zoo::na.locf(na.rm = F) |> 
      zoo::na.locf(fromLast = T, na.rm = F)
    return(sum(win_i_filled))
  }) |> 
    unlist()
  
  x_accum <- rep(NA,length(x))
  x_accum[length(x):36] <- win_accum
  x_accum[which(is.na(x))] <- NA
  
  return(eddi(x_accum)$EDDI)
}
climAnomaly <- function(x) {
  if (all(is.na(x))) return(rep(NA, length(x)))
  
  meses <- lubridate::month(fechas)
  
  df <- tibble(x = x, mes = meses)
  
  df |> 
    group_by(mes) |> 
    mutate(
      pe = {
        r <- rank(x, na.last = 'keep')
        (r - 0.33) / (sum(!is.na(x)) + 0.33)
      },
      W = case_when(pe <= .5 ~ sqrt(-2 * log(1 - pe)),
                    pe > .5  ~ sqrt(-2 * log(pe))),
      anomaly = W - (2.515517 + 0.802853 * W + 0.010328 * W^2) / 
        (1 + 1.432788 * W + 0.189269 * W^2 + 0.001308 * W^3),
      anomaly = case_when(pe <= .5 ~ anomaly,
                          pe > .5  ~ -anomaly)
    ) |> 
    pull(anomaly)
}

gldas_files <- dir_ls('data/processed/raster/GLDAS/')
fechas <- str_extract(gldas_files,'(?<=DAS_).*?(?=\\.tif)')
gldas_r <- gldas_files |> 
  rast() |> 
  setNames(fechas)

fechas <- as.Date(paste0(names(gldas_r), '-01'))
gldas_climAnom <- terra::app(gldas_r, fun = climAnomaly)

exportRast(gldas_climAnom,'data/processed/raster/GLDAS_anomaly/monthlyAnom_GLDAS_')

v <- c(setNames(setValues(gldas_r[[1]],1:ncell(gldas_r)),'id'),gldas_r)[sample(1:ncell(gldas_r),300)] |> 
  as_tibble() |> 
  filter(!if_all(c(everything(),-id), is.na)) |> 
  slice(1) |> 
  as.numeric()

id <- v[1]
v <- v[-1]

v_eddi <- gldas_climAnom[id] |> as_tibble() |> filter(!if_all(everything(), is.na)) |> slice(1) |> 
  as.numeric()

data_plot <- tibble(fecha = fechas,gws = v, gws_anomaly = v_eddi) |> 
  pivot_longer(cols= c(gws,gws_anomaly),names_to = 'var',values_to = 'value') |> 
  mutate(fecha = as.Date(paste0(fecha,'-01')))

data_plot |> 
  ggplot(aes(fecha,value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~var, ncol = 1,scales = 'free_y', strip.position = "right") +
  labs(x = NULL,y = NULL) +
  scale_x_date(breaks = seq(as.Date('2003-01-01'),as.Date('2025-12-01'),by='2 year'),
               minor_breaks = seq(as.Date('2003-01-01'),as.Date('2025-12-01'),by='1 year'),
               date_labels = '%Y') +
  theme_bw()

data_plot |> 
  filter(var == 'gws_anomaly') |> 
  ggplot(aes(fecha,value)) +
  geom_point() +
  geom_line()
