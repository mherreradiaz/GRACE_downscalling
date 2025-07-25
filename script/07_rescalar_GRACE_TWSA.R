library(fs)
library(terra)
library(rnaturalearth)
library(ncdf4)
library(sf)
library(tidyverse)
library(glue)

chl <- ne_countries(country='chile',returnclass = 'sv') |> 
  project('EPSG:4326')

# preparar inputs

grace_files <- dir_ls('data/processed/raster/GRACE_monthly/TWSA/')
fechas_grace <- str_extract(grace_files,'\\d{4}-\\d{2}')

pred_files <- dir_ls('data/processed/raster/downscaling_preds/')
DEM <- rast(pred_files[grep('DEM',pred_files)])
pred_files <- pred_files[-grep('DEM',pred_files)]
pred_names <- str_extract(pred_files,'(?<=preds/).*?(?=.tif)')

minmax_dates <- lapply(pred_files, \(pred) {
  min_date <- rast(pred) |> names() |> min()
  max_date <- rast(pred) |> names() |> max()
  data.frame(min_date = min_date, max_date = max_date)
}) |> 
  setNames(pred_names)

min_date_preds <- map_dfr(minmax_dates, 'min_date') |> unlist() |> max()
max_date_preds <- map_dfr(minmax_dates, 'max_date') |> unlist() |> min()

min_date <- max(c(min_date_preds,min(fechas_grace)))
max_date <- min(c(max_date_preds,max(fechas_grace)))

fechas_modelo <- fechas_grace[between(fechas_grace,min_date,max_date)]

r_preds <- lapply(pred_files, \(pred) {
  pred_r <- rast(pred) |> 
    subset(fechas_modelo)
}) |> 
  setNames(pred_names)

grace_r <- rast(grace_files) |> 
  mask(chl) |> 
  subset(fechas_modelo)

#3. Generar tabla con datos

grid <- grace_r[[1]]
values(grid) <- 1:ncell(grid)
grid_nonNA <- mask(grid, grace_r[[1]])
grid_pol <- as.polygons(grid_nonNA)

sample_pts <- grid_pol |> st_as_sf() |> st_centroid()

library(purrr)
library(bonsai)
library(stacks)
library(rsample)
library(recipes)
library(tidymodels)

# evaluación del desempeño con menores puntos

reduce_dataset <- function(df, p) {
  n <- nrow(df)
  num_remove <- floor(n * p / 100)
  
  if (num_remove == 0) {
    message('Porcentaje demasiado pequeño, no se elimina ninguna fila')
    return(df)
  }
  
  # Generar índices a eliminar distribuidos equitativamente
  idx_remove <- round(seq(1, n, length.out = num_remove + 2))[-c(1, num_remove + 2)]
  
  df_reduced <- df[-idx_remove, ]
  return(df_reduced)
}
sample_evenly <- function(x, n) {
  idx <- round(seq(1, length(x), length.out = n))
  return(x[idx])
}

resultados_pts <- map(sample_evenly(1:nlyr(twsa),4), \(i) {
  
  out_preds <- c(
    setNames(twsa[[i]], 'TWSA'),
    rast(lapply(r_preds, \(pred) pred[[i]])),
    DEM
  )
  
  # Extraer datos a puntos y filtrar NA
  data_modelo <- terra::extract(out_preds, vect(sample_pts)) |>
    select(-ID) |>
    drop_na()
  
  if (nrow(data_modelo) < 10 || length(unique(data_modelo$TWSA)) <= 1) {
    return(met = tibble(
      fecha = as.Date(glue('{fechas_modelo[i]}-01')),
      rsq = NA,
      rmse = NA,
      mae = NA
    ))
  }
  
  # División de datos
  set.seed(456)
  splits <- initial_split(data_modelo)
  data_train <- training(splits)
  data_test <- testing(splits)
  
  # Especificación modelo RF
  rf_spec <- rand_forest(
    trees = 1000,
    mtry = tune(),
    min_n = tune()
  ) |> 
    set_engine('ranger', importance = 'impurity') |> 
    set_mode('regression')
  
  # Recipe de preprocesamiento
  model_rec <- recipe(TWSA ~ ., data = data_train) |>
    step_impute_knn(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_zv(all_numeric_predictors())
  
  # Validación cruzada
  set.seed(453)
  vb_folds <- vfold_cv(data_train)
  ctrl <- control_stack_grid()
  
  twsa_res <- workflow_set(
    preproc = list(rec1 = model_rec),
    models = list(RF = rf_spec)
  ) |> 
    workflow_map(
      seed = 1603,
      resamples = vb_folds,
      grid = 10,
      metrics = metric_set(rsq, rmse, mae),
      control = ctrl,
    )
  
  # Resultados de folds
  fold_results <- twsa_res |>
    extract_workflow_set_result('rec1_RF') |>
    collect_metrics() |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01')))
  
  # Ajuste final
  model_last_fit <- twsa_res |>
    extract_workflow('rec1_RF') |>
    finalize_workflow(
      twsa_res |>
        extract_workflow_set_result('rec1_RF') |>
        select_best(metric = 'rsq')
    ) |>
    last_fit(split = splits, metrics = metric_set(rsq, rmse, mae))
  
  # Guardar métricas
  met <- model_last_fit |> collect_metrics() |>
    select(.metric, .estimate) |>
    pivot_wider(names_from = .metric, values_from = .estimate) |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01'))) |>
    relocate(fecha)
  
  return(list(
    met = met,
    folds = fold_results
  ))
})
resultados_pts_10 <- map(sample_evenly(1:nlyr(twsa),4), \(i) {
  
  out_preds <- c(
    setNames(twsa[[i]], 'TWSA'),
    rast(lapply(r_preds, \(pred) pred[[i]])),
    DEM
  )
  
  # Extraer datos a puntos y filtrar NA
  data_modelo <- terra::extract(out_preds, vect(sample_pts)) |>
    select(-ID) |>
    drop_na() |> 
    reduce_dataset(10)
  
  if (nrow(data_modelo) < 10 || length(unique(data_modelo$TWSA)) <= 1) {
    return(met = tibble(
      fecha = as.Date(glue('{fechas_modelo[i]}-01')),
      rsq = NA,
      rmse = NA,
      mae = NA
    ))
  }
  
  # División de datos
  set.seed(456)
  splits <- initial_split(data_modelo)
  data_train <- training(splits)
  data_test <- testing(splits)
  
  # Especificación modelo RF
  rf_spec <- rand_forest(
    trees = 1000,
    mtry = tune(),
    min_n = tune()
  ) |> 
    set_engine('ranger', importance = 'impurity') |> 
    set_mode('regression')
  
  # Recipe de preprocesamiento
  model_rec <- recipe(TWSA ~ ., data = data_train) |>
    step_impute_knn(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_zv(all_numeric_predictors())
  
  # Validación cruzada
  set.seed(453)
  vb_folds <- vfold_cv(data_train)
  ctrl <- control_stack_grid()
  
  twsa_res <- workflow_set(
    preproc = list(rec1 = model_rec),
    models = list(RF = rf_spec)
  ) |> 
    workflow_map(
      seed = 1603,
      resamples = vb_folds,
      grid = 10,
      metrics = metric_set(rsq, rmse, mae),
      control = ctrl,
    )
  
  # Resultados de folds
  fold_results <- twsa_res |>
    extract_workflow_set_result('rec1_RF') |>
    collect_metrics() |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01')))
  
  # Ajuste final
  model_last_fit <- twsa_res |>
    extract_workflow('rec1_RF') |>
    finalize_workflow(
      twsa_res |>
        extract_workflow_set_result('rec1_RF') |>
        select_best(metric = 'rsq')
    ) |>
    last_fit(split = splits, metrics = metric_set(rsq, rmse, mae))
  
  # Guardar métricas
  met <- model_last_fit |> collect_metrics() |>
    select(.metric, .estimate) |>
    pivot_wider(names_from = .metric, values_from = .estimate) |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01'))) |>
    relocate(fecha)
  
  return(list(
    met = met,
    folds = fold_results
  ))
})
resultados_pts_20 <- map(sample_evenly(1:nlyr(twsa),4), \(i) {
  
  out_preds <- c(
    setNames(twsa[[i]], 'TWSA'),
    rast(lapply(r_preds, \(pred) pred[[i]])),
    DEM
  )
  
  # Extraer datos a puntos y filtrar NA
  data_modelo <- terra::extract(out_preds, vect(sample_pts)) |>
    select(-ID) |>
    drop_na() |> 
    reduce_dataset(20)
  
  if (nrow(data_modelo) < 10 || length(unique(data_modelo$TWSA)) <= 1) {
    return(met = tibble(
      fecha = as.Date(glue('{fechas_modelo[i]}-01')),
      rsq = NA,
      rmse = NA,
      mae = NA
    ))
  }
  
  # División de datos
  set.seed(456)
  splits <- initial_split(data_modelo)
  data_train <- training(splits)
  data_test <- testing(splits)
  
  # Especificación modelo RF
  rf_spec <- rand_forest(
    trees = 1000,
    mtry = tune(),
    min_n = tune()
  ) |> 
    set_engine('ranger', importance = 'impurity') |> 
    set_mode('regression')
  
  # Recipe de preprocesamiento
  model_rec <- recipe(TWSA ~ ., data = data_train) |>
    step_impute_knn(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_zv(all_numeric_predictors())
  
  # Validación cruzada
  set.seed(453)
  vb_folds <- vfold_cv(data_train)
  ctrl <- control_stack_grid()
  
  twsa_res <- workflow_set(
    preproc = list(rec1 = model_rec),
    models = list(RF = rf_spec)
  ) |> 
    workflow_map(
      seed = 1603,
      resamples = vb_folds,
      grid = 10,
      metrics = metric_set(rsq, rmse, mae),
      control = ctrl,
    )
  
  # Resultados de folds
  fold_results <- twsa_res |>
    extract_workflow_set_result('rec1_RF') |>
    collect_metrics() |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01')))
  
  # Ajuste final
  model_last_fit <- twsa_res |>
    extract_workflow('rec1_RF') |>
    finalize_workflow(
      twsa_res |>
        extract_workflow_set_result('rec1_RF') |>
        select_best(metric = 'rsq')
    ) |>
    last_fit(split = splits, metrics = metric_set(rsq, rmse, mae))
  
  # Guardar métricas
  met <- model_last_fit |> collect_metrics() |>
    select(.metric, .estimate) |>
    pivot_wider(names_from = .metric, values_from = .estimate) |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01'))) |>
    relocate(fecha)
  
  return(list(
    met = met,
    folds = fold_results
  ))
})
resultados_pts_30 <- map(sample_evenly(1:nlyr(twsa),4), \(i) {
  
  out_preds <- c(
    setNames(twsa[[i]], 'TWSA'),
    rast(lapply(r_preds, \(pred) pred[[i]])),
    DEM
  )
  
  # Extraer datos a puntos y filtrar NA
  data_modelo <- terra::extract(out_preds, vect(sample_pts)) |>
    select(-ID) |>
    drop_na() |> 
    reduce_dataset(30)
  
  if (nrow(data_modelo) < 10 || length(unique(data_modelo$TWSA)) <= 1) {
    return(met = tibble(
      fecha = as.Date(glue('{fechas_modelo[i]}-01')),
      rsq = NA,
      rmse = NA,
      mae = NA
    ))
  }
  
  # División de datos
  set.seed(456)
  splits <- initial_split(data_modelo)
  data_train <- training(splits)
  data_test <- testing(splits)
  
  # Especificación modelo RF
  rf_spec <- rand_forest(
    trees = 1000,
    mtry = tune(),
    min_n = tune()
  ) |> 
    set_engine('ranger', importance = 'impurity') |> 
    set_mode('regression')
  
  # Recipe de preprocesamiento
  model_rec <- recipe(TWSA ~ ., data = data_train) |>
    step_impute_knn(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_zv(all_numeric_predictors())
  
  # Validación cruzada
  set.seed(453)
  vb_folds <- vfold_cv(data_train)
  ctrl <- control_stack_grid()
  
  twsa_res <- workflow_set(
    preproc = list(rec1 = model_rec),
    models = list(RF = rf_spec)
  ) |> 
    workflow_map(
      seed = 1603,
      resamples = vb_folds,
      grid = 10,
      metrics = metric_set(rsq, rmse, mae),
      control = ctrl,
    )
  
  # Resultados de folds
  fold_results <- twsa_res |>
    extract_workflow_set_result('rec1_RF') |>
    collect_metrics() |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01')))
  
  # Ajuste final
  model_last_fit <- twsa_res |>
    extract_workflow('rec1_RF') |>
    finalize_workflow(
      twsa_res |>
        extract_workflow_set_result('rec1_RF') |>
        select_best(metric = 'rsq')
    ) |>
    last_fit(split = splits, metrics = metric_set(rsq, rmse, mae))
  
  # Guardar métricas
  met <- model_last_fit |> collect_metrics() |>
    select(.metric, .estimate) |>
    pivot_wider(names_from = .metric, values_from = .estimate) |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01'))) |>
    relocate(fecha)
  
  return(list(
    met = met,
    folds = fold_results
  ))
})
resultados_pts_40 <- map(sample_evenly(1:nlyr(twsa),4), \(i) {
  
  out_preds <- c(
    setNames(twsa[[i]], 'TWSA'),
    rast(lapply(r_preds, \(pred) pred[[i]])),
    DEM
  )
  
  # Extraer datos a puntos y filtrar NA
  data_modelo <- terra::extract(out_preds, vect(sample_pts)) |>
    select(-ID) |>
    drop_na() |> 
    reduce_dataset(40)
  
  if (nrow(data_modelo) < 10 || length(unique(data_modelo$TWSA)) <= 1) {
    return(met = tibble(
      fecha = as.Date(glue('{fechas_modelo[i]}-01')),
      rsq = NA,
      rmse = NA,
      mae = NA
    ))
  }
  
  # División de datos
  set.seed(456)
  splits <- initial_split(data_modelo)
  data_train <- training(splits)
  data_test <- testing(splits)
  
  # Especificación modelo RF
  rf_spec <- rand_forest(
    trees = 1000,
    mtry = tune(),
    min_n = tune()
  ) |> 
    set_engine('ranger', importance = 'impurity') |> 
    set_mode('regression')
  
  # Recipe de preprocesamiento
  model_rec <- recipe(TWSA ~ ., data = data_train) |>
    step_impute_knn(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_zv(all_numeric_predictors())
  
  # Validación cruzada
  set.seed(453)
  vb_folds <- vfold_cv(data_train)
  ctrl <- control_stack_grid()
  
  twsa_res <- workflow_set(
    preproc = list(rec1 = model_rec),
    models = list(RF = rf_spec)
  ) |> 
    workflow_map(
      seed = 1603,
      resamples = vb_folds,
      grid = 10,
      metrics = metric_set(rsq, rmse, mae),
      control = ctrl,
    )
  
  # Resultados de folds
  fold_results <- twsa_res |>
    extract_workflow_set_result('rec1_RF') |>
    collect_metrics() |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01')))
  
  # Ajuste final
  model_last_fit <- twsa_res |>
    extract_workflow('rec1_RF') |>
    finalize_workflow(
      twsa_res |>
        extract_workflow_set_result('rec1_RF') |>
        select_best(metric = 'rsq')
    ) |>
    last_fit(split = splits, metrics = metric_set(rsq, rmse, mae))
  
  # Guardar métricas
  met <- model_last_fit |> collect_metrics() |>
    select(.metric, .estimate) |>
    pivot_wider(names_from = .metric, values_from = .estimate) |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01'))) |>
    relocate(fecha)
  
  return(list(
    met = met,
    folds = fold_results
  ))
})

lista <-list(resultados_pts,resultados_pts_10,
             resultados_pts_20,resultados_pts_30,
             resultados_pts_40)

folds_pts_df <- lapply(1:5, \(i) {
  map_dfr(lista[[i]], 'folds') |> 
    mutate(pts = paste0(c(10,9,8,7,6)[i],'0%'),
           .before = 1)
}) |> bind_rows() |> 
  mutate(pts = factor(pts,levels = paste0(c(10,9,8,7,6),'0%')),
         .metric = factor(.metric,levels = c('rsq','rmse','mae')))

rsq_pts_df <- lapply(1:5, \(i) {
  map_dfr(lista[[i]], 'met') |> 
    mutate(pts = paste0(c(10,9,8,7,6)[i],'0%'),
           .before = 1)
}) |> bind_rows() |> 
  mutate(pts = factor(pts,levels = paste0(c(10,9,8,7,6),'0%'))) |> 
  pivot_longer(cols = c(rsq,rmse,mae), names_to = '.metric',values_to = 'mean') |> 
  mutate(.metric = factor(.metric,levels = c('rsq','rmse','mae')))

folds_pts_df |> 
  ggplot(aes(as.factor(fecha), mean, fill = pts)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(data = rsq_pts_df,
             aes(x = as.factor(fecha), y = mean, fill = pts),
             position = position_dodge(width = 0.75),
             shape = 21,
             size = 3,
             inherit.aes = FALSE) +
  scale_fill_viridis_d(option = 'D') +
  facet_wrap(~ .metric, scales = 'free_y', ncol = 1) +
  labs(x = NULL, y = "metric's values") +
  theme_bw()

# rescalado

resultados <- map(seq_along(fechas_modelo), \(i) {
  
  # Cargar variable objetivo (TWSA) como raster independiente
  twsa_layer <- setNames(grace_r[[i]], 'TWSA')
  
  # Crear stack de predictores (excluyendo TWSA)
  out_preds <- c(
    rast(lapply(r_preds, \(pred) pred[[i]])),
    setNames(DEM, 'DEM')
  )
  
  # Extraer TWSA a puntos
  twsa_vals <- terra::extract(twsa_layer, vect(sample_pts), cores = 10) |> 
    select(-ID)
  
  # Extraer predictores a puntos
  preds_vals <- terra::extract(out_preds, vect(sample_pts), cores = 10) |>
    select(-ID)
  
  # Combinar en data_modelo
  data_modelo <- bind_cols(twsa_vals, preds_vals) |> 
    drop_na()
  
  if (nrow(data_modelo) < 10 || length(unique(data_modelo$TWSA)) <= 1) {
    return(met = tibble(
      fecha = as.Date(glue('{fechas_modelo[i]}-01')),
      rsq = NA,
      rmse = NA,
      mae = NA
    ))
  }
  
  stats_pred <- data_modelo |>
    summarise(across(everything(), list(
      mean = \(x) mean(x, na.rm = TRUE),
      sd   = \(x) sd(x, na.rm = TRUE)
    ))) |> 
    pivot_longer(cols = everything(), names_to = c('variable', '.value'), names_sep = '_') |> 
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01'))) |> 
    relocate(fecha)
  
  # División de datos
  set.seed(456)
  splits <- initial_split(data_modelo)
  data_train <- training(splits)
  data_test <- testing(splits)
  
  # Especificación modelo RF
  rf_spec <- rand_forest(
    trees = 1000,
    mtry = tune(),
    min_n = tune()
  ) |> 
    set_engine('ranger', importance = 'impurity') |> 
    set_mode('regression')
  
  # Recipe de preprocesamiento
  model_rec <- recipe(TWSA ~ ., data = data_train) |>
    step_impute_knn(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_zv(all_numeric_predictors())
  
  # Validación cruzada
  set.seed(453)
  vb_folds <- vfold_cv(data_train)
  ctrl <- control_stack_grid()
  
  twsa_res <- workflow_set(
    preproc = list(rec1 = model_rec),
    models = list(RF = rf_spec)
  ) |> 
    workflow_map(
      seed = 1603,
      resamples = vb_folds,
      grid = 10,
      metrics = metric_set(rsq, rmse, mae),
      control = ctrl,
    )
  
  # Resultados de folds
  fold_results <- twsa_res |>
    extract_workflow_set_result('rec1_RF') |>
    collect_metrics() |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01')))
  
  # Ajuste final
  model_last_fit <- twsa_res |>
    extract_workflow('rec1_RF') |>
    finalize_workflow(
      twsa_res |>
        extract_workflow_set_result('rec1_RF') |>
        select_best(metric = 'rsq')
    ) |>
    last_fit(split = splits, metrics = metric_set(rsq, rmse, mae))
  
  # Importancia de variables
  vip <- model_last_fit |>
    extract_workflow() |>
    extract_fit_parsnip() |>
    vip::vi() |> 
    select(variable = Variable, importance = Importance) |> 
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01')))
  
  # Guardar métricas
  met <- model_last_fit |> collect_metrics() |>
    select(.metric, .estimate) |>
    pivot_wider(names_from = .metric, values_from = .estimate) |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01'))) |> 
    relocate(fecha)
  
  # Predicción final sobre raster de predictores
  out_preds[is.na(out_preds)] <- -9999
  twsa_preds <- predict(model_last_fit |> extract_workflow(), as.data.frame(out_preds))
  
  # Crear raster de predicciones con la geometría del primer predictor
  out <- out_preds[[1]]
  values(out) <- twsa_preds$.pred
  out <- mask(out, chl)
  
  # Exportar raster de predicciones
  writeRaster(out, glue('data/processed/raster/GRACE_TWSA_highres/GRACE_TWSA_highres_{fechas_modelo[i]}.tif'), 
              overwrite = TRUE)
  
  return(list(
    met = met,
    vip = vip,
    stats = stats_pred,
    folds = fold_results
  ))
})

ds_metrics <- map_dfr(resultados, 'met') |> 
  select(fecha,everything())
ds_folds <- map_dfr(resultados, 'folds') |> 
  select(fecha,everything())
ds_vip <- map_dfr(resultados, 'vip') |> 
  select(fecha,everything())
ds_stats <- map_dfr(resultados, 'stats') |> 
  select(fecha,everything())

write_rds(ds_metrics,'data/processed/rds/GRACE_TWSA_downscaling_metrics.rds')
write_rds(ds_folds,'data/processed/rds/GRACE_TWSA_downscaling_folds.rds')
write_rds(ds_vip,'data/processed/rds/GRACE_TWSA_downscaling_VIP.rds')
write_rds(ds_stats,'data/processed/rds/GRACE_TWSA_downscaling_preds_stats.rds')

# visualizar

library(Kendall)
complete.dates <- function(data, date_col, group_cols = NULL, by = "month") {
  library(dplyr)
  library(tidyr)
  library(rlang)
  
  date_col <- enquo(date_col)
  group_cols <- enquos(group_cols)
  
  dates <- pull(data, !!date_col)
  
  seq_dates <- seq(min(dates, na.rm = TRUE),
                   max(dates, na.rm = TRUE),
                   by = by)
  
  if (length(unique(dates)) == nrow(data) & length(group_cols) > 0) {
    group_cols <- NULL
    cat('✅ No es necesario añadir group_cols. Solo hay un valor por fecha.\n')
  }
  
  if (is.null(group_cols) || length(group_cols) == 0) {
    data_completed <- tibble(!!date_col := seq_dates) |> 
      left_join(data, by = quo_name(date_col)) |> 
      arrange(!!date_col)
    
  } else {
    expanded <- data |> 
      distinct(!!!group_cols) |> 
      crossing(!!date_col := seq_dates)
    
    data_completed <- expanded |> 
      left_join(data, by = c(setNames(sapply(group_cols, quo_name), sapply(group_cols, quo_name)),
                             quo_name(date_col))) |> 
      arrange(!!!group_cols, !!date_col)
  }
  
  return(data_completed)
}

data_metrics <- read_rds('data/processed/rds/GRACE_TWSA_downscaling_metrics.rds')
data_folds <- read_rds('data/processed/rds/GRACE_TWSA_downscaling_folds.rds')
data_vip <- read_rds('data/processed/rds/GRACE_TWSA_downscaling_VIP.rds')
data_stats <- read_rds('data/processed/rds/GRACE_TWSA_downscaling_preds_stats.rds')

# métricas

data_metrics |> 
  reframe(mean_rsq = mean(rsq, na.rm=T),
          sd_rsq = sd(rsq, na.rm=T),
          mean_rmse = mean(rmse, na.rm=T),
          sd_rmse = sd(rmse, na.rm=T),
          mean_mae = mean(mae, na.rm=T),
          sd_mae = sd(mae, na.rm=T),
          trend_rsq = as.numeric(MannKendall(rsq)$tau),
          trend_rmse = as.numeric(MannKendall(rmse)$tau),
          trend_mae = as.numeric(MannKendall(mae)$tau))

data_metrics |>
  rename(RSQ = rsq, RMSE = rmse,MAE = mae) |> 
  pivot_longer(cols = c(MAE, RMSE, RSQ), names_to = 'metrica', values_to = 'valor') |>
  mutate(metrica = factor(metrica,levels=c('RSQ','RMSE','MAE'))) |> 
  complete.dates(fecha,metrica) |> 
  ggplot(aes(x = as.Date(fecha), y = valor)) +
  geom_line(na.rm = T,linewidth = .5,alpha = .8, color = 'darkcyan') +
  geom_line(stat = "smooth", method = "loess", span = 0.3,
            linewidth = .5, alpha = 0.8, color = 'brown2') +
  facet_wrap(~metrica, scales = 'free_y', ncol = 1,strip.position = 'left') +
  labs(x = NULL,y = NULL) +
  scale_x_date(breaks = seq(as.Date('2005-01-01'), as.Date('2023-01-01'), by = '2 years'),
               minor_breaks = seq(as.Date('2000-01-01'), as.Date('2024-01-01'), by = '1 year'),
               date_labels = '%Y',
               expand = c(.01,0)) +
  theme_bw() +
  theme(strip.background = element_rect(fill='white'))

ggsave('output/fig/downscaling/GRACE_TWSA/model_metrics.png',width = 10, height = 6)

# folds

data_folds |> 
  rename(value = mean) |> 
  group_by(fecha, .metric) |> 
  reframe(mean = mean(value),
          se = sd(value) / sqrt(length(value))) |> 
  mutate(.metric = factor(toupper(.metric), levels = c('RSQ','RMSE','MAE'))) |> 
  complete.dates(fecha,.metric) |> 
  ggplot(aes(x = as.Date(fecha), y = mean)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), fill = 'grey70', alpha = 0.5) +  # sombra con SE
  geom_line(na.rm = TRUE, linewidth = 0.5, alpha = 0.8, color = 'darkcyan') +
  geom_line(stat = "smooth", method = "loess", span = 0.3,
            linewidth = .5, alpha = 0.8, color = 'brown2') +
  facet_wrap(~ .metric, scales = 'free_y', ncol = 1, strip.position = 'left') +
  labs(x = NULL, y = NULL) +
  scale_x_date(breaks = seq(as.Date('2005-01-01'), as.Date('2023-01-01'), by = '2 years'),
               minor_breaks = seq(as.Date('2000-01-01'), as.Date('2024-01-01'), by = '1 year'),
               date_labels = '%Y',
               expand = c(.01, 0)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = 'white'))

ggsave('output/fig/downscaling/GRACE_TWSA/model_folds_metrics.png',width = 10, height = 6)

# vip

data_vip |>
  group_by(fecha) |> 
  slice_max(importance, n = 1) |> 
  group_by(variable) |> 
  summarise(freq = n()) |> 
  ggplot(aes(x = reorder(variable, freq), y = freq)) +
  geom_col(fill = 'darkcyan',color = 'darkslategrey') +
  coord_flip() +
  scale_y_continuous(expand = c(0,0),limits = c(0,128)) +
  scale_x_discrete(expand = c(.062,0)) +
  labs(x = NULL, y = 'count of months as top predictor') +
  theme_bw()

ggsave('output/fig/downscaling/GRACE_TWSA/model_vip_count.png',width = 10, height = 6)

imp_var <- data_vip |>
  group_by(fecha) |> 
  slice_max(importance, n = 1) |>
  group_by(variable) |> 
  reframe(n = n()) |> 
  arrange(n) |> 
  pull(variable)

pal <- c('orange2','navajowhite','brown2','darkslategrey','darkcyan')

data_vip |>
  group_by(fecha) |> 
  slice_max(importance, n = 1) |>
  mutate(variable = factor(variable,levels = imp_var)) |> 
  ggplot(aes(x = fecha, y = variable, fill = variable)) +
  # geom_point(size = 2, shape = 21) +
  geom_tile(fill = 'darkcyan') +
  labs(x = NULL, y = NULL) +
  scale_x_date(breaks = seq(as.Date('2005-01-01'), as.Date('2023-01-01'), by = '2 years'),
               minor_breaks = seq(as.Date('2000-01-01'), as.Date('2024-01-01'), by = '1 year'),
               date_labels = '%Y',
               expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(legend.position = 'none')

ggsave('output/fig/downscaling/GRACE_TWSA/model_vip_month.png',width = 12, height = 5)

# stats predictores

data_stats |>
  filter(variable %in% c(imp_var,'TWSA')) |> 
  mutate(variable = factor(variable,levels = rev(c(imp_var,'TWSA')))) |> 
  complete.dates(fecha,variable) |> 
  ggplot(aes(as.Date(fecha),mean)) +
  geom_line(na.rm = T,linewidth = .5,alpha = .8, color = 'darkcyan') +
  geom_line(stat = "smooth", method = "loess", span = 0.3,
            linewidth = .5, alpha = 0.8, color = 'brown2') +
  facet_wrap(~variable, scales = 'free_y',ncol = 2) +
  scale_x_date(breaks = seq(as.Date('2005-01-01'), as.Date('2023-01-01'), by = '2 years'),
               minor_breaks = seq(as.Date('2000-01-01'), as.Date('2024-01-01'), by = '1 year'),
               date_labels = '%Y',
               expand = c(.01,0)) +
  labs(x = NULL) +
  theme_bw() +
  theme(strip.background = element_rect(fill='white'))

ggsave('output/fig/downscaling/GRACE_TWSA/model_stats_mean.png',width = 13, height = 7)

data_stats |>
  filter(variable %in% c(imp_var,'TWSA')) |> 
  mutate(variable = factor(variable,levels = rev(c(imp_var,'TWSA')))) |> 
  complete.dates(fecha,variable) |> 
  ggplot(aes(as.Date(fecha),sd)) +
  geom_line(na.rm = T,linewidth = .5,alpha = .8, color = 'darkcyan') +
  geom_line(stat = "smooth", method = "loess", span = 0.3,
            linewidth = .5, alpha = 0.8, color = 'brown2') +
  facet_wrap(~variable, scales = 'free_y',ncol = 2) +
  scale_x_date(breaks = seq(as.Date('2005-01-01'), as.Date('2023-01-01'), by = '2 years'),
               minor_breaks = seq(as.Date('2000-01-01'), as.Date('2024-01-01'), by = '1 year'),
               date_labels = '%Y',
               expand = c(.01,0)) +
  labs(x = NULL) +
  theme_bw() +
  theme(strip.background = element_rect(fill='white'))

ggsave('output/fig/downscaling/GRACE_TWSA/model_stats_sd.png',width = 13, height = 7)
