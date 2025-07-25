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

grace_files <- dir_ls('data/processed/raster/GRACE_monthly/GWSA/')
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

# generar tabla con datos

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

# rescalar

resultados <- map(seq_along(fechas_modelo), \(i) {
  
  # Cargar variable objetivo (GWSA) como raster independiente
  gwsa_layer <- setNames(grace_r[[i]], 'GWSA')
  
  # Crear stack de predictores (excluyendo GWSA)
  out_preds <- c(
    rast(lapply(r_preds, \(pred) pred[[i]])),
    setNames(DEM, 'DEM')
  )
  
  # Extraer GWSA a puntos
  gwsa_vals <- terra::extract(gwsa_layer, vect(sample_pts), cores = 10) |> 
    select(-ID)
  
  # Extraer predictores a puntos
  preds_vals <- terra::extract(out_preds, vect(sample_pts), cores = 10) |>
    select(-ID)
  
  # Combinar en data_modelo
  data_modelo <- bind_cols(gwsa_vals, preds_vals) |> 
    drop_na()
  
  if (nrow(data_modelo) < 10 || length(unique(data_modelo$GWSA)) <= 1) {
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
  model_rec <- recipe(GWSA ~ ., data = data_train) |>
    step_impute_knn(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_zv(all_numeric_predictors())
  
  # Validación cruzada
  set.seed(453)
  vb_folds <- vfold_cv(data_train)
  ctrl <- control_stack_grid()
  
  gwsa_res <- workflow_set(
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
  fold_results <- gwsa_res |>
    extract_workflow_set_result('rec1_RF') |>
    collect_metrics() |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01')))
  
  # Ajuste final
  model_last_fit <- gwsa_res |>
    extract_workflow('rec1_RF') |>
    finalize_workflow(
      gwsa_res |>
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
  gwsa_preds <- predict(model_last_fit |> extract_workflow(), as.data.frame(out_preds))
  
  # Crear raster de predicciones con la geometría del primer predictor
  out <- out_preds[[1]]
  values(out) <- gwsa_preds$.pred
  out <- mask(out, chl)
  
  # Exportar raster de predicciones
  writeRaster(out, glue('data/processed/raster/GRACE_GWSA_highres/GRACE_GWSA_highres_{fechas_modelo[i]}.tif'), 
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

write_rds(ds_metrics,'data/processed/rds/GRACE_GWSA_downscaling_metrics.rds')
write_rds(ds_folds,'data/processed/rds/GRACE_GWSA_downscaling_folds.rds')
write_rds(ds_vip,'data/processed/rds/GRACE_GWSA_downscaling_VIP.rds')
write_rds(ds_stats,'data/processed/rds/GRACE_GWSA_downscaling_preds_stats.rds')

# evaluar desempeño con menores puntos

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

resultados_pts <- map(sample_evenly(1:nlyr(gwsa),4), \(i) {
  
  out_preds <- c(
    setNames(gwsa[[i]], 'GWSA'),
    rast(lapply(r_preds, \(pred) pred[[i]])),
    DEM
  )
  
  # Extraer datos a puntos y filtrar NA
  data_modelo <- terra::extract(out_preds, vect(sample_pts)) |>
    select(-ID) |>
    drop_na()
  
  if (nrow(data_modelo) < 10 || length(unique(data_modelo$GWSA)) <= 1) {
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
  model_rec <- recipe(GWSA ~ ., data = data_train) |>
    step_impute_knn(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_zv(all_numeric_predictors())
  
  # Validación cruzada
  set.seed(453)
  vb_folds <- vfold_cv(data_train)
  ctrl <- control_stack_grid()
  
  gwsa_res <- workflow_set(
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
  fold_results <- gwsa_res |>
    extract_workflow_set_result('rec1_RF') |>
    collect_metrics() |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01')))
  
  # Ajuste final
  model_last_fit <- gwsa_res |>
    extract_workflow('rec1_RF') |>
    finalize_workflow(
      gwsa_res |>
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
resultados_pts_10 <- map(sample_evenly(1:nlyr(gwsa),4), \(i) {
  
  out_preds <- c(
    setNames(gwsa[[i]], 'GWSA'),
    rast(lapply(r_preds, \(pred) pred[[i]])),
    DEM
  )
  
  # Extraer datos a puntos y filtrar NA
  data_modelo <- terra::extract(out_preds, vect(sample_pts)) |>
    select(-ID) |>
    drop_na() |> 
    reduce_dataset(10)
  
  if (nrow(data_modelo) < 10 || length(unique(data_modelo$GWSA)) <= 1) {
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
  model_rec <- recipe(GWSA ~ ., data = data_train) |>
    step_impute_knn(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_zv(all_numeric_predictors())
  
  # Validación cruzada
  set.seed(453)
  vb_folds <- vfold_cv(data_train)
  ctrl <- control_stack_grid()
  
  gwsa_res <- workflow_set(
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
  fold_results <- gwsa_res |>
    extract_workflow_set_result('rec1_RF') |>
    collect_metrics() |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01')))
  
  # Ajuste final
  model_last_fit <- gwsa_res |>
    extract_workflow('rec1_RF') |>
    finalize_workflow(
      gwsa_res |>
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
resultados_pts_20 <- map(sample_evenly(1:nlyr(gwsa),4), \(i) {
  
  out_preds <- c(
    setNames(gwsa[[i]], 'GWSA'),
    rast(lapply(r_preds, \(pred) pred[[i]])),
    DEM
  )
  
  # Extraer datos a puntos y filtrar NA
  data_modelo <- terra::extract(out_preds, vect(sample_pts)) |>
    select(-ID) |>
    drop_na() |> 
    reduce_dataset(20)
  
  if (nrow(data_modelo) < 10 || length(unique(data_modelo$GWSA)) <= 1) {
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
  model_rec <- recipe(GWSA ~ ., data = data_train) |>
    step_impute_knn(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_zv(all_numeric_predictors())
  
  # Validación cruzada
  set.seed(453)
  vb_folds <- vfold_cv(data_train)
  ctrl <- control_stack_grid()
  
  gwsa_res <- workflow_set(
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
  fold_results <- gwsa_res |>
    extract_workflow_set_result('rec1_RF') |>
    collect_metrics() |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01')))
  
  # Ajuste final
  model_last_fit <- gwsa_res |>
    extract_workflow('rec1_RF') |>
    finalize_workflow(
      gwsa_res |>
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
resultados_pts_30 <- map(sample_evenly(1:nlyr(gwsa),4), \(i) {
  
  out_preds <- c(
    setNames(gwsa[[i]], 'GWSA'),
    rast(lapply(r_preds, \(pred) pred[[i]])),
    DEM
  )
  
  # Extraer datos a puntos y filtrar NA
  data_modelo <- terra::extract(out_preds, vect(sample_pts)) |>
    select(-ID) |>
    drop_na() |> 
    reduce_dataset(30)
  
  if (nrow(data_modelo) < 10 || length(unique(data_modelo$GWSA)) <= 1) {
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
  model_rec <- recipe(GWSA ~ ., data = data_train) |>
    step_impute_knn(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_zv(all_numeric_predictors())
  
  # Validación cruzada
  set.seed(453)
  vb_folds <- vfold_cv(data_train)
  ctrl <- control_stack_grid()
  
  gwsa_res <- workflow_set(
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
  fold_results <- gwsa_res |>
    extract_workflow_set_result('rec1_RF') |>
    collect_metrics() |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01')))
  
  # Ajuste final
  model_last_fit <- gwsa_res |>
    extract_workflow('rec1_RF') |>
    finalize_workflow(
      gwsa_res |>
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
resultados_pts_40 <- map(sample_evenly(1:nlyr(gwsa),4), \(i) {
  
  out_preds <- c(
    setNames(gwsa[[i]], 'GWSA'),
    rast(lapply(r_preds, \(pred) pred[[i]])),
    DEM
  )
  
  # Extraer datos a puntos y filtrar NA
  data_modelo <- terra::extract(out_preds, vect(sample_pts)) |>
    select(-ID) |>
    drop_na() |> 
    reduce_dataset(40)
  
  if (nrow(data_modelo) < 10 || length(unique(data_modelo$GWSA)) <= 1) {
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
  model_rec <- recipe(GWSA ~ ., data = data_train) |>
    step_impute_knn(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_zv(all_numeric_predictors())
  
  # Validación cruzada
  set.seed(453)
  vb_folds <- vfold_cv(data_train)
  ctrl <- control_stack_grid()
  
  gwsa_res <- workflow_set(
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
  fold_results <- gwsa_res |>
    extract_workflow_set_result('rec1_RF') |>
    collect_metrics() |>
    mutate(fecha = as.Date(glue('{fechas_modelo[i]}-01')))
  
  # Ajuste final
  model_last_fit <- gwsa_res |>
    extract_workflow('rec1_RF') |>
    finalize_workflow(
      gwsa_res |>
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

