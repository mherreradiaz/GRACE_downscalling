library(tidyverse)
library(glue)
library(Kendall)

complete.dates <- function(data, date_col, ..., by = "month") {
  library(dplyr)
  library(tidyr)
  library(rlang)
  
  date_col <- enquo(date_col)
  group_cols <- enquos(...)
  
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

vars_comp <- c("GLDAS_GWSA", "GRACE_TWSA", "GRACE_GWSA")

resultados <- lapply(vars_comp,\(var) {
  metrics <- read_rds(glue('data/processed/rds/{var}_downscaling_metrics.rds')) |> 
    mutate(variable = var,
           .before = 1)
  folds <- read_rds(glue('data/processed/rds/{var}_downscaling_folds.rds')) |> 
    mutate(variable = var,
           .before = 1)
  vip <- read_rds(glue('data/processed/rds/{var}_downscaling_VIP.rds')) |> 
    mutate(variable = var,
           .before = 1)
  stats <- read_rds(glue('data/processed/rds/{var}_downscaling_preds_stats.rds')) |> 
    mutate(variable = var,
           .before = 1)
  return(list(metrics = metrics,folds = folds,vip = vip,stats = stats))
})

data_metrics <- map_dfr(resultados, 'metrics') |> 
  select(fecha,everything()) |> 
  mutate(variable == ifelse(variable == 'GLDAS_GWS','GLDAS_GWSA',variable),
         rmse = ifelse(variable =='GLDAS_GWSA',rmse/10,rmse),
         mae = ifelse(variable =='GLDAS_GWSA',mae/10,mae))
data_folds <- map_dfr(resultados, 'folds') |> 
  select(fecha,everything()) |> 
  mutate(variable == ifelse(variable == 'GLDAS_GWS','GLDAS_GWSA',variable),,
         mean = ifelse(variable =='GLDAS_GWSA' & (.metric == 'rmse' | .metric == 'mae'),mean/10,mean))
data_vip <- map_dfr(resultados, 'vip') |> 
  select(fecha,everything())

# metricas generales

data_metrics |> 
  group_by(variable) |> 
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
  rename(RSQ = rsq, RMSE = rmse, MAE = mae) |> 
  pivot_longer(cols = c(MAE, RMSE, RSQ), names_to = 'metrica', values_to = 'valor') |>
  mutate(metrica = factor(metrica, levels = c('RSQ', 'RMSE', 'MAE'))) |> 
  complete.dates(fecha, variable, metrica) |> 
  ggplot(aes(x = as.Date(fecha), y = valor)) +
  geom_line(na.rm = TRUE, linewidth = .4, alpha = .8, color = 'darkcyan') +
  geom_line(stat = "smooth", method = "loess", span = 0.3,
            linewidth = .4, alpha = 0.8, color = 'brown2') +
  facet_grid2(rows = vars(metrica), cols = vars(variable),
              scales = 'free') +
  labs(x = NULL, y = NULL) +
  scale_x_date(date_labels = "%Y", breaks = seq(as.Date('2003-01-01'), as.Date('2023-01-01'), by = '4 years'),
               limits = c(as.Date('2002-07-01'),as.Date('2024-01-01')),expand = c(.01, 0)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = 'white'))

ggsave('output/fig/downscaling/model_metrics.png',width = 8, height = 5)

pal <- c('orange2','navajowhite','brown2','darkslategrey','darkcyan')

data_metrics |>
  rename(RSQ = rsq, RMSE = rmse, MAE = mae) |> 
  pivot_longer(cols = c(MAE, RMSE, RSQ), names_to = 'metrica', values_to = 'valor') |>
  mutate(metrica = factor(metrica, levels = c('RSQ', 'RMSE', 'MAE'))) |> 
  complete.dates(fecha, variable, metrica) |> 
  ggplot(aes(fecha,valor, color = variable)) +
  geom_line(na.rm = TRUE, linewidth = .6, alpha = .5) +
  geom_line(stat = "smooth", method = "loess", span = 0.2,
            linewidth = .6, linetype = 'dashed') +
  facet_wrap(~metrica,scales = 'free_y',ncol=1, strip.position = 'left') +
  labs(x = NULL, y = NULL, color = NULL) +
  scale_x_date(date_labels = "%Y", breaks = seq(as.Date('2003-01-01'), as.Date('2023-01-01'), by = '4 years'),
               limits = c(as.Date('2002-07-01'),as.Date('2024-01-01')),expand = c(.01, 0)) +
  # scale_color_manual(values = pal[c(1,3,5)]) +
  theme_bw() +
  theme(strip.background = element_rect(fill = 'white'),
        legend.position = 'bottom')

ggsave('output/fig/downscaling/model_metrics_alt.png',width = 8, height = 6)

# métricas

data_metrics |>
  rename(RSQ = rsq, RMSE = rmse, MAE = mae) |> 
  pivot_longer(cols = c(MAE, RMSE, RSQ), names_to = 'metrica', values_to = 'valor') |>
  mutate(metrica = factor(metrica, levels = c('RSQ', 'RMSE', 'MAE'))) |> 
  complete.dates(fecha, variable, metrica) |> 
  ggplot(aes(x = as.Date(fecha), y = valor)) +
  geom_line(na.rm = TRUE, linewidth = .5, alpha = .8, color = 'darkcyan') +
  geom_line(stat = "smooth", method = "loess", span = 0.3,
            linewidth = .5, alpha = 0.8, color = 'brown2') +
  facet_wrap(metrica ~ variable, scales = "free_y",axes = 'all_y',axis.labels = 'margins') +
  labs(x = NULL, y = NULL) +
  scale_x_date(breaks = seq(as.Date('2003-01-01'), as.Date('2023-01-01'), by = '2 years'),
               minor_breaks = seq(as.Date('2000-01-01'), as.Date('2024-01-01'), by = '1 year'),
               date_labels = '%Y',
               expand = c(.01, 0)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = 'white'))

ggsave('output/fig/downscaling/GLDAS_GWSA/model_metrics.png',width = 10, height = 6)

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

ggsave('output/fig/downscaling/GLDAS_GWSA/model_folds_metrics.png',width = 10, height = 6)

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
  scale_x_discrete(expand = c(.11,0)) +
  labs(x = NULL, y = 'count of months as top predictor') +
  theme_bw()

ggsave('output/fig/downscaling/GLDAS_GWSA/model_vip_count.png',width = 10, height = 6)

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

ggsave('output/fig/downscaling/GLDAS_GWSA/model_vip_month.png',width = 12, height = 5)

# stats predictores

data_stats |>
  filter(variable %in% c(imp_var,'GWS')) |> 
  mutate(variable = factor(variable,levels = rev(c(imp_var,'GWS')))) |> 
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

ggsave('output/fig/downscaling/GLDAS_GWSA/model_stats_mean.png',width = 13, height = 7)

data_stats |>
  filter(variable %in% c(imp_var,'GWS')) |> 
  mutate(variable = factor(variable,levels = rev(c(imp_var,'GWS')))) |> 
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

ggsave('output/fig/downscaling/GLDAS_GWSA/model_stats_sd.png',width = 13, height = 7)
