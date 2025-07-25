library(terra)
library(tidyterra)
library(tidyverse)
library(glue)

data_GWL <- read_rds('data/processed/rds/GWL_chile.rds')
pozos <- vect('data/processed/vectorial/pozos_chile.gpkg')

product_names <- c('GLDAS_GWS','GRACE_TWSA','GRACE_GWSA')

product_r <- lapply(product_names,\(product) {
  files <-  dir_ls(glue('data/processed/raster/{product}_highres/'))
  fechas <- str_extract(files,'\\d{4}-\\d{2}')
  r <- rast(files) |> 
    setNames(fechas)
}) |> 
  setNames(product_names)

data_products <- lapply(product_names, \(name) {
  terra::extract(product_r[[name]],pozos) |> 
    mutate(codigo = pozos$codigo,
           .before = ID) |> 
    select(-ID) |> 
    pivot_longer(c(everything(),-codigo), names_to = 'fecha',
                 values_to = name) |> 
    mutate(fecha = as.Date(paste0(fecha,'-01'))) |> 
    select(fecha,codigo,all_of(name))
}) |> 
  reduce(full_join,by = c('fecha', 'codigo'))

data <- data_products |> 
  left_join(data_GWL) |> 
  group_by(cuenca,codigo) |> 
  mutate(GWL = as.numeric(scale(GWL))) |> 
  ungroup() |> 
  select(fecha,cuenca,codigo,GWL,everything()) |> 
  arrange(fecha,codigo)

write_rds(data,'data/processed/rds/validation_data.rds')

# correlacion

data_raw <- read_rds('data/processed/rds/validation_data.rds')
vars_comp <- c('GLDAS_GWS','GRACE_TWSA','GRACE_GWSA')

cor_general <- map_dfr(vars_comp, ~ {
  tibble(
    variable = .x,
    r = cor(data$GWL, data[[.x]], method = "pearson", use = "complete.obs")
  )
})

cor_well <- map_dfr(vars_comp, ~ {
  data |> 
    group_by(codigo) |> 
    reframe({
      n <- sum(complete.cases(GWL, .data[[.x]]))
      if (n > 3) {
        test <- cor.test(GWL, .data[[.x]], method = 'pearson', use = 'complete.obs')
        tibble(variable = .x,
               r = test$estimate,
               p_value = test$p.value,
               n = n)
      } else {
        tibble(variable = .x, r = NA, p_value = NA, n = n)
      }
    })
}) |> 
  arrange(variable,desc(r)) |> 
  mutate(variable = ifelse(variable == 'GLDAS_GWS','GLDAS_GWSA',variable))

write_rds(cor_well,'data/processed/rds/downscaled_products_correlation.rds')

left_join(pozos, y = cor_mes_pozo) |> 
  writeVector('data/processed/vectorial/cor_GLDAS.gpkg',
              overwrite=T)

# visualizar

data_raw <- read_rds('data/processed/rds/downscaled_products_correlation.rds')

n_umbral <- data_raw |> # filtar el 75% con más datos y r significativos p-value < 0.05
  group_by(variable) |> 
  reframe(n_umbral = quantile(n,.25))

data <- data_raw |> 
  left_join(n_umbral) |> 
  filter(n > n_umbral,
         p_value < 0.05)

data_raw |> 
  filter(n > n_umbral,
         p_value < 0.05) |> 
  pull(codigo) |> 
  unique() |> 
  length()

data |>
  ggplot(aes(x = r)) +
  geom_histogram(binwidth = 0.2, boundary = 0, fill = 'darkcyan', color = 'black',linewidth = .2) +
  scale_x_continuous(limits = c(-1.1,1.1), expand = c(0,0), breaks = seq(-1, 1, by = 0.2)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,135), breaks = seq(0,140,by = 20)) +
  labs(x = 'Pearson r',y = 'n° of wells') +
  facet_wrap(~variable) +
  theme_bw() +
  theme(strip.background = element_rect(fill='white'))

ggsave('output/fig/validation/r_histogram.png',width = 8, height = 5)

data |>
  ggplot(aes(x = n)) +
  geom_histogram(binwidth = 10, boundary = 0, fill = 'darkcyan', color = 'black',linewidth = .2) +
  scale_x_continuous(expand = c(0.01,0), breaks = seq(0, 211, by = 20)) +
  labs(x = 'n° of monthly observation for correlation by well',y = 'frequency') +
  scale_y_continuous(limits = c(0,45), expand = c(0,0),
                     breaks = seq(0, 100, by = 10),
                     minor_breaks = seq(0, 150, by = 5)) +
  facet_wrap(~variable) +
  theme_bw()

ggsave('output/fig/validation/n_histogram.png',width = 8, height = 6)
