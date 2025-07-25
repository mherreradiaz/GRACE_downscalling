library(terra)
library(tidyterra)
library(tidyverse)
library(glue)

data_GWD <- read_rds('data/processed/rds/GWD_chile.rds')
pozos <- vect('data/processed/vectorial/pozos_chile.gpkg')

grace_files <- dir_ls('data/processed/raster/GRACE_highres/')
fechas <- str_extract(grace_files,'(?<=res_).*?(?=\\.tif)')
grace_r <- grace_files |> 
  rast() |> 
  setNames(fechas)

data_grace <- terra::extract(grace_r,pozos) |> 
  mutate(codigo = pozos$codigo,
         .before = ID) |> 
  select(-ID) |> 
  pivot_longer(c(everything(),-codigo), names_to = 'fecha',
               values_to = 'lwe') |> 
  mutate(fecha = as.Date(paste0(fecha,'-01'))) |> 
  select(fecha,codigo,'lwe')

data <- data_grace |> 
  left_join(data_GWD) |> 
  group_by(cuenca,codigo) |> 
  mutate(GWD = as.numeric(scale(GWD))) |> 
  ungroup() |> 
  select(fecha,cuenca,codigo,GWD,lwe)

cor_general <- data |> # correlación general
  reframe(r = cor(GWD, lwe, method = 'pearson', use = 'complete.obs'))

cor_mes_pozo <- data |> # correlacion temporal a nivel de pozo
  group_by(codigo) |> 
  reframe({
    n <- sum(complete.cases(GWD, lwe))
    if (n > 3) {
      test <- cor.test(GWD, lwe, method = 'pearson', use = 'complete.obs')
      tibble(r = test$estimate,
             p_value = test$p.value,
             n = n)
    } else {
      tibble(r = NA, p_value = NA, n = n)
    }
  }) |> 
  na.omit()

cor_mes_pozo |> 
  arrange(desc(r))

write_rds(cor_mes_pozo,'data/processed/rds/cor_GRACE_hr_mes_pozo.rds')

left_join(pozos, y = cor_mes_pozo) |> 
  writeVector('data/processed/vectorial/pozos/cor_GRACE_hr.shp',
              overwrite=T)

cor_mes <- data |> # correlacion de cada mes
  group_by(mes = month(fecha)) |> 
  reframe({
    n <- sum(complete.cases(GWD, lwe))
    if (n > 3) {
      test <- cor.test(GWD, lwe, method = 'pearson', use = 'complete.obs')
      tibble(r = test$estimate,
             p_value = test$p.value,
             n = n)
    } else {
      tibble(r = NA, p_value = NA, n = n)
    }
  }) |> 
  na.omit()

write_rds(cor_mes,'data/processed/rds/cor_GRACE_hr_mes.rds')

# visualizar

data_raw <- read_rds('data/processed/rds/cor_GRACE_hr_mes_pozo.rds')

n_umbral <- data_raw |> # filtar el 75% con más datos y r significativos p-value < 0.05
  pull(n) |> 
  quantile(.25)

data <- data_raw |> 
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
  geom_histogram(binwidth = 0.2, boundary = 0, fill = 'steelblue', color = 'black',linewidth = .2) +
  scale_x_continuous(expand = c(0.01,0), breaks = seq(-1, 1, by = 0.2)) +
  labs(x = 'Pearson r',y = 'frequency') +
  scale_y_continuous(limits = c(0,150), expand = c(0,0),
                     breaks = seq(0, 140, by = 20),
                     minor_breaks = seq(0, 150, by = 10)) +
  theme_bw()

ggsave('output/fig/correlation/GRACE_highres_GWD/r_histogram.png',width = 8, height = 6)

data |>
  ggplot(aes(x = n)) +
  geom_histogram(binwidth = 10, boundary = 0, fill = 'steelblue', color = 'black',linewidth = .2) +
  scale_x_continuous(expand = c(0.01,0), breaks = seq(0, 211, by = 20)) +
  labs(x = 'n° of monthly observation for correlation by well',y = 'frequency') +
  scale_y_continuous(limits = c(0,100), expand = c(0,0),
                     breaks = seq(0, 100, by = 40),
                     minor_breaks = seq(0, 150, by = 10)) +
  theme_bw()

ggsave('output/fig/correlation/GRACE_highres_GWD/n_histogram.png',width = 8, height = 6)

data |> # porcentaje de r
  mutate(grupo_r = cut(r, breaks = seq(-1, 1, by = 0.2), include.lowest = TRUE)) |> 
  group_by(grupo_r) |> 
  reframe(n = n()/nrow(data)*100)


