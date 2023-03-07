# Script to test GOTM model parameters
# Starts by reading in the yml file 
# sets the parameter values to different combinations
# writes the yml and then runs GOTM with new parameter values
# read in the output - assess fit using RMSE and NSE
# assign to output table
library(tidyverse)
lake <- 'CRAM'
model <- 'GOTM'
dir <- here::here()

working_dir <- file.path(dir, lake, model)

setwd(working_dir)

# depths to model
depth_n <- 19
# Observations to test model fit
max_depth <- read.table('hypsograph.dat', skip = 1, col.names = c('depth', 'area')) |> 
  summarise(max_depth = max(-depth)) |> 
  pull()

cuts <- tibble::tibble(depth = round(seq(0, max_depth, length.out = depth_n), 2),
                       cuts = as.integer(factor(depth)))

obs <- readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-expanded-observations.csv.gz",
                       show_col_types = FALSE) |> 
  mutate(site_id = ifelse(site_id == 'TOOK', 'TOOL', site_id)) |> 
  filter(site_id == lake) |> 
  dplyr::mutate(cuts = cut(depth, breaks = cuts$depth, 
                           include.lowest = TRUE, right = FALSE, labels = FALSE)) |>
  dplyr::filter(lubridate::hour(datetime) == 0) |>
  dplyr::group_by(cuts, variable, datetime, site_id) |>
  dplyr::summarize(observation = mean(observation, na.rm = TRUE), .groups = "drop") |>
  dplyr::left_join(cuts, by = "cuts") |>
  dplyr::select(site_id, datetime, variable, depth, observation) |> 
  na.omit()

# Parameters to test
shf <- round(seq(0.5, 1.5, length.out = 5), 2)
wsf <- round(seq(0.5, 1.5, length.out = 5), 2)
swr <- round(seq(0.5, 1.5, length.out = 5), 2)
params <- expand.grid(shf = shf, wsf = wsf, swr = swr)


# table to fill
params <- params |> 
  mutate(RMSE = NA,
         NSE = NA)

for (i in 1:nrow(params)) {
  # Read in yaml
  LakeEnsemblR::input_yaml_multiple(file = 'gotm.yaml',key1 = 'surface', key2 = 'meteo', key3 = 'swr', key4 = 'scale_factor', 
                                    value = params$swr[i])
  LakeEnsemblR::input_yaml_multiple(file = 'gotm.yaml',key1 = 'surface', key2 = 'fluxes', key3 = 'heat', key4 = 'scale_factor', 
                                    value = params$shf[i])
  LakeEnsemblR::input_yaml_multiple(file = 'gotm.yaml',key1 = 'surface', key2 = 'meteo', key3 = 'u10', key4 = 'scale_factor', 
                                    value = params$wsf[i])
  # Run simstrat with updated parameters
  GOTMr::run_gotm()
  
  # read in the output
  nc <- ncdf4::nc_open(file.path("output", "output.nc"))
  tim = ncdf4::ncvar_get(nc, "time")
  tunits = ncdf4::ncatt_get(nc, "time")
  lnam = tunits$long_name
  tustr <- strsplit(tunits$units, " ")
  step = tustr[[1]][1]
  tdstr <- strsplit(unlist(tustr)[3], "-")
  tmonth <- as.integer(unlist(tdstr)[2])
  tday <- as.integer(unlist(tdstr)[3])
  tyear <- as.integer(unlist(tdstr)[1])
  tdstr <- strsplit(unlist(tustr)[4], ":")
  thour <- as.integer(unlist(tdstr)[1])
  tmin <- as.integer(unlist(tdstr)[2])
  origin <- as.POSIXct(paste0(tyear, "-", tmonth,
                              "-", tday, " ", thour, ":", tmin),
                       format = "%Y-%m-%d %H:%M", tz = "UTC")
  if (step == "hours") {
    tim <- tim * 60 * 60
  }
  if (step == "minutes") {
    tim <- tim * 60
  }
  
  time <- as.POSIXct(tim, origin = origin, tz = "UTC")
  
  
  temp <- ncdf4::ncvar_get(nc, "temp")
  z <- ncdf4::ncvar_get(nc, "z")
  
  gotm_output <- as.data.frame(t(temp))
  gotm_output$datetime <- time
  gotm_output <- gotm_output[, c(ncol(gotm_output), 1:(ncol(gotm_output) - 1))]
  
  colnames(gotm_output) <- c("datetime", rev(cuts$depth))
  
  obs_pred_matrix <- gotm_output |> 
    pivot_longer(cols = -datetime,
                 names_to = 'depth',
                 values_to = 'prediction') |>
    mutate(depth = round(as.numeric(depth),2)) |> 
    inner_join(obs, by = c('depth', 'datetime')) 
  
  # calculate the metrics of goodness of fit
  params$RMSE[i] <- hydroGOF::rmse(sim = obs_pred_matrix$prediction, obs = obs_pred_matrix$observation)
  params$NSE[i] <- hydroGOF::NSE(sim = obs_pred_matrix$prediction, obs = obs_pred_matrix$observation)
  
  message('Parameter set ', i, '/',nrow(params), ' at ', Sys.time(), 
          ' GOTM run with parameters wsf = ', params$wsf[i], ',',
          ' shf = ', params$shf[i], ',',
          ' and swr = ', params$swr[i])
}

# plot parameter fits
RMSE_p <- params |> 
  mutate(best_NSE = if_else(NSE == max(NSE), T, NA),
         best_RMSE = if_else(RMSE == min(RMSE), T, NA)) |> 
  ggplot(aes(x=wsf, y=shf, fill = RMSE)) +  
  geom_tile() +
  geom_point(aes(colour=best_RMSE)) +
  scale_fill_continuous(guide = 'none', na.value = NA, breaks = 'reverse') +
  scale_colour_manual(breaks = c(T), values=c('red'), na.value = NA, guide="none") +
  facet_wrap(~swr, nrow = 1, labeller = label_both)

NSE_p <- params |> 
  mutate(best_NSE = if_else(NSE == max(NSE), T, NA),
         best_RMSE = if_else(RMSE == min(RMSE), T, NA),
         NSE = ifelse(NSE <0, NA, NSE)) |> 
  ggplot(aes(x=wsf, y=shf, fill = NSE)) +  
  geom_tile() +
  geom_point(aes(colour=best_NSE)) +
  scale_fill_continuous(breaks = 'reverse', na.value = NA ) +
  scale_colour_manual(breaks = c(T), values=c('red'), na.value = NA, guide="none")+
  facet_wrap(~swr, nrow = 1, labeller = label_both)

ggpubr::ggarrange(RMSE_p, NSE_p, ncol = 1)

# write parameter fits
write_csv(params, file = paste0(lake, '_params.csv'))

# extract best params and run again to plot obs v pred
best_wsf <-  params %>% slice_min(RMSE) |> select(wsf) |> pull()
best_shf <-  params %>% slice_min(RMSE) |> select(shf) |> pull()
best_swr <-  params %>% slice_min(RMSE) |> select(swr) |> pull()


# Read in yaml
LakeEnsemblR::input_yaml_multiple(file = 'gotm.yaml',key1 = 'surface', key2 = 'meteo', key3 = 'swr', key4 = 'scale_factor', 
                                  value = best_swr)
LakeEnsemblR::input_yaml_multiple(file = 'gotm.yaml',key1 = 'surface', key2 = 'fluxes', key3 = 'heat', key4 = 'scale_factor', 
                                  value = best_shf)
LakeEnsemblR::input_yaml_multiple(file = 'gotm.yaml',key1 = 'surface', key2 = 'meteo', key3 = 'u10', key4 = 'scale_factor', 
                                  value = best_wsf)
# Run gotm with updated parameters
GOTMr::run_gotm()

# read in the output
nc <- ncdf4::nc_open(file.path("output", "output.nc"))
tim = ncdf4::ncvar_get(nc, "time")
tunits = ncdf4::ncatt_get(nc, "time")
lnam = tunits$long_name
tustr <- strsplit(tunits$units, " ")
step = tustr[[1]][1]
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth <- as.integer(unlist(tdstr)[2])
tday <- as.integer(unlist(tdstr)[3])
tyear <- as.integer(unlist(tdstr)[1])
tdstr <- strsplit(unlist(tustr)[4], ":")
thour <- as.integer(unlist(tdstr)[1])
tmin <- as.integer(unlist(tdstr)[2])
origin <- as.POSIXct(paste0(tyear, "-", tmonth,
                            "-", tday, " ", thour, ":", tmin),
                     format = "%Y-%m-%d %H:%M", tz = "UTC")
if (step == "hours") {
  tim <- tim * 60 * 60
}
if (step == "minutes") {
  tim <- tim * 60
}

time <- as.POSIXct(tim, origin = origin, tz = "UTC")


temp <- ncdf4::ncvar_get(nc, "temp")
z <- ncdf4::ncvar_get(nc, "z")

gotm_output <- as.data.frame(t(temp))
gotm_output$datetime <- time
gotm_output <- gotm_output[, c(ncol(gotm_output), 1:(ncol(gotm_output) - 1))]

colnames(gotm_output) <- c("datetime", rev(cuts$depth))

obs_pred_matrix <- gotm_output |> 
  pivot_longer(cols = -datetime,
               names_to = 'depth',
               values_to = 'prediction') |>
  mutate(depth = round(as.numeric(depth),2)) |> 
  inner_join(obs, by = c('depth', 'datetime')) 

# calculate the metrics of goodness of fit
hydroGOF::rmse(sim = obs_pred_matrix$prediction, obs = obs_pred_matrix$observation)
hydroGOF::NSE(sim = obs_pred_matrix$prediction, obs = obs_pred_matrix$observation)

# Comparison with observations
gotm_output |>
  pivot_longer(cols = -datetime,
               names_to = 'depth',
               values_to = 'prediction')|>
  mutate(depth = round(as.numeric(depth),2)) |> 
  inner_join(obs, by = c('depth', 'datetime')) |>
  ggplot(aes(x=datetime)) +
  geom_point(aes(y=observation), alpha = 0.2) +
  geom_line(aes(y=prediction)) +
  facet_wrap(~depth) +
  coord_cartesian(xlim = lubridate::as_datetime(c('2022-01-01', '2022-12-31'))) +
  labs(title = paste('best_wsf=', best_wsf, 'best_shf=', best_shf,'best_swr=', best_swr))


obs_pred_matrix %>%
  group_by(depth) |> 
  summarise(rmse = hydroGOF::rmse(prediction, observation))

