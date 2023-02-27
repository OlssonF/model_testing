# Script to test Simstrat model parameters
# Starts by reading in the JSON file 
# sets the parameter values to different combinations
# writes the JSON and then runs Simstrat with new parameter values
# read in the output - assess fit using RMSE and NSE
# assign to output table
library(tidyverse)
lake <- 'CRAM'
model <- 'Simstrat'
dir <- here::here()

working_dir <- file.path(dir, lake, model)

setwd(working_dir)

# depths to model
depth_v <- 1
# Observations to test model fit
max_depth <- read.table('hypsograph.dat', header = T) |> 
  summarise(max_depth = max(-Depth..m.)) |> 
  mutate(max_depth = round(max_depth/depth_v)*depth_v) |> 
  pull()

cuts <- tibble::tibble(depth = seq(0, max_depth, depth_v),
                       cuts = as.integer(factor(depth)))

obs <- readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-expanded-observations.csv.gz",
                show_col_types = FALSE) |> 
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
f_wind <- round(seq(0.5, 1.5, length.out = 12), 2)
p_lw <- round(seq(0.5, 1.5, length.out = 12), 2)

params <- expand.grid(f_wind = f_wind, p_lw = p_lw)

# for model output
timestep <- 300
reference_year <- 2022

# table to fill
params <- params |> 
  mutate(RMSE = NA,
         NSE = NA)

for (i in 1:nrow(params)) {
  # Read in json
  LakeEnsemblR::input_json(file = 'simstrat.par', label = 'ModelParameters', key = 'f_wind', value = params$f_wind[i])
  LakeEnsemblR::input_json(file = 'simstrat.par', label = 'ModelParameters', key = 'p_lw', value = params$p_lw[i]) 
  # Run simstrat with updated parameters
  SimstratR::run_simstrat()
  
  # read in the output
  temp <- read.table(file.path("output", "T_out.dat"), header = TRUE, sep = ",", 
                     check.names = FALSE)
  temp[, 1] <- as.POSIXct(temp[, 1] * 3600 * 24, 
                          origin = paste0(reference_year, "-01-01"))
  temp[, 1] <- lubridate::round_date(temp[, 1], 
                                     unit = lubridate::seconds_to_period(timestep))
  temp <- temp[, c(1, ncol(temp):2)]
  temp <- temp[, colSums(is.na(temp)) < nrow(temp)]
  
  colnames(temp) <- c('datetime', -as.numeric(colnames(temp)[-1]))
  
  # Plotting
  # temp |>  pivot_longer(cols = -datetime, names_to = 'depth', values_to = 'temp') |> 
  #   ggplot(aes(x=Datetime, y=temp, colour = depth)) +
  #   geom_line() + 
  #   geom_point
  # 
  # # Comparison with observations
  # temp |> 
  #   pivot_longer(cols = -datetime,
  #                names_to = 'depth',
  #                values_to = 'prediction') |>
  #   mutate(depth = as.numeric(depth)) |> 
  #   inner_join(obs, by = c('depth', 'datetime')) |> 
  #   ggplot(aes(x=datetime)) +
  #   geom_point(aes(y=observation), alpha = 0.2) +
  #   geom_line(aes(y=prediction)) +
  #   facet_wrap(~depth)
  
  
  obs_pred_matrix <- temp |> 
    pivot_longer(cols = -datetime,
                 names_to = 'depth',
                 values_to = 'prediction') |>
    mutate(depth = as.numeric(depth)) |> 
    inner_join(obs, by = c('depth', 'datetime')) 
  
  # calculate the metrics of goodness of fit
  params$RMSE[i] <- hydroGOF::rmse(sim = obs_pred_matrix$prediction, obs = obs_pred_matrix$observation)
  params$NSE[i] <- hydroGOF::NSE(sim = obs_pred_matrix$prediction, obs = obs_pred_matrix$observation)

  message('Parameter set ', i, '/',nrow(params), ' at ', Sys.time(), 
          ' Simstrat fit with parameters f_wind = ', params$f_wind[i], ' and p_lw = ', params$p_lw[i])
}

write_csv(params, file = paste0(lake, '_params.csv'))

RMSE_p <- params |> 
  mutate(best_NSE = if_else(NSE == max(NSE), T, NA),
         best_RMSE = if_else(RMSE == min(RMSE), T, NA)) |> 
  ggplot(aes(x=f_wind, y=p_lw, fill = RMSE)) +  
  geom_tile() +
  geom_point(aes(colour=best_RMSE)) +
  scale_colour_manual(breaks = c(T), values=c('red'), na.value = NA, guide="none")

NSE_p <- params |> 
  mutate(best_NSE = if_else(NSE == max(NSE), T, NA),
         best_RMSE = if_else(RMSE == min(RMSE), T, NA),
         NSE = ifelse(NSE <0, NA, NSE)) |> 
  ggplot(aes(x=f_wind, y=p_lw, fill = NSE)) +  
  geom_tile() +
  geom_point(aes(colour=best_NSE)) +
  scale_fill_continuous(breaks = 'reverse', na.value = NA ) +
  scale_colour_manual(breaks = c(T), values=c('red'), na.value = NA, guide="none")


# extract best params and run again to plot obs v pred
best_fwind <-  params %>% slice_min(RMSE) |> select(f_wind) |> pull()
# best_fwind <-  1.4
best_plw <-  params %>% slice_min(RMSE) |> select(p_lw) |> pull()
# best_plw <-  0.95


# Read in j.son
LakeEnsemblR::input_json(file = 'simstrat.par', label = 'ModelParameters', key = 'f_wind', value = best_fwind)
LakeEnsemblR::input_json(file = 'simstrat.par', label = 'ModelParameters', key = 'p_lw', value = best_plw) 
# Run simstrat with updated parameters
SimstratR::run_simstrat()

# read in the output
temp <- read.table(file.path("output", "T_out.dat"), header = TRUE, sep = ",", 
                   check.names = FALSE)
temp[, 1] <- as.POSIXct(temp[, 1] * 3600 * 24, 
                        origin = paste0(reference_year, "-01-01"))
temp[, 1] <- lubridate::round_date(temp[, 1], 
                                   unit = lubridate::seconds_to_period(timestep))
temp <- temp[, c(1, ncol(temp):2)]
temp <- temp[, colSums(is.na(temp)) < nrow(temp)]

colnames(temp) <- c('datetime', -as.numeric(colnames(temp)[-1]))

# Plotting
temp |>  pivot_longer(cols = -datetime, names_to = 'depth', values_to = 'temp') |>
  ggplot(aes(x=datetime, y=temp, colour = depth)) +
  geom_line() 


# Comparison with observations
temp |>
  pivot_longer(cols = -datetime,
               names_to = 'depth',
               values_to = 'prediction') |>
  mutate(depth = as.numeric(depth)) |>
  right_join(obs, by = c('depth', 'datetime')) |>
  filter(depth %in% c(0,1,2,3,4,5,6,8,10)) |> 
  ggplot(aes(x=datetime)) +
  geom_point(aes(y=observation), alpha = 0.2) +
  geom_line(aes(y=prediction)) +
  facet_wrap(~depth) +
  coord_cartesian(xlim = as_datetime(c('2022-01-01', '2022-12-31')))


obs_pred_matrix <- temp |> 
  pivot_longer(cols = -datetime,
               names_to = 'depth',
               values_to = 'prediction') |>
  mutate(depth = as.numeric(depth)) |> 
  inner_join(obs, by = c('depth', 'datetime')) 

# calculate the metrics of goodness of fit
hydroGOF::rmse(sim = obs_pred_matrix$prediction, obs = obs_pred_matrix$observation)
hydroGOF::NSE(sim = obs_pred_matrix$prediction, obs = obs_pred_matrix$observation)

obs_pred_matrix %>%
  group_by(depth) |> 
  summarise(rmse = hydroGOF::rmse(prediction, observation))
