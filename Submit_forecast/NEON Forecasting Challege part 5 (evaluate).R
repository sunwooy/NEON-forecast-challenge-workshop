library(tidyverse)
library(lubridate)
library(slider)
library(scoringRules)
library(duckdbfs)
library(dplyr)
library(tidyr)
library(ggplot2)
# THIS IS START DATE (t=0 where we have observation)
reference_date <- "2025-03-24 00:00:00"

###### Steve's forecast ############
#------- Read data --------
# read in the targets data
targets <- read_csv("https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D/aquatics-targets.csv.gz")

# Filter the targets, add water temp lag and past 3 day air temperature mean
targets <- targets %>%
  filter(site_id == "BARC",
         variable == "temperature",
         datetime <= as_datetime(reference_date)) %>%
  mutate(wtemp_yday = lag(observation))



# Past stacked weather -----
weather_past_s3 <- neon4cast::noaa_stage3()
met_variables <- c("air_temperature")
weather_past <- weather_past_s3  |> 
  dplyr::filter(site_id == "BARC",
                datetime >= ymd('2017-01-01'),
                variable %in% met_variables,
                datetime <= as_datetime(reference_date)) |> 
  dplyr::collect()

# aggregate the past to mean values
weather_past_daily <- weather_past |> 
  mutate(datetime = as_date(datetime)) |> 
  group_by(datetime, site_id, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  # convert air temperature to Celsius if it is included in the weather data
  mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
  pivot_wider(names_from = variable, values_from = prediction)

targets_lm <- targets |> 
  pivot_wider(names_from = 'variable', values_from = 'observation') |> 
  left_join(weather_past_daily, by = c("datetime","site_id")) |>
  arrange(site_id, datetime) |>
  group_by(site_id) |>
  mutate(
    wtemp_yday = lag(temperature, 1),
    airtemp_yday = slide_dbl(
      air_temperature,
      ~ mean(.x, na.rm = FALSE), # Remove if any of 3 days  are missing (today and 2 days previous)
      .before = 2,
      .after = 0,
      .complete = TRUE
    )
  ) |>
  ungroup()|>
  drop_na()


# Future weather forecast --------
# New forecast only available at 5am UTC the next day
noaa_date <- as_datetime(reference_date) + days(1)

weather_future_s3 <- neon4cast::noaa_stage2(start_date = as.character(noaa_date))

weather_future <- weather_future_s3 |> 
  dplyr::filter(datetime >= reference_date,
                site_id=='BARC',
                variable %in% met_variables) |> 
  collect()

weather_future_daily <- weather_future |> 
  mutate(datetime = as_date(datetime)) |> 
  # mean daily forecasts at each site per ensemble
  group_by(datetime, site_id, parameter, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  # convert air temperature to Celsius if it is included in the weather data
  mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
  pivot_wider(names_from = variable, values_from = prediction) |> 
  select(any_of(c('datetime', 'site_id', met_variables, 'parameter')))



# Loop through each site to fit the model
forecast_df <- NULL

n_members = 31

# # Get horizon length by looking at available forecast dates that is not today. 

# forecasted_dates <- unique(weather_future_daily$datetime)

forecasted_dates <- as_datetime(reference_date) + 
  days(seq_along(unique(weather_future_daily$datetime)) - 1)

forecast_start_date <- forecasted_dates[1]




#Fit linear model based on past data: 
fit <- lm(targets_lm$temperature ~ targets_lm$air_temperature + targets_lm$wtemp_yday + targets_lm$airtemp_yday)
fit_summary <- summary(fit)

coeffs <- fit$coefficients
params_se <- fit_summary$coefficients[,2]

mod <- predict(fit, data=targets_lm$air_temperature)

r2 <- fit_summary$r.squared
residuals <- mod - targets_lm$temperature
err <- mean(residuals, na.rm = TRUE) 
rmse <- round(sqrt(mean((mod - targets_lm$temperature)^2, na.rm = TRUE)), 2) 



# Parameter Uncertanity
param_df <- data.frame(beta1 = rnorm(n_members, coeffs[1], params_se[1]),
                       beta2 = rnorm(n_members, coeffs[2], params_se[2]),
                       beta3 = rnorm(n_members, coeffs[3], params_se[3]),
                       beta4 = rnorm(n_members, coeffs[4], params_se[4]))

# Process Uncertainty
sigma <- sd(residuals, na.rm = TRUE)




future_air <- weather_future_daily |>
  arrange(parameter, datetime)

past_air <- weather_past_daily |>
  filter(
    datetime >= forecast_start_date - days(3),
    datetime < forecast_start_date
  ) |>
  select(datetime, site_id, air_temperature) |>
  crossing(future_air |> distinct(parameter)) |>
  arrange(parameter, datetime)

air_for_lags <- bind_rows(
  past_air,
  future_air |> select(datetime, site_id, air_temperature, parameter)
) |>
  arrange(parameter, datetime) |>
  group_by(parameter) |>
  mutate(
    airtemp_lag3 = slide_dbl(
      air_temperature,
      ~ mean(.x, na.rm = FALSE),
      .before = 2,
      .after = 0,
      .complete = TRUE
    )
  ) |>
  ungroup()





curr_wt <- targets_lm %>%
  filter(datetime == as_datetime(reference_date)) %>%
  pull(temperature)

air_init <- weather_future_daily %>%
  filter(datetime == forecast_start_date) %>%
  pull(air_temperature)

init_mean <- curr_wt[1]
ic_sd <- 0.1

ic_uc <- rnorm(n = n_members, mean = init_mean, sd = ic_sd)

ic_df <- tibble(
  forecast_date = rep(as.Date(reference_date), times = n_members),
  ensemble_member = 1:n_members,
  forecast_variable = "water_temperature",
  value = ic_uc
)

# Now we're ready to forecast
forecast_total_unc <- tibble(forecast_date = rep(forecasted_dates, times = n_members),
                             ensemble_member = rep(1:n_members, each = length(forecasted_dates)),
                             forecast_variable = "water_temperature",
                             value = as.double(NA),
                             uc_type = "total") %>%
  rows_update(ic_df, by = c("forecast_date","ensemble_member","forecast_variable")) 

for(j in 2:length(forecasted_dates)) {
  
  #pull prediction dataframe for relevant date
  temp_pred <- forecast_total_unc %>%
    filter(forecast_date == forecasted_dates[j])
  
  #pull driver ensemble for the relevant date; here we are using all 31 NOAA ensemble members
  air_lag <- air_for_lags %>%
    filter(datetime == forecasted_dates[j])
  
  #pull lagged water temp values
  temp_lag <- forecast_total_unc %>%
    filter(forecast_date == forecasted_dates[j-1])
  
  
  
  #run model using initial conditions and parameter distributions instead of fixed values, and adding process uncertainty. Initial condition uncertanity from forecast_total_unc where rows updated by ic_df
  temp_pred$value <- param_df$beta1 + air_lag$air_temperature * param_df$beta2 + 
    temp_lag$value * param_df$beta3 + air_lag$airtemp_lag3 * param_df$beta4 + rnorm(n = n_members, mean = 0, sd = sigma) 
  
  #insert values back into the forecast dataframe
  forecast_total_unc <- forecast_total_unc %>%
    rows_update(temp_pred, by = c("forecast_date","ensemble_member","forecast_variable"))
}


###### 11.10 Part 1: a single forecast

# 1, 2.  

all_results <- open_dataset("s3://bio230014-bucket01/challenges/forecasts/bundled-parquet/project_id=neon4cast/duration=P1D/variable=temperature", 
                            s3_endpoint = "sdsc.osn.xsede.org", 
                            anonymous = TRUE)

my_forecasts <- all_results |> 
  filter(reference_datetime == as_datetime(reference_date),
         site_id == "BARC",
         model_id %in% c("climatology", "persistenceRW")) |> 
  collect()


clim <- my_forecasts |> 
  filter(model_id == "climatology") |> 
  pivot_wider(names_from = parameter, values_from = prediction)

# Plot an ensemble forecast
persis <- my_forecasts |> 
  filter(model_id == "persistenceRW") 


forecast_total_unc <- forecast_total_unc %>%
  rename('datetime' = 'forecast_date',
         'parameter' = 'ensemble_member',
         'prediction' = 'value')


# 3. 

obs <- read_csv("https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D/aquatics-targets.csv.gz")

obs <- obs %>%
  filter(
    site_id == "BARC",
    variable == "temperature",
    datetime %in% unique(forecast_total_unc$datetime)) %>%
  mutate(
    model = "observation",
    center=observation
  )

  
clim_plot <- my_forecasts |> 
  filter(model_id == "climatology") |> 
  pivot_wider(names_from = parameter, values_from = prediction) |> 
  mutate(
    model = "climatology",
    lower = mu - 1.96 * sigma,
    upper = mu + 1.96 * sigma,
    center = mu
  ) |> 
  select(datetime, model, lower, upper, center)

persis_plot <- my_forecasts |> 
  filter(model_id == "persistenceRW") |> 
  group_by(datetime) |> 
  summarise(
    model = "persistenceRW",
    lower = quantile(prediction, 0.05, na.rm = TRUE),
    upper = quantile(prediction, 0.95, na.rm = TRUE),
    center = mean(prediction, na.rm = TRUE),
    .groups = "drop"
  )


total_unc_plot <- forecast_total_unc |> 
  group_by(datetime) |> 
  summarise(
    model = "forecast_total_unc",
    lower = quantile(prediction, 0.05, na.rm = TRUE),
    upper = quantile(prediction, 0.95, na.rm = TRUE),
    center = mean(prediction, na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- bind_rows(clim_plot, persis_plot, total_unc_plot, obs)

ggplot(plot_df, aes(x = datetime, ymin = lower, ymax = upper, fill = model)) +
  geom_ribbon(alpha = 0.25) +
  geom_line(aes(y = center, color = model), linewidth = 1)



#4. Based on visual inspection, the persistence model is the worst compared to the observation. Additionally, while my model does a great job with prediction in the beginning,
# it tends to fail further out in time relative to climatology. 


#5. Based on visual inspection, the persistence model stands out as a model with the most uncertainty with little changes in the mean temperature prediction across forecast, suggesting it's
# not a great model for water temperature. This makes sense since today's water temp isn't likely related to tmrw's water temperature. When comparing my model with climatology, I see that 
# although I start with much lower uncertainty in the beginning, uncertainty increases with increasing horizon relative to climatology. Perhaps most concerning is that 
# all three models fail to capture peak temperature observed on April 07th in the uncertanity range. 



# 6, 7. 
my_forecasts <- all_results |> 
  filter(reference_datetime == as_datetime(reference_date),
         site_id == "BARC",
         model_id %in% c("climatology", "persistenceRW")) |> 
  collect()



multi_forecast_ensemble <- my_forecasts |> 
  filter(model_id == "persistenceRW") |> 
  left_join(obs, by = c("datetime", "site_id", "variable")) |> 
  select(model_id, datetime, prediction, observation, parameter)

multi_forecast_parametric <- my_forecasts |> 
  filter(model_id == "climatology") |> 
  left_join(obs, by = c("datetime", "site_id", "variable")) |> 
  select(model_id, datetime, prediction, observation, parameter) |> 
  pivot_wider(names_from = parameter, values_from = prediction)

multi_forecast_ensemble2 <- forecast_total_unc |> 
  mutate(model_id = "forecast_total_unc") |> 
  left_join(obs, by = c("datetime")) |> 
  select(model_id, datetime, prediction, observation, parameter)



normal_crps <- multi_forecast_parametric |> 
  group_by(datetime, model_id) |> 
  summarize(crps = crps_norm(observation, mean = mu, sd = sigma))


# Scoring multiple forecaasts using crps_sample in a group_by summarize required
# making a small function to get a single observation
crps_tidy_ensemble <- function(prediction, observation) {
  obs <- observation[1] #all the observations are the same for each row, but we only need one
  crps_sample(obs, prediction)
}

#Use the function
ensemble_crps <- multi_forecast_ensemble |> 
  group_by(datetime, model_id) |> 
  summarize(crps = crps_tidy_ensemble(prediction, observation)) 

#Use the function
ensemble_crps2 <- multi_forecast_ensemble2 |> 
  group_by(datetime, model_id) |> 
  summarize(crps = crps_tidy_ensemble(prediction, observation)) |>
  filter(datetime > as_datetime(reference_date))

#Plot the histogram of the scores

bind_rows(normal_crps, ensemble_crps, ensemble_crps2) |> 
  ggplot(aes(x = crps, fill = model_id)) +
  geom_histogram()


crps_all <- bind_rows(normal_crps, ensemble_crps, ensemble_crps2) |>
  ungroup() |>
  mutate(model_id = factor(
    model_id,
    levels = c("climatology", "persistenceRW", "forecast_total_unc")
  ))

ggplot(crps_all, aes(x = datetime, y = crps, color = model_id)) +
  geom_line() +
  geom_point() 


#  Similar to Question 4 and 5, the climatology and our model performs the best. Again, while my model seems to perform better in the beginning, it  is outperformed by climatology  later.






###### 11.10 Part 2: multiple forecast

all_results <- open_dataset("s3://bio230014-bucket01/challenges/forecasts/bundled-parquet/project_id=neon4cast/duration=P1D/variable=temperature", 
                            s3_endpoint = "sdsc.osn.xsede.org", 
                            anonymous = TRUE)
my_forecasts <- all_results |> 
  filter(
    site_id == "BARC",
    model_id %in% c("climatology", "persistenceRW")) |> 
  collect()

forecast_dates = unique(my_forecasts$datetime)

obs <- read_csv("https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D/aquatics-targets.csv.gz")


multi_forecast_ensemble <- my_forecasts |> 
  filter(model_id == "persistenceRW") |> 
  left_join(obs, by = c("datetime", "site_id", "variable")) |> 
  mutate(horizon = as.numeric(datetime - reference_datetime, units="days"))|>
  select(model_id, datetime, prediction, observation, parameter, horizon) |>
  arrange(datetime)

multi_forecast_parametric <- my_forecasts |> 
  filter(model_id == "climatology") |>
  left_join(obs, by = c("datetime", "site_id", "variable")) |> 
  mutate(horizon = as.numeric(datetime - reference_datetime, units="days"))|>
  select(model_id, datetime, prediction, observation, parameter, horizon) |>
  arrange(datetime)|>
  pivot_wider(
    names_from = parameter,
    values_from = prediction,
    values_fn = mean
  )



normal_crps <- multi_forecast_parametric |> 
  group_by(datetime, model_id, horizon) |> 
  summarize(crps = crps_norm(observation, mean = mu, sd = sigma)) |>
  filter(!is.na(crps))



crps_tidy_ensemble <- function(prediction, observation) {
  obs <- observation[1] #all the observations are the same for each row, but we only need one
  crps_sample(obs, prediction)
}

#Use the function
ensemble_crps <- multi_forecast_ensemble |> 
  filter(!is.na(observation), !is.na(prediction)) |> 
  group_by(datetime, model_id, horizon) |> 
  summarize(crps = crps_tidy_ensemble(prediction, observation)) 


mean(ensemble_crps$crps)
mean(normal_crps$crps)



# Question 9. 

# Climatology has the lowest mean crps with 0.88083C


# Question 10.

# No, mean CRPS is not enough to determine best model because certain models might perform better at certain horizons. We would want CRPS at horizon intervals. 




