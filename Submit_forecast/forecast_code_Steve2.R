## install.packages('remotes')
## install.packages('tidyverse') # collection of R packages for data manipulation, analysis, and visualisation
## install.packages('lubridate') # working with dates and times
## remotes::install_github('eco4cast/neon4cast') # package from NEON4cast challenge organisers to assist with forecast building and submission

# ------ Load packages -----
library(tidyverse)
library(lubridate)
library(slider)

#--------------------------#

# Change this for your model ID
# Include the word "example" in my_model_id for a test submission
# Don't include the word "example" in my_model_id for a forecast that you have registered (see neon4cast.org for the registration form)
my_model_id <- 'example_model'

# --Model description--- #

# Add a brief description of your modeling approach

# It's a linear model that uses air temperature, a 3-day mean air temperature, and a 1-day lag water temperature as predictor for water temperature. 

# -- Uncertainty representation -- #

# Describe what sources of uncertainty are included in your forecast and how you estimate each source.

# Driver, parameter, and process uncertanity are included in the model. In addition, because initial water temperature data are not available at forecast start date, we initiated 
# water temperature as air temperature with uncertainty of 2C as initial conditions (normal distribution). 

#------- Read data --------
# read in the targets data
targets <- read_csv("https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D/aquatics-targets.csv.gz")

# read in the sites data
aquatic_sites <- read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-ci/refs/heads/main/neon4cast_field_site_metadata.csv") |>
  dplyr::filter(aquatics == 1)

focal_sites <- aquatic_sites |> 
  filter(field_site_subtype == 'Lake') |> 
  pull(field_site_id)

# Filter the targets
# Date up to 2026-03-31
targets <- targets %>%
  filter(site_id %in% focal_sites,
         variable == 'temperature')

targets %>%
  group_by(site_id) %>%
  summarise(max_val = max(datetime, na.rm = TRUE))
#--------------------------#



# ------ Weather data ------
met_variables <- c("air_temperature")

# Past stacked weather -----
weather_past_s3 <- neon4cast::noaa_stage3()

weather_past <- weather_past_s3  |> 
  dplyr::filter(site_id %in% focal_sites,
                datetime >= ymd('2017-01-01'),
                variable %in% met_variables) |> 
  dplyr::collect()

# aggregate the past to mean values
weather_past_daily <- weather_past |> 
  mutate(datetime = as_date(datetime)) |> 
  group_by(datetime, site_id, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  # convert air temperature to Celsius if it is included in the weather data
  mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
  pivot_wider(names_from = variable, values_from = prediction)


  

# Future weather forecast --------
# New forecast only available at 5am UTC the next day
forecast_date <- Sys.Date() 
noaa_date <- forecast_date - days(1)

weather_future_s3 <- neon4cast::noaa_stage2(start_date = as.character(noaa_date))

weather_future <- weather_future_s3 |> 
  dplyr::filter(datetime >= forecast_date,
                site_id %in% focal_sites,
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

#--------------------------#




# ----- Fit model & generate forecast----

# Generate a dataframe to fit the model to 
targets_lm <- targets |> 
  pivot_wider(names_from = 'variable', values_from = 'observation') |> 
  left_join(weather_past_daily, by = c("datetime","site_id")) |>
  arrange(site_id, datetime) |>
  group_by(site_id) |>
  mutate(
    wtemp_yday = lag(temperature, 1),
    airtemp_yday = slide_dbl(
      air_temperature,
      ~ mean(.x, na.rm = TRUE),
      .before = 3,
      .after = -1,
      .complete = TRUE
    )
  ) |>
  ungroup()|>
  drop_na()



# Loop through each site to fit the model
forecast_df <- NULL


n_members = 31

# # Get horizon length by looking at available forecast dates that is not today. 

forecast_dates_avail <- unique(weather_future_daily$datetime)
forecast_start_date <- forecast_dates_avail[1]
forecasted_dates <- forecast_dates_avail[1:length(forecast_dates_avail)]

for(i in 1:length(focal_sites)) {  
  
  curr_site <- focal_sites[i]
  
  site_target <- targets_lm |>
    filter(site_id == curr_site)

  
  #Fit linear model based on past data: 
  fit <- lm(site_target$temperature ~ site_target$air_temperature + site_target$wtemp_yday + site_target$airtemp_yday)
  fit_summary <- summary(fit)
  
  coeffs <- fit$coefficients
  params_se <- fit_summary$coefficients[,2]
  
  mod <- predict(fit, data=site_target$air_temperature)
  
  r2 <- fit_summary$r.squared
  residuals <- mod - site_target$temperature
  err <- mean(residuals, na.rm = TRUE) 
  rmse <- round(sqrt(mean((mod - site_target$temperature)^2, na.rm = TRUE)), 2) 
  
  
  
  # Parameter Uncertanity
  param_df <- data.frame(beta1 = rnorm(n_members, coeffs[1], params_se[1]),
                         beta2 = rnorm(n_members, coeffs[2], params_se[2]),
                         beta3 = rnorm(n_members, coeffs[3], params_se[3]),
                         beta4 = rnorm(n_members, coeffs[4], params_se[4]))
  
  # Process Uncertainty
  sigma <- sd(residuals, na.rm = TRUE)
  
  
  # Since we have a 3 day lag air temperature, we are going to get past 3 day air temp from the forecast_start_date (2 days from forecast_date)
  future_air <- weather_future_daily |>
    filter(site_id == curr_site) |>
    arrange(parameter, datetime)
  
  past_air <- weather_past_daily |>
    filter(
      site_id == curr_site,
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
        ~ mean(.x, na.rm = TRUE),
        .before = 3,
        .after = -1,
        .complete = TRUE
      )
    ) |>
    ungroup()
  
  
  # Most of our data is not up-to-date. Thus, we will have initial conditions as air temperature with 2C uncertanity.
  
  curr_wt <- weather_past_daily %>%
    filter(site_id == curr_site,
           datetime == forecast_start_date) %>%
    pull(air_temperature)
    
  ic_sd <- 2
  ic_uc <- rnorm(n = n_members, mean = curr_wt, sd = ic_sd)
  ic_df <- tibble(forecast_date = rep(as.Date(forecast_start_date), times = n_members),
                  ensemble_member = c(1:n_members),
                  forecast_variable = "water_temperature",
                  value = ic_uc,
                  uc_type = "initial_conditions")
  
  
  
  
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
  
  

  # # put all the relevant information into a tibble that we can bind together
  # curr_site_df <- tibble(datetime = noaa_future_site$datetime,
  #                        site_id = curr_site,
  #                        parameter = noaa_future_site$parameter,
  #                        prediction = forecasted_temperature,
  #                        variable = "temperature") #Change this if you are forecasting a different variable
  curr_site_df <- forecast_total_unc %>%
    transmute(
      datetime = forecast_date,
      site_id = curr_site,
      parameter = ensemble_member,
      prediction = value,
      variable = "temperature",
      uc_type = uc_type
    )
  
  forecast_df <- bind_rows(forecast_df, curr_site_df)
  message(curr_site, 'forecast run')
  
  
  
}


#--------------------------#


#---- Covert to EFI standard ----

# Make forecast fit the EFI standards
forecast_df_EFI <- forecast_df %>%
  filter(datetime > forecast_date) %>%
  mutate(model_id = my_model_id,
         reference_datetime = forecast_date,
         family = 'ensemble',
         duration = 'P1D',
         parameter = as.character(parameter),
         project_id = 'neon4cast') %>%
  select(datetime, reference_datetime, duration, site_id, family, parameter, variable, prediction, model_id, project_id)
#---------------------------#



# ----- Submit forecast -----
# Write the forecast to file
theme <- 'aquatics'
date <- forecast_df_EFI$reference_datetime[1]
forecast_name <- paste0(forecast_df_EFI$model_id[1], ".csv")
forecast_file <- paste(theme, date, forecast_name, sep = '-')

write_csv(forecast_df_EFI, forecast_file)

neon4cast::forecast_output_validator(forecast_file)


neon4cast::submit(forecast_file =  forecast_file, ask = FALSE) # if ask = T (default), it will produce a pop-up box asking if you want to submit

#--------------------------#

forecast_df_EFI |> 
  ggplot(aes(x=datetime, y=prediction, group = parameter)) +
  geom_line() +
  facet_wrap(~site_id) +
  labs(title = paste0('Forecast generated for ', forecast_df_EFI$variable[1], ' on ', forecast_df_EFI$reference_datetime[1]))

plot_file_name <- paste0("Submit_forecast/", forecast_df_EFI$variable[1], '-', forecast_df_EFI$reference_datetime[1], ".png")
ggsave(plot_file_name)


