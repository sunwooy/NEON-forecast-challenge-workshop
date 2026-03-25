## install.packages('remotes')
## install.packages('tidyverse') # collection of R packages for data manipulation, analysis, and visualisation
## install.packages('lubridate') # working with dates and times
## install.packages('zoo')
## remotes::install_github('eco4cast/neon4cast') # package from NEON4cast challenge organisers to assist with forecast building and submission

# ------ Load packages -----
library(tidyverse)
library(lubridate)
library(zoo)

#--------------------------#

# Change this for your model ID
# Include the word "example" in my_model_id for a test submission
# Don't include the word "example" in my_model_id for a forecast that you have registered (see neon4cast.org for the registration form)
my_model_id <- 'example_Model_steve_1'

# --Model description--- #

# Add a brief description of your modeling approach

# -- Uncertainty representation -- #

# Describe what sources of uncertainty are included in your forecast and how you estimate each source.

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
targets <- targets %>%
  filter(site_id %in% focal_sites,
         variable == 'temperature')
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
  select(any_of(c('datetime', 'site_id', met_variables, 'parameter'))) |>
  filter(datetime <= min(datetime) + days(6)) 



targets_lm <- targets |> 
  pivot_wider(names_from = 'variable', values_from = 'observation') |> 
  left_join(weather_past_daily, 
            by = c("datetime","site_id"))|>
  mutate(wtemp_yday = lag(temperature))|>
  mutate(airtemp_yday = rollapplyr(air_temperature, width=3, FUN=mean, fill=NA, align="right", na.rm=TRUE))


forecast_df <- NULL


for(i in 1:length(focal_sites)) {  
  
  curr_site <- focal_sites[i]
  
  site_target <- targets_lm |> 
    filter(site_id == curr_site) |> 
    arrange(datetime)
  # Future weather Forecasts from NOAA
  noaa_future_site <- weather_future_daily |> 
    mutate(parameter = as.character(parameter)) |>
    filter(site_id == curr_site) |> 
    arrange(parameter, datetime)
  
  # Remove NAs to remove gap for lag air temp and water temp
  fit_data <- site_target |> 
    filter(
      !is.na(temperature),
      !is.na(air_temperature),
      !is.na(wtemp_yday),
      !is.na(airtemp_yday)
    )
  
  fit <- lm(temperature ~ air_temperature + wtemp_yday + airtemp_yday, data = fit_data)
  coeffs <- fit$coefficients
  fit_summary <- summary(fit)
  params_se <- fit_summary$coefficients[,2]
  
  residuals <- fit$residuals
  sigma <- sd(residuals, na.rm = TRUE)
  
  n_members <- length(unique(noaa_future_site$parameter))
  
  # Parameter Uncertainty
  param_df <- tibble(
    parameter = unique(noaa_future_site$parameter),
    beta0 = rnorm(n_members, coeffs[1], params_se[1]),
    beta1 = rnorm(n_members, coeffs[2], params_se[2]),
    beta2 = rnorm(n_members, coeffs[3], params_se[3]),
    beta3 = rnorm(n_members, coeffs[4], params_se[4])
  )
  
  
  
  # Pull the most recent observed lake temp before start of forecast for the lag in water temp
  last_obs_wtemp <- site_target |> 
    filter(datetime < forecast_date, !is.na(temperature)) |> 
    arrange(datetime) |> 
    slice_tail(n = 1) |> 
    pull(temperature)
  
  # Build df of Initial condition uncertainty
  params <- unique(noaa_future_site$parameter)
  ic_sd <- 0.1   
  ic_df <- tibble(
    parameter = params,
    ic_value = rnorm(n = n_members, mean = last_obs_wtemp, sd = ic_sd)
  )
  
  # Historical air temp for 3 day rolling air temp for later
  air_hist <- site_target |> 
    filter(!is.na(air_temperature)) |> 
    arrange(datetime) |> 
    select(datetime, air_temperature)
  
  site_forecasts <- list()
  # For each NOAA ensemble member
  for(p in unique(noaa_future_site$parameter)) {
    
    # Pull future air temperature from NOAA
    future_member <- noaa_future_site |> 
      filter(parameter == p) |> 
      arrange(datetime)
    
    # Pull uncertainty calculated earlier
    betas <- param_df |> 
      filter(parameter == p)
    
    beta0 <- betas$beta0
    beta1 <- betas$beta1
    beta2 <- betas$beta2
    beta3 <- betas$beta3
    
    # Pull  yesterday's water temp WITH initial condition uncertanity
    prev_wtemp <- ic_df |> 
      filter(parameter == p) |> 
      pull(ic_value)
    
    # Since we're taking 3 day previous air temp INCLUDING NOAA's forecasted temp, merge past and future air temp
    full_air <- bind_rows(
      air_hist,
      future_member |> 
      select(datetime, air_temperature)) |> 
      arrange(datetime)
    
    preds <- numeric(nrow(future_member))
    
    for(j in 1:nrow(future_member)) {
      
      curr_date <- future_member$datetime[j]
      curr_air <- future_member$air_temperature[j]
      
      # Previous 3 days air temp mean then calculate mean
      prev3_air <- full_air |> 
        filter(datetime < curr_date) |> 
        arrange(datetime) |> 
        slice_tail(n = 3) |> 
        pull(air_temperature)
      
      airtemp_yday <- mean(prev3_air, na.rm = TRUE)
      
      pred_val <- beta0 +
        beta1 * curr_air +
        beta2 * prev_wtemp +
        beta3 * airtemp_yday +
        rnorm(1, mean = 0, sd = sigma)   # process uncertainty
      
      preds[j] <- pred_val
      
      # update lagged water temp for next day with forecasted water temperature
      prev_wtemp <- pred_val
    }
    
    member_df <- future_member |> 
      mutate(prediction = preds, variable = "temperature") |> 
      select(datetime, site_id, parameter, variable, prediction)
    site_forecasts[[as.character(p)]] <- member_df
  }
  
  curr_site_df <- bind_rows(site_forecasts)
  forecast_df <- bind_rows(forecast_df, curr_site_df)
  
  message(curr_site, " forecast complete")
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
         parameter = parameter,
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

