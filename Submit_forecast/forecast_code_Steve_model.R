## install.packages('remotes')
## install.packages('tidyverse') # collection of R packages for data manipulation, analysis, and visualisation
## install.packages('lubridate') # working with dates and times
## remotes::install_github('eco4cast/neon4cast') # package from NEON4cast challenge organisers to assist with forecast building and submission

# ------ Load packages -----
library(tidyverse)
library(lubridate)
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
            by = c("datetime","site_id"))

forecast_df <- NULL


for(i in 1:length(focal_sites)) {  
  
  curr_site <- focal_sites[i]
  
  site_target <- targets_lm |> filter(site_id == curr_site)
  noaa_future_site <- weather_future_daily |> filter(site_id == curr_site)
  
  
  fit <- lm(temperature ~ air_temperature, data = site_target)
  coeffs <- fit$coefficients
  fit_summary <- summary(fit)
  params_se <- fit_summary$coefficients[,2]
  
  
  residuals <- site_target$temperature - predict(fit)
  sigma <- sd(residuals, na.rm = TRUE)
  
  
  n_members <- 31
  
  # Parameter Uncertanity
  param_df <- tibble(
    parameter = unique(noaa_future_site$parameter),
    beta0 = rnorm(n_members, coeffs[1], params_se[1]),
    beta1 = rnorm(n_members, coeffs[2], params_se[2])
  )
  
  
  forecast_total_unc <- noaa_future_site %>%
    left_join(param_df, by = "parameter") %>%
    mutate(
      forecast_date = datetime,
      value = as.double(NA)
    )
  
 
  forecasted_dates <- unique(forecast_total_unc$forecast_date)
  
  for(j in 1:length(forecasted_dates)) {
    
    temp_pred <- forecast_total_unc %>%
      filter(forecast_date == forecasted_dates[j])
    
    temp_pred <- temp_pred %>%
      mutate(
        value = beta0 +
          beta1 * air_temperature +
          rnorm(n(), mean = 0, sd = sigma) #Process uncertanity
      )
    
    forecast_total_unc <- forecast_total_unc %>%
      rows_update(temp_pred, by = c("forecast_date","parameter"))
  }
  
  curr_site_df <- forecast_total_unc %>%
    transmute(
      datetime = forecast_date,
      site_id = curr_site,
      parameter = as.character(parameter),
      prediction = value,
      variable = "temperature")

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

