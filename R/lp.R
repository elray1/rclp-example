# Example of using distfromq to calculate a linear pool ensemble

# tidyverse packages
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# covidHubUtils for loading component model forecast files
# https://github.com/reichlab/covidHubUtils
library(covidHubUtils)

# distfromq for approximating a distribution's pdf and cdf from a collection
# of quantiles
# https://github.com/reichlab/distfromq
library(distfromq)

# function for loading truth data
source("./R/load_flu_hosp_data.R")

# function to do a grid search to get quantiles from a cdf
source("./R/get_rclp_quantiles.R")


# The reference_date is the date of the Saturday relative to which week-ahead
# targets are defined.
reference_date <- "2022-04-02"

# Get the list of required locations
required_locations <-
  readr::read_csv(file = "./data/locations.csv") %>%
  dplyr::select("location", "abbreviation")

# The forecast_date is the Monday of forecast creation.
forecast_date <- as.character(as.Date(reference_date) + 2)

# Load component forecasts
component_forecasts <- covidHubUtils::load_forecasts_repo(
  file_path = paste0("component-forecasts/"),
  forecast_dates = forecast_date,
  locations = NULL,
  types = "quantile",
  targets = NULL,
  hub = "FluSight",
  verbose = TRUE
)


# compute ensemble forecasts as linear pool
# three possible approaches:
#  1. from each forecaster, draw random samples from their distribution; take quantiles
#  2. from each forecaster, extract many quantiles from their distribution; take quantiles
#     basically, we're regarding the quantiles from the component forecasters' distributions
#     as random samples. this strategy is essentially equivalent to 1, but with less
#     Monte Carlo variability
#  3. construct a mixture model and solve for quantiles
#     There is related code in get_rclp_quantiles.R, but it would need to be updated;
#     not doing that for now.

# Strategy 1: random samples
n_samples <- 1e4
quantile_levels <- unique(component_forecasts$quantile)

ensemble_r_forecasts <- component_forecasts |>
  dplyr::group_by(model, forecast_date, location, horizon,
                  temporal_resolution, target_variable, target_end_date) |>
  dplyr::summarize(
    pred_samples = list(distfromq::make_r_fn(
      ps = quantile,
      qs = value)(n_samples)),
    .groups = "drop"
  ) |>
  tidyr::unnest(pred_samples) |>
  dplyr::group_by(forecast_date, location, horizon,
                  temporal_resolution, target_variable, target_end_date) |>
  dplyr::summarize(
    quantile = list(quantile_levels),
    value = list(quantile(pred_samples, probs = quantile_levels)),
    .groups = "drop"
  ) |>
  tidyr::unnest(cols = tidyselect::all_of(c("quantile", "value"))) |>
  dplyr::mutate(
    model = "ensemble_r"
  )


# Strategy 1: quantiles
n_samples <- 1e4
quantile_levels <- unique(component_forecasts$quantile)

ensemble_q_forecasts <- component_forecasts |>
  dplyr::group_by(model, forecast_date, location, horizon,
                  temporal_resolution, target_variable, target_end_date) |>
  dplyr::summarize(
    pred_qs = list(distfromq::make_q_fn(
      ps = quantile,
      qs = value)(seq(from = 0, to = 1, length.out = n_samples + 2)[2:n_samples])),
    .groups = "drop"
  ) |>
  tidyr::unnest(pred_qs) |>
  dplyr::group_by(forecast_date, location, horizon,
                  temporal_resolution, target_variable, target_end_date) |>
  dplyr::summarize(
    quantile = list(quantile_levels),
    value = list(quantile(pred_qs, probs = quantile_levels)),
    .groups = "drop"
  ) |>
  tidyr::unnest(cols = tidyselect::all_of(c("quantile", "value"))) |>
  dplyr::mutate(
    model = "ensemble_q"
  )


# make a plot
all_forecasts <- dplyr::bind_rows(
  component_forecasts,
  ensemble_r_forecasts,
  ensemble_q_forecasts)

ggplot() +
  geom_ribbon(
    data = all_forecasts |>
      filter(quantile %in% c(0.025, 0.975), location == "US") |>
      mutate(interval_bound = ifelse(quantile < 0.5, "lower", "upper")) |>
      select(-quantile) |>
      pivot_wider(names_from = "interval_bound", values_from = "value"),
    mapping = aes(x = target_end_date, ymin = lower, ymax = upper),
    fill = "cornflowerblue"
  ) +
  geom_line(
    data = all_forecasts |>
      filter(quantile == 0.5, location == "US"),
    mapping = aes(x = target_end_date, y = value)
  ) +
  facet_wrap( ~ model)

