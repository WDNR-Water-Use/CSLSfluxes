#' Process weather data for calcs
#'
#' Processes CSLS weather data for water/solute balance calculations. Transforms
#' data from CSLSdata package, calculates lake evaporation using CSLSevap
#' package, and summarizes at a daily or monthly time step for the desired time
#' period.
#'
#' If using with non-CSLS data, will need an alternate way to get weather data
#' (including lake evaporation) into the format outputted by this function.
#'
#' @param start_date start date of analysis period (POSIX).
#' @param end_date end date of analysis period (POSIX).
#' @param dt time step at which to summarize data. Defaults to "day" for daily
#'           time step, can also be "month" for monthly time step.
#'
#' @return df, a data frame with the following columns:
#'   \item{lake}{name of lake, i.e., "Pleasant", "Long", and "Plainfield"}
#'   \item{date}{date of observation (POSIX). If monthly time step, monthly
#'               summary is assigned to first day of the month}
#'   \item{day}{day of year of observation (1-366)}
#'   \item{P_mm}{precipitation (mm)}
#'   \item{E_mm}{lake evaporation (mm)}
#'   \item{atmp_C}{mean air temperature (deg C)}
#'   \item{RH_pct}{mean relative humidity (percent)}
#'   \item{irr_factor}{irradiance factor, fraction from 0-1 representing
#'                     relative intensity of solar radiation based on day of
#'                     year. For use in dynamic lake model for some solute
#'                     reactions}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr filter mutate select
#' @importFrom lubridate yday floor_date
#'
#' @export

process_weather <- function(start_date, end_date, dt = "day"){
  # * Use CSLSevap package to calculate lake evaporation for all 3 CSLS lakes
  # * Returns daily atmp_min/atmp_max (C), RH_min/RH_max (%), P (mm) and E (mm)
  # * CSLS_daily_met uses CSLSdata::weather, which has been pre-processed to
  #   ensure no data gaps.
  daily_weather <- CSLSevap::CSLS_daily_met()

  # Summarize min/max atmp and RH as mean values
  # Add in day of year and irradiance factor
  if (dt == "day") {
    # Daily time step
    df  <- daily_weather %>%
           filter(.data$date >= start_date,
                  .data$date <= end_date) %>%
           group_by(.data$lake, .data$date) %>%
           mutate(P_mm = .data$P,
                  E_mm = .data$E,
                  atmp_C = mean(c(.data$atmp_min, .data$atmp_max),
                                na.rm = TRUE),
                  RH_pct = mean(c(.data$RH_min, .data$RH_max),
                                na.rm = TRUE),
                  day = yday(.data$date),
                  irr_factor = ifelse(sin(.data$day/100) > 0,
                                      sin(.data$day/100), 0)) %>%
           ungroup() %>%
           select(.data$lake, .data$date, .data$day, .data$P_mm, .data$E_mm,
                  .data$atmp_C, .data$RH_pct, .data$irr_factor)
  } else if (dt == "month") {
    # Monthly time step
    df <- daily_weather %>%
          filter(.data$date >= start_date,
                 .data$date <= end_date) %>%
          group_by(lake = .data$lake,
                   date = floor_date(.data$date, unit = "month")) %>%
          summarise(P_mm = sum(.data$P),
                    E_mm = sum(.data$E),
                    atmp_C = mean(c(mean(.data$atmp_min, na.rm = TRUE),
                                    mean(.data$atmp_max, na.rm = TRUE))),
                    RH_pct = mean(c(mean(.data$RH_min, na.rm = TRUE),
                                    mean(.data$RH_max, na.rm = TRUE)))) %>%
          ungroup() %>%
          mutate(day = yday(.data$date),
                 irr_factor = ifelse(sin(.data$day/100) > 0,
                                     sin(.data$day/100), 0)) %>%
          select(.data$lake, .data$date, .data$day, .data$P_mm, .data$E_mm,
                 .data$atmp_C, .data$RH_pct, .data$irr_factor)
  }
  return(df)
}
