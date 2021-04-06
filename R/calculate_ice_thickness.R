#' Calculate ice thickness
#'
#' Calculates the thickness of ice on a lake. Assumes a linear increase from a
#' static ice-on date to a static ice-off date with a static date of max
#' thickness and a static max thickness.
#'
#' @param date date to calculate ice thickness (POSIX)
#' @param ice_on day of year when ice starts forming. Defaults to day 320 (Nov
#'               16).
#' @param ice_off day of year when ice is all gone. Defaults to day 105 (Apr 15)
#' @param ice_max day of year when ice thickness is at a maxiumum. Defaults to
#'                day 60 (Mar 1).
#' @param ice_max_depth max thickness of ice (m). Defaults to 24 inches.
#'
#' @return ice_depth, a numeric value of ice thickness (m).
#'
#' @importFrom lubridate yday ceiling_date floor_date days
#' @importFrom NISTunits NISTinchTOmeter
#'
#' @export

calculate_ice_thickness <- function(date,
                                    ice_on = 320,
                                    ice_off = 105,
                                    ice_max = 60,
                                    ice_max_depth = NISTunits::NISTinchTOmeter(24)*1000) {
  day            <- yday(date)
  days_in_year   <- yday(ceiling_date(date, unit = "year") - days(1))
  days_last_year <- yday(floor_date(date, unit = "year") - days(1))
  if (day > ice_on){
    ice_depth <- ice_max_depth*(day - ice_on)/
                 (days_in_year - ice_on + ice_max)
  } else if (day <= ice_max){
    ice_depth <- ice_max_depth*(days_last_year - ice_on + day)/
                 (days_last_year - ice_on + ice_max)
  } else if (day < ice_off) {
    ice_depth <- ice_max_depth*(ice_off - day)/
                 (ice_off - ice_max)
  } else {
    ice_depth <- 0
  }
  return(ice_depth)
}
