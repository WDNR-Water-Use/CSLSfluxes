#' Add precipitation isotopes to dataset
#'
#' Adds precipitation stable isotope measurements from Maribeth Kniffin (May
#' 2016 - April 2017) and incorporates them in the current timeseries, but only
#' for months without CSLS precipitation measurements.
#'
#' @param tracer data frame with water chemistry information for stable isotope.
#' @param start_date start date of analysis period, POSIXct
#' @param end_date end date of analysis period, POSIXct
#'
#' @return tracer, same data frame as provided, but with additional
#'   precipitation measurements during the analysis period for months with
#'   missing precipitation measurements.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr filter summarise mutate select
#' @import lubridate
#'
#' @export

add_pcpn_isotopes <- function(tracer, start_date, end_date){

  # Separate CSLS and Kniffin precipitation samples
  kniffin <- tracer %>%
             filter(.data$site_type == "precipitation",
                    year(.data$date) < 2018)
  csls    <- tracer %>%
             filter(.data$site_type == "precipitation",
                    year(.data$date) > 2017)

  # Loop through all months in csls timeseries
  add_pcpn  <- NULL
  this_date <- start_date + months(1) - days(1)
  while (this_date <= end_date) {
    csls_month <- csls %>% filter(month(.data$date) == month(this_date))
    # If no CSLS samples this month, add Kniffin measurement (by adjusting date)
    if (nrow(csls_month) == 0) {
      kniffin_month <- kniffin %>%
                       filter(month(.data$date) == month(this_date)) %>%
                       summarise(result = mean(.data$result, na.rm = TRUE)) %>%
                       mutate(date = this_date,
                              lake = "Precip",
                              site_type = "precipitation") %>%
                       select(.data$lake, .data$date, .data$site_type, .data$result)
      add_pcpn      <- rbind(add_pcpn, kniffin_month)
    }
    this_date <- floor_date(this_date %m+% months(1), unit = "month") +
                 months(1) - days(1)
  }

  # If no Kniffin data in CSLS gaps, revert to NA
  add_pcpn$result[is.nan(add_pcpn$result)] <- NA
  tracer <- rbind(tracer, add_pcpn)

  return(tracer)
}
