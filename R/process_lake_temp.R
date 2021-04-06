#' Process lake temperature data for calcs
#'
#' Processes CSLS lake temperature data for water/solute balance calculations.
#' Filters HOBO temperature data from CSLSdata::water_chem for surface and
#' bottom temperatures, and summarizes at a daily or monthly time step for the
#' desired time period.
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
#'   \item{ltmp_bot_C}{lake temperature at bottom of lake, deg C. Used for some
#'                     solute reaction rates.}
#'   \item{ltmp_surf_C}{lake temperature at surface of lake, deg C. Used for
#'                      stable isotope in evaporation calculations}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr group_by summarise
#' @importFrom lubridate floor_date
#'
#' @export

process_lake_temp <- function(start_date, end_date, dt = "day"){

  # Lake bottom temperature
  ltmp_bot   <- filter_ltmp(depth = "bottom")  %>%
                group_by(date = floor_date(.data$date, unit = "day"),
                         lake = .data$lake) %>%
                summarise(ltmp_bot_C = mean(.data$result, na.rm = TRUE)) %>%
                ungroup()
  ltmp_bot   <- interpolate_values(ltmp_bot, "lake", "ltmp_bot_C",
                                   start_date, end_date)

  ltmp_surf  <- filter_ltmp(depth = "surface")  %>%
                group_by(date = floor_date(.data$date, unit = "day"),
                         lake = .data$lake) %>%
                summarise(ltmp_surf_C = mean(.data$result, na.rm = TRUE)) %>%
                ungroup()
  ltmp_surf  <- interpolate_values(ltmp_surf, "lake", "ltmp_surf_C",
                                   start_date, end_date)

  if (dt == "day") {
    df  <- full_join(ltmp_bot, ltmp_surf, by = c("lake", "date"))
  } else if (dt == "month") {
    df  <- full_join(ltmp_bot, ltmp_surf, by = c("lake", "date")) %>%
           group_by(lake = .data$lake,
                    date = floor_date(.data$date, unit = "month")) %>%
           summarise(ltmp_bot_C = mean(.data$ltmp_bot_C, na.rm = TRUE),
                     ltmp_surf_C = mean(.data$ltmp_surf_C, na.rm = TRUE)) %>%
           ungroup()
  }
  return(df)
}
