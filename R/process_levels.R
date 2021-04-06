#' Process lake elevation, area, volume data for calcs
#'
#' Processes CSLS lake level data for water/solute balance calculations. Ensures
#' no gaps in daily lake elevation data from CSLSdata::lake_levels, uses
#' CSLSdata::bathymetry to convert daily elevations to areas and volumes,
#' then summarizes to desired time step (daily or monthly).
#'
#' If using with non-CSLS data, will need an alternate way to get lake level,
#' lake area, and lake volume data into the format outputted by this function.
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
#'   \item{level_m}{mean lake elevation (mamsl) for the time step (daily or
#'                  monthly)}
#'   \item{area_m2}{mean lake area (m2) for the time step (daily or monthly)}
#'   \item{vol_m3}{mean lake volume (m3) for the time step (daily or monthly)}
#'   \item{dV_m3}{change in lake volume (m3) for the time step (daily or
#'                monthly)}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr filter mutate select group_by arrange summarise lag
#' @importFrom lubridate floor_date
#' @importFrom stats approxfun
#'
#' @export

process_levels <- function(start_date, end_date, dt = "month"){
  # Get lake levels (elevations)
  lake_levels  <- CSLSdata::lake_levels %>%
                  select(.data$date, .data$lake, .data$level_m)

  # Interpolate any missing records (to daily time step)
  lake_levels <- interpolate_values(lake_levels, "lake", "level_m",
                                    start_date, end_date)

  # Convert daily elevations to areas and volumes
  bathymetry       <- CSLSdata::bathymetry
  daily_bathy <- NULL
  for (lake in unique(lake_levels$lake))  {
    this_elev_area_vol <- bathymetry %>% filter(.data$lake == !!lake)
    f_elev_area        <- approxfun(x = this_elev_area_vol$elev_m,
                                    y = this_elev_area_vol$area_m2)
    f_elev_vol         <- approxfun(x = this_elev_area_vol$elev_m,
                                    y = this_elev_area_vol$vol_m3)
    these_levels       <- lake_levels %>%
                          filter(.data$lake == !!lake) %>%
                          mutate(area_m2 = f_elev_area(.data$level_m),
                                 vol_m3 = f_elev_vol(.data$level_m))
    daily_bathy <- rbind(daily_bathy, these_levels)
  }
  daily_bathy   <- daily_bathy %>%
                   group_by(.data$lake) %>%
                   arrange(.data$date) %>%
                   mutate(dV_m3 = .data$vol_m3 - lag(.data$vol_m3)) %>%
                   ungroup()

  # Finalize data frame based on time step
  if (dt == "day") {
    # Daily time step
    df <- daily_bathy
  } else if (dt == "month") {
    # Monthly time step
    df <- daily_bathy %>%
          group_by(lake = .data$lake,
                   date = floor_date(.data$date, unit = "month")) %>%
          summarise(level_m = mean(.data$level_m),
                    area_m2 = mean(.data$area_m2),
                    vol_m3 = mean(.data$vol_m3),
                    dV_m3 = sum(.data$dV_m3, na.rm = TRUE)) %>%
          ungroup()
  }

  return(df)
}
