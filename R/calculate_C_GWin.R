#' Calculate lake concentration of parameter
#'
#' Calculates the concentration of a parameter in the lake given an initial lake
#' concentration, monthly fluxes of precipitation, groundwater, and lake volume,
#' and median concentrations of groundwater and precipitation.
#'
#' @param C_vals data frame with columns for "lake", "C_pcpn", "C1_lake",
#'               "C2_lake", and "C_GWout"
#' @param fluxes data frame with lake, date, level_m, GWin_m3, and P_m3 from
#'               MODFLOW simulations to use for mass balance calculations.
#' @param scenario which MODFLOW scenario is being evaluated here, defaults to "cal".
#' @param sim which MODFLOW simulation is being evaluated here, defaults to 0
#' @param start_date start date for evaluation, defaults to 2017-10-01
#' @param end_date end date for evaluation, defaults to 2018-09-30.
#' @param lakes name of lakes to include in analysis, defaults to
#'              c("Pleasant", "Long", "Plainfield")
#'
#' @return C_GWin, a data frame with the following columns:
#'   \item{lake}{name of lakes, i.e., "Pleasant", "Long", and "Plainfield"}
#'   \item{C_GWin}{calculated concentration in groundwater inflow}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr filter mutate select left_join group_by summarise
#' @importFrom lubridate as_datetime
#'
#' @export

calculate_C_GWin <- function(C_vals,
                             fluxes = CSLSdata::MODFLOW,
                             scenario = "cal",
                             sim = 0,
                             start_date = as_datetime("2017-10-01"),
                             end_date = as_datetime("2018-09-30"),
                             lakes = c("Pleasant", "Long", "Plainfield")) {

  # Water Volumes (m3) ---------------------------------------------------------
  fluxes <- fluxes %>%
            filter(.data$scenario == !!scenario,
                   .data$sim == !!sim,
                   .data$date >= start_date,
                   .data$date <= end_date,
                   .data$lake %in% lakes) %>%
            mutate(GWin_m3 = abs(.data$GWin_m3),
                   min_date = min(.data$date),
                   max_date = max(.data$date)) %>%
            group_by(lake = .data$lake) %>%
            summarise(vol1_m3 = .data$vol_m3[.data$date == .data$min_date],
                      vol2_m3 = .data$vol_m3[.data$date == .data$max_date],
                      P_m3 = sum(.data$P_m3, na.rm = TRUE),
                      GWin_m3 = sum(.data$GWin_m3, na.rm = TRUE),
                      GWout_m3 = sum(.data$GWout_m3, na.rm = TRUE),
                      .groups = "drop")


  C_GWin <- left_join(fluxes, C_vals, by = "lake") %>%
            mutate(C_GWin = round((.data$C2_lake*.data$vol2_m3 -
                                     .data$C1_lake*.data$vol1_m3 -
                                     .data$C_pcpn*.data$P_m3 +
                                     .data$C_GWout*.data$GWout_m3)/
                                   (.data$GWin_m3),
                                  2)) %>%
            select(.data$lake, .data$C_pcpn, .data$C1_lake, .data$C2_lake,
                   .data$C_GWout, .data$C_GWin)

  return(C_GWin)

}
