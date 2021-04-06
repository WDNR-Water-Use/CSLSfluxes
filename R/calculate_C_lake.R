#' Calculate lake concentration of parameter
#'
#' Calculates the concentration of a parameter in the lake given an initial lake
#' concentration, monthly fluxes of precipitation, groundwater, and lake volume,
#' and median concentrations of groundwater and precipitation.
#'
#' @param C_vals data frame with columns for "lake", "C0_lake", "C_pcpn" and
#'                 "C_GWin".
#' @param fluxes data frame with lake, date, level_m, GWin_m3, and P_m3 from
#'               MODFLOW simulations to use for mass balance calculations.
#' @param scenario which MODFLOW scenario is being evaluated here (e.g., "cal",
#'                 "no_irr", "cur_irr")
#' @param sim which MODFLOW simulation is being evaluated here (a number)
#' @param start_date start date for evaluation, defaults to 1985-10-01
#' @param end_date end date for evaluation, defaults to 2018-09-30.
#' @param lakes a list of lakes being evaluated. Defaults to c("Pleasant",
#'              "Long", "Plainfield")
#' @param dt time step of calculations. Defaults to "year", can also be "month",
#'           or "day".
#'
#' @return output, a data frame with the following columns:
#'   \item{scenario}{name of MODFLOW scenario evaluated (e.g., "cal", "no_irr")}
#'   \item{sim}{number of MODFLOW simulation evaluated}
#'   \item{lake}{name of lake, i.e., "Pleasant", "Long", and "Plainfield"}
#'   \item{date}{date of observation (POSIX) used for weather inputs in MODFLOW}
#'   \item{level_m}{lake level (m) from MODFLOW}
#'   \item{vol_m3}{lake volume (m3) from MODFLOW}
#'   \item{P_m3}{precipitation (m3) from MODFLOW}
#'   \item{C_pcpn}{median precipitation concentration (mg/L) from CSLSdata
#'                 observations}
#'   \item{M_P}{calculated mass of parameter in precipitation (g)}
#'   \item{GWin_m3}{groundwater inflow (m3) from MODFLOW}
#'   \item{C_GWin}{median groundwater inflow concentration (mg/L) from CSLSdata
#'                 observations}
#'   \item{M_GWin}{calculated mass of parameter in groundwater inflow (g)}
#'   \item{GWout_m3}{groundwater outflow (m3) from MODFLOW}
#'   \item{M_GWin}{calculated mass of parameter in groundwater outflow (g)}
#'   \item{C0_lake}{initial lake concentration (mg/L) provided to function}
#'   \item{M_lake}{calculated mass of parameter in lake (g)}
#'   \item{C_lake}{calculated concentration of parameter in lake (mg/L)}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr filter mutate select left_join group_by summarise arrange
#' @importFrom stats approxfun median
#' @importFrom lubridate as_datetime month year
#'
#' @export

calculate_C_lake <- function(C_vals,
                             fluxes = CSLSdata::MODFLOW,
                             scenario = "no_irr",
                             sim = 1,
                             start_date = as_datetime("1985-10-01"),
                             end_date = as_datetime("2018-09-30"),
                             lakes = c("Pleasant", "Long", "Plainfield"),
                             dt = "year") {

  # Water Volumes (m3) ---------------------------------------------------------
  if (dt == "day") {
    dates        <- seq(start_date, end_date, "1 day")
    dates        <- data.frame(lake = c(rep("Pleasant", length(dates)),
                                        rep("Long", length(dates)),
                                        rep("Plainfield", length(dates))),
                               day = c(dates, dates, dates)) %>%
                    mutate(date = floor_date(day, unit = "month"))
    fluxes       <- fluxes %>%
                    filter(.data$scenario == !!scenario,
                           .data$sim == !!sim,
                           .data$date >= start_date,
                           .data$date <= end_date + months(1)) %>%
                    mutate(GWin_m3_d = abs(.data$GWin_m3_d)) %>%
                    select(lake = .data$lake,
                           date = .data$date,
                           level_m = .data$level_m,
                           vol_m3 = .data$vol_m3,
                           P_m3 = .data$P_m3_d,
                           GWin_m3 = .data$GWin_m3_d,
                           GWout_m3 = .data$GWout_m3_d)
    monthly_vol  <- fluxes %>% select(.data$lake, .data$date, .data$vol_m3,
                                      .data$level_m)
    daily_vol    <- interpolate_values(monthly_vol,
                                       group_vars = c("lake"),
                                       val_var = "vol_m3",
                                       start_date = start_date,
                                       end_date = end_date)
    daily_level  <- interpolate_values(monthly_vol,
                                       group_vars = c("lake"),
                                       val_var = "level_m",
                                       start_date = start_date,
                                       end_date = end_date)
    daily_fluxes <- left_join(dates, fluxes, by = c("lake", "date")) %>%
                    select(lake = .data$lake,
                           date = .data$day,
                           P_m3 = .data$P_m3,
                           GWin_m3 = .data$GWin_m3,
                           GWout_m3 = .data$GWout_m3)
    fluxes       <- left_join(daily_fluxes, daily_vol,
                              by = c("lake", "date")) %>%
                    left_join(daily_level,
                              by = c("lake", "date")) %>%
                    filter(.data$lake %in% lakes)
  } else {
    fluxes       <- fluxes %>%
                    filter(.data$scenario == !!scenario,
                           .data$sim == !!sim,
                           .data$date >= start_date,
                           .data$date <= end_date + months(1),
                           .data$lake %in% lakes) %>%
                    mutate(GWin_m3 = abs(.data$GWin_m3),
                           water_year = ifelse(month(.data$date) >= 10,
                                               year(.data$date) + 1,
                                               year(.data$date))) %>%
                    select(lake = .data$lake,
                           date = .data$date,
                           water_year = .data$water_year,
                           level_m = .data$level_m,
                           vol_m3 = .data$vol_m3,
                           P_m3 = .data$P_m3,
                           GWin_m3 = .data$GWin_m3,
                           GWout_m3 = .data$GWout_m3)
  }
  if (dt == "year") {
    fluxes <- fluxes %>%
              group_by(lake = .data$lake,
                       water_year = .data$water_year) %>%
              mutate(min_date = min(.data$date),
                     max_date = max(.data$date)) %>%
              summarise(level_m = mean(.data$level_m),
                        vol_m3 = .data$vol_m3[.data$date == .data$min_date],
                        P_m3 = sum(.data$P_m3, na.rm = TRUE),
                        GWin_m3 = sum(.data$GWin_m3, na.rm = TRUE),
                        GWout_m3 = sum(.data$GWout_m3, na.rm = TRUE),
                        date = min(.data$min_date),
                        .groups = "drop") %>%
             select(lake = .data$lake,
                    date = .data$date,
                    water_year = .data$water_year,
                    level_m = .data$level_m,
                    vol_m3 = .data$vol_m3,
                    P_m3 = .data$P_m3,
                    GWin_m3 = .data$GWin_m3,
                    GWout_m3 = .data$GWout_m3)
  }

  # Calculate Masses -----------------------------------------------------------
  input <- left_join(fluxes, C_vals, by = "lake") %>%
            mutate(M_GWin = .data$GWin_m3*.data$C_GWin,
                   M_P = .data$P_m3*.data$C_pcpn)

  # Loop through calcs ---------------------------------------------------------
  input$M_GWout <- NA
  input$M_lake  <- NA
  input$C_lake  <- NA
  output        <- NULL
  for (lake in lakes) {
    this_lake <- input %>%
                 filter(.data$lake == !!lake) %>%
                 arrange(.data$date)
    this_lake$C_lake[1]  <- this_lake$C0_lake[1]
    this_lake$M_lake[1]  <- this_lake$C_lake[1]*this_lake$vol_m3[1]
    for (i in 1:(nrow(this_lake)-1)) {
      this_lake$M_GWout[i]  <- this_lake$C_lake[i]*this_lake$GWout_m3[i]
      this_lake$M_lake[i+1] <- this_lake$M_lake[i] + this_lake$M_P[i] +
                               this_lake$M_GWin[i] - this_lake$M_GWout[i]
      this_lake$C_lake[i+1]  <- this_lake$M_lake[i+1]/this_lake$vol_m3[i+1]
    }
    output <- rbind(output, this_lake)
  }

  output <- output %>%
            mutate(scenario = scenario,
                   sim = sim) %>%
            select(.data$scenario, .data$sim, .data$lake, .data$date,
                   .data$level_m, .data$vol_m3, .data$P_m3, .data$C_pcpn,
                   .data$M_P, .data$GWin_m3, .data$C_GWin, .data$M_GWin,
                   .data$GWout_m3, .data$M_GWout, .data$C0_lake, .data$M_lake,
                   .data$C_lake)

  return(output)

}
