% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_C_lake.R
\name{calculate_C_lake}
\alias{calculate_C_lake}
\title{Calculate lake concentration of parameter}
\usage{
calculate_C_lake(
  C_vals,
  fluxes = CSLSdata::MODFLOW,
  scenario = "no_irr",
  sim = 1,
  start_date = as_datetime("1985-10-01"),
  end_date = as_datetime("2018-09-30"),
  lakes = c("Pleasant", "Long", "Plainfield"),
  dt = "year"
)
}
\arguments{
\item{C_vals}{data frame with columns for "lake", "C0_lake", "C_pcpn" and
"C_GWin".}

\item{fluxes}{data frame with lake, date, level_m, GWin_m3, and P_m3 from
MODFLOW simulations to use for mass balance calculations.}

\item{scenario}{which MODFLOW scenario is being evaluated here (e.g., "cal",
"no_irr", "cur_irr")}

\item{sim}{which MODFLOW simulation is being evaluated here (a number)}

\item{start_date}{start date for evaluation, defaults to 1985-10-01}

\item{end_date}{end date for evaluation, defaults to 2018-09-30.}

\item{lakes}{a list of lakes being evaluated. Defaults to c("Pleasant",
"Long", "Plainfield")}

\item{dt}{time step of calculations. Defaults to "year", can also be "month",
or "day".}
}
\value{
output, a data frame with the following columns:
  \item{scenario}{name of MODFLOW scenario evaluated (e.g., "cal", "no_irr")}
  \item{sim}{number of MODFLOW simulation evaluated}
  \item{lake}{name of lake, i.e., "Pleasant", "Long", and "Plainfield"}
  \item{date}{date of observation (POSIX) used for weather inputs in MODFLOW}
  \item{level_m}{lake level (m) from MODFLOW}
  \item{vol_m3}{lake volume (m3) from MODFLOW}
  \item{P_m3}{precipitation (m3) from MODFLOW}
  \item{C_pcpn}{median precipitation concentration (mg/L) from CSLSdata
                observations}
  \item{M_P}{calculated mass of parameter in precipitation (g)}
  \item{GWin_m3}{groundwater inflow (m3) from MODFLOW}
  \item{C_GWin}{median groundwater inflow concentration (mg/L) from CSLSdata
                observations}
  \item{M_GWin}{calculated mass of parameter in groundwater inflow (g)}
  \item{GWout_m3}{groundwater outflow (m3) from MODFLOW}
  \item{M_GWin}{calculated mass of parameter in groundwater outflow (g)}
  \item{C0_lake}{initial lake concentration (mg/L) provided to function}
  \item{M_lake}{calculated mass of parameter in lake (g)}
  \item{C_lake}{calculated concentration of parameter in lake (mg/L)}
}
\description{
Calculates the concentration of a parameter in the lake given an initial lake
concentration, monthly fluxes of precipitation, groundwater, and lake volume,
and median concentrations of groundwater and precipitation.
}
