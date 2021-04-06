#d18O_evap_functions.R
# All functions required to calculate the isotopic composition of evaporation.
# Includes:
# - calculate_Cevap
# - Cevap_sat_vapor_press
# - Cevap_normalized_humidity
# - Cevap_isotope_frac
# - Cevap_kinetic_frac
# - Cevap_total_frac
# - Cevap_Catm

# ------------------------------------------------------------------------------
#' Evaporation stable isotope
#'
#' Calculates the stable isotope composition of evaporation based on equation 5
#' in Krabbenhoft et al. (1990). Note that alpha* is equivalent to alpha^-1.
#'
#' @references Krabbenhoft, D. P., C. J. Bowser, M. P. Anderson, and J. W.
#'   Valley. (1990). Estimating Groundwater Exchange with Lakes: 1. The Stable
#'   Isotope Mass Balance Method. Water Resources Research, 26(10):2445-2453.
#'   https://doi.org/10.1029/WR026i010p02445
#'
#' @param atmp air temperature (C)
#' @param ltmp lake surface temperature (C)
#' @param RH relative humidity (percent)
#' @param Cpcpn stable isotope composition of precipitation
#' @param Clake stable isotope composition of the lake
#' @param parameter name of stable isotope, "d18O" or "d2H"
#'
#' @return Cevap - the isotope composition of evaporation
#'
#' @importFrom NISTunits NISTdegCtOk
#'
#' @export

calculate_Cevap <- function(atmp, ltmp, RH, Cpcpn, Clake, parameter) {

  # Required parameters
  es_a          <- Cevap_sat_vapor_press(atmp)
  es_l          <- Cevap_sat_vapor_press(ltmp)
  h             <- Cevap_normalized_humidity(RH, es_a, es_l)
  alpha         <- Cevap_isotope_frac(NISTunits::NISTdegCtOk(ltmp))

  if (parameter == "d18O") {
    K = 14.3
  } else if (parameter == "d2H") {
    K = 12.5
  }

  # Seasonal correction factor based on difference between calculated and
  # reported values of d18O_atm in Krabbenhoft et al., 1990.
  # if (seasonal_correction) {
  #   corr <- data.frame(month = c(1, 2, 3, 4,
  #                                5, 6, 7,
  #                                8, 9,
  #                                10, 11, 12),
  #                      factor = c(1, 1, 1, 1,
  #                                 1.0319431, 1.0474563, 1.1181465,
  #                                 1.1587506, 1.0335792,
  #                                 1, 1, 1))
  #   factor <- corr$factor[corr$month == month]
  # } else {
  #   factor <- 1
  # }

  delta_epsilon <- Cevap_kinetic_frac(h, K)
  epsilon       <- Cevap_total_frac(alpha, delta_epsilon)
  Catm          <- Cevap_Catm(Cpcpn, alpha)

  # Evaporation
  Cevap <- ((1/alpha)*Clake - h*Catm - epsilon)/
               (1 - h + delta_epsilon*10^(-3))
  return(Cevap)
}
# ------------------------------------------------------------------------------
#' Saturation Vapor Pressure
#'
#' Calculates the saturation vapor pressure based on temperature (of air or
#' water) based on equations 11 and 12 of Allen et al. (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param tmp temperature of air or water (degrees C)
#'
#' @return es - saturation vapor pressure (kPa)
#'
#' @export

Cevap_sat_vapor_press <- function(tmp) {
  es <- 0.6108*exp(17.27*tmp/(237.3 + tmp))
  return(es)
}
# ------------------------------------------------------------------------------
#' Normalized Relative Humidity
#'
#' Calculates the relative humidity normalized to the temperature of the surface
#' water. Based on Equation 1.8 of Mook (2000).
#'
#' @references Mook, W.G. (ed.) 2000. Environmental Isotopes in the Hydrologic
#'   Cycle: Volume III: Surface Water. UNESCO. Paris, France.
#'
#' @param RH relative humidity (percent)
#' @param es_a saturation vapor pressure for the air (kPa)
#' @param es_l saturation vapor pressure for the lake (kPa)
#'
#' @return h - the relative humidity normalized to the temperature of the
#'         surface water (-)
#'
#' @export

Cevap_normalized_humidity <- function(RH, es_a, es_l) {
  h <- (RH/100) * (es_a/es_l)
  return(h)
}
# ------------------------------------------------------------------------------
#' Equilibrium Isotope Fractionation Factor
#'
#' Calculates the equilibrium isotope fractionation factor at the temperature of
#' the air-water interface based on Eq. 16a in Gibson et al. (2016) or Eq. 1.6
#' in Mook (2000). Returns this value as the ratio in liquid vs. the ratio in
#' vapor (i.e., LV form, alpha > 1). Gibson et al. refer to this formulation as
#' alpha+ and to the VL form (i.e., alpha < 1) as alpha*. Krabbenhoft et al.
#' (1990) use alpha* (i.e, VL form, alpha < 1) in their equations.
#'
#' @references Gibson, J.J., S.J. Birks, and Y. Yi. 2016. Stable isotope mass
#'   balance of lakes: a contemporary perspective. Quaternary Science Reviews,
#'   131:316-328. https://doi.org/10.1016/j.quascirev.2015.04.013
#'
#' @references Mook, W.G. (ed.) 2000. Environmental Isotopes in the Hydrologic
#'   Cycle: Volume III: Surface Water. UNESCO. Paris, France.
#'
#' @references Krabbenhoft, D. P., C. J. Bowser, M. P. Anderson, and J. W.
#'   Valley. (1990). Estimating Groundwater Exchange with Lakes: 1. The Stable
#'   Isotope Mass Balance Method. Water Resources Research, 26(10):2445-2453.
#'   https://doi.org/10.1029/WR026i010p02445
#'
#' @param ltmp lake surface temperature (K)
#' @param method equation to use, defaults to "Gibson" but can also be "Mook".
#'
#' @return alpha (-), the equilibrium isotope fractionation factor in LV form
#'        (i.e., alpha > 1)
#'
#' @export

Cevap_isotope_frac <- function(ltmp, method = "Mook"){
  if (method == "Gibson") {
    alpha <- exp(-7.685e-3 + 6.7123/ltmp - 1666.4/(ltmp^2) + 350410/(ltmp^3))
  } else if (method == "Mook") {
    alpha <- 1/exp(2.0667e-3 + (0.4156/ltmp) - (1.137e3/(ltmp^2)))
  }
  return(alpha)
}
# ------------------------------------------------------------------------------
#' Kinetic Fractionation Factor
#'
#' Calculates the kinetic fractionation factor based on equation 6 in
#' Krabeenhoft et al. (1990).
#'
#' @references Krabbenhoft, D. P., C. J. Bowser, M. P. Anderson, and J. W.
#'   Valley. (1990). Estimating Groundwater Exchange with Lakes: 1. The Stable
#'   Isotope Mass Balance Method. Water Resources Research, 26(10):2445-2453.
#'   https://doi.org/10.1029/WR026i010p02445
#'
#'
#' @param h relative humidity normalized to the temperature of the surface water
#'          (-)
#' @param K constant determinded by wind tunnel experiments for different
#'          isotopes, defaults to K(18O) = 14.3.
#'
#' @return delta_epsilon - the kinetic fractionation factor (-)
#'
#' @export

Cevap_kinetic_frac <- function(h, K = 14.3) {
  delta_epsilon <- K*(1-h)
  return(delta_epsilon)
}
# ------------------------------------------------------------------------------
#' Total Fractionation Factor
#'
#' Calculates the total fractionation factor based on the definition of epsilon
#' for equation 5 in Krabbenhoft et al. (1990). Note that while here we
#' calculate calculates alpha in LV form (i.e, alpha+ > 1), the equation in
#' Krabbenhoft et al. (1990) assumes alpha is in V/L form (i.e., alpha* < 1).
#'
#' @references Krabbenhoft, D. P., C. J. Bowser, M. P. Anderson, and J. W.
#'   Valley. (1990). Estimating Groundwater Exchange with Lakes: 1. The Stable
#'   Isotope Mass Balance Method. Water Resources Research, 26(10):2445-2453.
#'   https://doi.org/10.1029/WR026i010p02445
#'
#' @param alpha equilibrium isotope fractionation factor (-) at the temperature
#'              of the air-water interface in LV form (i.e., alpha > 1).
#' @param delta_epsilon kinetic fractionation factor (-)
#'
#' @return epsilon - the total fractionation factor (-)
#'
#' @export

Cevap_total_frac <- function(alpha, delta_epsilon) {
  epsilon <- 1000*(1 - 1/alpha) + delta_epsilon
  return(epsilon)
}
# ------------------------------------------------------------------------------
#' Atmosphere d18O
#'
#' Calculates the d18O isotope composition of the atmosphere based on Equation
#' 18 and the definition for espilon+ in the explanation for Equation 3 in
#' Gibson et al. (2016). Alternatively, can instead calcuate this value based on
#' Equation 1.10 and the definition for epsilon in the explanation for Equation
#' 1.4 in Mook (2000).
#'
#' @references Gibson, J.J., S.J. Birks, and Y. Yi. 2016. Stable isotope mass
#'   balance of lakes: a contemporary perspective. Quaternary Science Reviews,
#'   131:316-328. https://doi.org/10.1016/j.quascirev.2015.04.013
#'
#' @references Mook, W.G. (ed.) 2000. Environmental Isotopes in the Hydrologic
#'   Cycle: Volume III: Surface Water. UNESCO. Paris, France.
#'
#' @param Cpcpn isotopic composition of precipitation
#' @param alpha equilibrium isotope fractionation factor (-) at the temperature
#'              of the air-water interface in LV form (i.e., alpha > 1).
#' @param method defaults to "Gibson" to use those equations, can also be "Mook".
#' @param k weighted factor which reflects seasonality, ranging from 0.5 for
#'          highly seasonal climates to 1 for non-seasonal climates. Defaults
#'          to 0.8
#'
#' @return Catm - the d18O isotope composition of the atmosphere
#'
#' @export

Cevap_Catm <- function(Cpcpn, alpha, method = "Mook_corrected", k = 1) {
  if (method == "Gibson") {
    epsilon_plus <- (alpha - 1)*1000
    Catm     <- (Cpcpn - k*epsilon_plus)/(1 + k*1e-3*epsilon_plus)
  } else if (method == "Mook") {
    epsilon_star <- (1/alpha) - 1
    Catm     <- (1/alpha)*Cpcpn + epsilon_star
  }  else if (method == "Mook_corrected") {
    epsilon_star <- ((1/alpha) - 1)*1000
    Catm     <- (1/alpha)*Cpcpn + epsilon_star
  }
  return(Catm)
}
