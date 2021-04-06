#' Calculate precipitation concentration
#'
#' Calculates the median concentration of a parameter in precipitation across
#' all samples (default) or between given dates.
#'
#' @param param description of parameter to evaluate, defaults to
#'              "MAGNESIUM TOTAL RECOVERABLE"
#' @param water_chem data frame with water chemistry to evaluate in
#'                   "filter_parameter", defaults toCSLSdata::water_chem.
#' @param start_date start date for evaluation, defaults to NULL to start with
#'                   first sample.
#' @param end_date end date for evaluation, defaults to NULL to end with last
#'                 sample.
#' @param lakes lakes to map precipitation values to.
#'
#' @return C_pcpn, a data frame with the following columns:
#'   \item{lake}{name of lake, i.e., "Pleasant", "Long", and "Plainfield"}
#'   \item{C_pcpn}{median value of parameter in precipitation to use}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr filter summarise
#' @importFrom stats median
#'
#' @export

get_C_pcpn <- function(param = "MAGNESIUM TOTAL RECOVERABLE",
                       water_chem = CSLSdata::water_chem,
                       start_date = NULL,
                       end_date = NULL,
                       lakes = c("Pleasant", "Long", "Plainfield")) {

  # Get precipitation chemistry
  chem   <- filter_parameter(water_chem, param) %>%
            filter(.data$site_type == "precipitation")

  # Limit dates, if desired
  if (!is.null(start_date)) {
    chem <- chem %>% filter(.data$date >= start_date)
  }
  if (!is.null(end_date)) {
    chem <- chem %>% filter(.data$date <= end_date)
  }

  # Take median value
  median <- median(chem$result, na.rm = TRUE)

  # Assign to all "lakes"
  C_pcpn <- data.frame(lake = lakes,
                       C_pcpn = median)
  return(C_pcpn)

}
