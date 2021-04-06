#' Calculate lake concentration
#'
#' Calculates the concentration of a parameter in the lake across
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
#' @param summary_type Either "median", "min", or "max". Defaults to "median".
#' @param value_name name of column with resulting values. Defaults to "C_lake".
#'
#' @return C_pcpn, a data frame with the following columns:
#'   \item{lake}{name of lakes, i.e., "Pleasant", "Long", and "Plainfield"}
#'   \item{C_lake}{summary value of parameter in lake to use}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr filter summarise group_by
#' @importFrom stats median
#'
#' @export

get_C_lake <- function(param = "MAGNESIUM TOTAL RECOVERABLE",
                       water_chem = CSLSdata::water_chem,
                       start_date = NULL,
                       end_date = NULL,
                       summary_type = "median",
                       value_name = "C_lake") {

  # Get lake chemistry
  chem   <- filter_parameter(water_chem, param) %>%
            filter(.data$site_type == "lake")

  # Limit dates, if desired
  if (!is.null(start_date)) {
    chem <- chem %>% filter(.data$date >= start_date)
  }
  if (!is.null(end_date)) {
    chem <- chem %>% filter(.data$date <= end_date)
  }

  if (summary_type == "median") {
    C_lake <- chem %>%
              group_by(.data$lake) %>%
              summarise(C_lake = median(.data$result, na.rm = TRUE),
                        .groups = "drop")
  } else if (summary_type == "min") {
    C_lake <- chem %>%
              group_by(.data$lake) %>%
              summarise(C_lake = min(.data$result, na.rm = TRUE),
                        .groups = "drop")
  } else if (summary_type == "max") {
    C_lake <- chem %>%
              group_by(.data$lake) %>%
              summarise(C_lake = max(.data$result, na.rm = TRUE),
                        .groups = "drop")
  } else {
    warning(sprintf("unrecognized summary_type %s", summary_type))
  }

  colnames(C_lake) <- c("lake", value_name)

  return(C_lake)

}
