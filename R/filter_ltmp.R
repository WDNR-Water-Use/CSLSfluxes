#' Filter to lake surface or bottom temperature
#'
#' Filters large dataset of water chemistry parameters to only lake surface or
#' lake bottom temperature measurements.
#'
#' @param chem_df data frame with water chemistry information for all sites.
#'                Defaults to CSLSdata::water_chem.
#' @param depth indicates which depth to use. "surface" pulls shallowest
#'              records, "bottom" pulls deepest records. Defaults to "shallow".
#' @param use_HOBO logical defaults to TRUE to use HOBO data for temperature
#'                 measurements. if false, uses field profile data.
#'
#' @return lst, same data frame as chem_df, but subset to only lake surface
#'         temperature measurements.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr filter group_by mutate
#'
#' @export

filter_ltmp <- function(chem_df = CSLSdata::water_chem,
                        depth = "surface",
                        use_HOBO = TRUE){

  # Grab HOBO or field temperature profile data
  if (use_HOBO){
    ltmp <- filter_parameter(chem_df, "TEMPERATURE HOBO")
  } else {
    ltmp <- filter_parameter(chem_df, "TEMPERATURE FIELD")
  }

  # Keep only shallowest measurement for each date/time
  if (depth == "surface") {
    ltmp <- ltmp %>%
            filter(!is.na(.data$result)) %>%
            group_by(.data$lake, .data$date) %>%
            mutate(min_depth = min(.data$depth1_m)) %>%
            ungroup() %>%
            filter(.data$depth1_m == .data$min_depth)
  } else if (depth == "bottom") {
    ltmp <- ltmp %>%
            filter(!is.na(.data$result)) %>%
            group_by(.data$lake, .data$date) %>%
            mutate(max_depth = max(.data$depth1_m)) %>%
            ungroup() %>%
            filter(.data$depth1_m == .data$max_depth)
  } else {
    stop("Must enter depth of 'surface' or 'bottom'.")
  }

  # Remove extraneous columns
  ltmp$min_depth <- NULL
  ltmp$max_depth <- NULL

  return(ltmp)
}
