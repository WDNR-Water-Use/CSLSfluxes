#' Filter SWIMS by parameter
#'
#' Subsets SWIMS dataset for a CSLS lake to only measurements for a parameter of
#' interest.
#'
#' @param water_chem data frame with water_chem information for sites of interest
#' @param param parameter to subset by, must exactly match DNR description
#'                    (e.g., "ALUMINUM,TOTAL RECOVERABLE"). Can be a vector of
#'                    multiple descriptors.
#' @param plotting_name - name of parameter to use for plotting, e.g., "Total
#'                        Recoverable Alumninum", no units. Defaults to "".
#' @param numeric logical defaults to TRUE to indicate results should be numeric
#' @param no_blanks logical defaults to TRUE to exclude blank samples
#' @param no_dups logical defaults to TRUE to exclude duplicate samples
#' @param no_age logical defaults to FALSE to include samples that were analyzed
#'               past the holding date.
#' @param no_bad_well logical defaults to FALSE to include samples from wells
#'                    known to behave weirdly.
#' @param no_bad_sample logical defaults to TRUE to exclude bad measurements
#'                      from QC.
#' @param no_lod logical defaults to FALSE to include samples with results below
#'               the limit of detection.
#' @param no_comment logical defaults to FALSE to include samples with other
#'                   comments from the lab
#' @param note_lake_bottom logical defaults to FALSE to leave all lake samples
#'                         as "lake". If true, notes samples > 5m deep as
#'                         "lake_bottom".
#'
#' @return df, same data frame, with only the results for the
#'         parameter of interest. Also includes the given plotting name ("name").
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr filter mutate
#'
#' @export

filter_parameter <- function(water_chem = CSLSdata::water_chem,
                             param, plotting_name = "",
                             numeric = TRUE, no_blanks = TRUE, no_dups = TRUE,
                             no_age = FALSE, no_bad_well = FALSE,
                             no_bad_sample = TRUE, no_lod = FALSE,
                             no_comment = FALSE, note_lake_bottom = FALSE){
  # Filter to desired parameters -----------------------------------------------
  df <- water_chem %>%
        filter(.data$description %in% param) %>%
        mutate(result = .data$result,
               name = plotting_name)
  if (numeric) {
    df$result <- suppressWarnings(as.numeric(df$result))
  }

  # Remove blanks, duplicates, and old samples ---------------------------------
  if (no_blanks) {
    df <- df %>% filter(.data$flag != "BLANK")
  }
  if (no_dups) {
    df <- df %>% filter(.data$flag != "DUPLICATE")
  }
  if (no_age) {
    df <- df %>% filter(.data$flag != "AGE")
  }
  if (no_bad_well) {
    df <- df %>% filter(.data$flag != "BAD_WELL")
  }
  if (no_bad_sample) {
    df <- df %>% filter(.data$flag != "BAD_SAMPLE")
  }
  if (no_lod) {
    df <- df %>% filter(.data$flag != "LOD")
  }
  if (no_comment) {
    df <- df %>% filter(.data$flag != "COMMENT")
  }

  # Convert units, if needed ---------------------------------------------------
  # Alkalinity: 20 ueq/L = 1 mg/L CaCO3
  if ("ALKALINITY TOTAL GRAN AVAL UEQ/L" %in% param) {
    df <- df %>%
          mutate(result = ifelse(.data$units == "ueq/L",
                                 .data$result/20,
                                 .data$result))
  }

  # Note hypolimnion -----------------------------------------------------------
  if (note_lake_bottom){
    if (!is.null(df$site_type[df$site_type == "lake" &
                              (df$depth1_m > 5 | df$depth2_m > 5)])) {
      df$site_type <- as.character(df$site_type)
      df$site_type[df$site_type == "lake" &
                     (df$depth1_m > 5 | df$depth2_m > 5)] <- "lake_bottom"
      df$site_type <- as.factor(df$site_type)
    }
  }

  return(df)
}
