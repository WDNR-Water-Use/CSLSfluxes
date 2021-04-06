#' Interpolate values
#'
#' Interpolates values to a daily time step given start and end date of desired
#' timeseries. Can summarize at a larger time step (e.g., monthly) if desired
#' using final_dt flag.
#'
#' @param df data frame with water chemistry information for tracer/solute.
#'           Includes:
#'   * **date:** date of measurement
#'   * **group1 and/or group2:** grouping for measurements
#'   * **value:** value of measurement
#' @param group_vars name of group1 and/or group2 columns, for identifying and
#'                   renaming them.
#' @param val_var name of value column, for identifying and renaming it
#' @param start_date start date of interpolated timeseries (POSIX). If earlier
#'                   than earliest date in df, values before earliest date in df
#'                   will have value equal to earliest date in df.
#' @param end_date end date of interpolated timeseries (POSIX). If later than
#'                 latest date in df, values after latest date in df will have
#'                 value equal to latest date in df.
#' @param final_dt unit of time to use for interpolation. Defaults to "day", can
#'           also be "month". Used as the \code{unit} for
#'          \code{lubridate::floor_date}
#'
#' @return df, a data frame with interpolated values.
#'
#' @import lubridate
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr arrange filter mutate
#' @importFrom zoo read.zoo na.approx
#'
#' @export

interpolate_values <- function(df,
                               group_vars,
                               val_var,
                               start_date,
                               end_date,
                               final_dt = "day"){

  # Rename columns w/generic group/value names
  for (i in 1:length(group_vars)) {
    colnames(df)[colnames(df) == group_vars[i]] <- sprintf("group%d", i)
  }
  colnames(df)[colnames(df) == val_var] <- "value"


  # Define dates timeseries
  dates <- data.frame(date = seq(start_date, end_date, by = '1 day'))

  # Interpolate
  interpolated <- NULL
  if (length(group_vars) == 1) {
    # Only one grouping variable (e.g., lake)
    for (group1 in unique(df$group1)) {
      subset        <- df %>%
                       filter(.data$group1 == !!group1) %>%
                       full_join(dates, by = "date") %>%
                       arrange(.data$date)%>%
                       mutate(group1 = !!group1)
      zoo_subset    <- read.zoo(select(subset, c("date", "value")),
                                index.name = "date")
      zoo_subset    <- data.frame(value = na.approx(zoo_subset, rule = 2))
      subset$value  <- zoo_subset$value
      interpolated  <- rbind(interpolated, subset)
    }
    df <- interpolated %>%
          filter(.data$date %in% dates$date) %>%
          group_by(date = floor_date(.data$date, unit = final_dt),
                   group1 = .data$group1) %>%
          summarise(value = mean(.data$value, na.rm = TRUE)) %>%
          ungroup()
  } else if (length(group_vars) == 2) {
    # Two grouping variables (e.g., lake and site_type)
    for (group1 in unique(df$group1)) {
      for (group2 in unique(filter(df, .data$group1 == !!group1)$group2)) {
      subset        <- df %>%
                       filter(.data$group1 == !!group1,
                              .data$group2 == !!group2) %>%
                       full_join(dates, by = "date") %>%
                       arrange(.data$date) %>%
                       mutate(group1 = !!group1,
                              group2 = !!group2)
      zoo_subset    <- read.zoo(select(subset, c("date", "value")),
                                index.name = "date")
      zoo_subset    <- data.frame(value = na.approx(zoo_subset, rule = 2))
      subset$value  <- zoo_subset$value
      interpolated  <- rbind(interpolated, subset)
      }
    }
    df <- interpolated %>%
          filter(.data$date %in% dates$date) %>%
          group_by(date = floor_date(.data$date, unit = final_dt),
                   group1 = .data$group1,
                   group2 = .data$group2) %>%
          summarise(value = mean(.data$value, na.rm = TRUE)) %>%
          ungroup()
  } else {
    stop("More than two groups supplied to interpolate values, function needs to be updated to handle this.")
  }

  # Rename columns w/actual group/value names
  for (i in 1:length(group_vars)) {
    colnames(df)[colnames(df) == sprintf("group%d", i)] <- group_vars[i]
  }
  colnames(df)[colnames(df) == "value"] <- val_var

  return(df)
}
