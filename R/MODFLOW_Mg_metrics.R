#' Dataset: MODFLOW Mg budget metrics
#'
#' Mg budget metrics calculated using all simulation results from monte carlo
#' MODFLOW scenarios
#'
#' @docType data
#'
#' @usage data(MODFLOW_Mg_metrics)
#'
#' @format A data frame with the following columns.
#' \describe{
#'   \item{lake}{name of lake, e.g., Pleasant, Long, Plainfield}
#'   \item{metric}{name of hydrologic metric, i.e. solute_budget}
#'   \item{variable}{name of variation on metrics, e.g. median, max, q10}
#'   \item{value}{value of hydrologic metric}
#'   \item{scenario}{MODFLOW scenario (e.g., "cur_irr" or "no_irr")}
#'   \item{sim}{id of MODFLOW simulation}
#' }
"MODFLOW_Mg_metrics"
