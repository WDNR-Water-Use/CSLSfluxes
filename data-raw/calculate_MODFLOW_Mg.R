# Evaluate Mg budget

# This generates a REALLY large list/data frame: MODFLOW_Mg. In order to get
# through this script, I need to exit (while saving my workspace) after the
# initial loop. Re-open, convert list to data frame, exit (save workspace).
# Re-open, save as Rda, exit (save workspace).

library(dplyr)
library(lubridate)
library(reshape2)
library(CSLSdata)
library(CSLSscenarios)
library(CSLSfluxes)

param     <- "MAGNESIUM TOTAL RECOVERABLE"
lakes     <- c("Pleasant", "Long", "Plainfield")
scenarios <- c("no_irr", "cur_irr", "wells_off")

# 1. Select upper/lower bound scenarios ----------------------------------------
MODFLOW_metrics  <- CSLSscenarios::MODFLOW_metrics %>%
                    filter(.data$series == "month")

use_sims <- select_bounds(MODFLOW_metrics,
                          base_scenario = "no_irr",
                          compare_scenario = "cur_irr")

# 2. Tune GWin concentrations --------------------------------------------------
C_pcpn  <- get_C_pcpn(param, lakes = lakes)
C_GWout <- get_C_lake(param,
                      summary_type = "median",
                      value_name = "C_GWout")
C1_lake <- get_C_lake(param,
                      start_date = as_datetime("2018-09-01"),
                      end_date = as_datetime("2018-10-31"),
                      summary_type = "median",
                      value_name = "C1_lake")
C2_lake <- get_C_lake(param,
                      start_date = as_datetime("2019-09-01"),
                      end_date = as_datetime("2019-10-31"),
                      summary_type = "median",
                      value_name = "C2_lake")
C_vals  <- left_join(C_pcpn, C1_lake, by = "lake") %>%
           left_join(C2_lake, by = "lake") %>%
           left_join(C_GWout, by = "lake")
C_GWin  <- calculate_C_GWin(C_vals)
C_vals  <- C_GWin %>%
           select(lake = .data$lake,
                  C_pcpn = .data$C_pcpn,
                  C0_lake = .data$C_GWout,
                  C_GWin = .data$C_GWin)

# 3. Calculate lake concentrations ---------------------------------------------
# No irr
MODFLOW_Mg <- list()
i <- 1
scenario <- "no_irr"
this_scenario <- CSLSdata::MODFLOW %>% filter(.data$scenario == !!scenario)
for (sim in unique(this_scenario$sim)) {
  if (sim %% 20 == 0) {
    message(sprintf("Starting sim %s", sim))
  }
  C_lake     <- calculate_C_lake(C_vals,
                                 scenario = scenario,
                                 sim = sim,
                                 dt = "day")
  MODFLOW_Mg[[i]] <- C_lake
  i <- i + 1
}
MODFLOW_Mg_no_irr <- bind_rows(MODFLOW_Mg)
save(MODFLOW_Mg_no_irr, file="data-raw/MODFLOW_Mg_no_irr.Rda")

# Cur irr
MODFLOW_Mg <- list()
i <- 1
scenario <- "cur_irr"
this_scenario <- CSLSdata::MODFLOW %>% filter(.data$scenario == !!scenario)
for (sim in unique(this_scenario$sim)) {
  if (sim %% 20 == 0) {
    message(sprintf("Starting sim %s", sim))
  }
  C_lake     <- calculate_C_lake(C_vals,
                                 scenario = scenario,
                                 sim = sim,
                                 dt = "day")
  MODFLOW_Mg[[i]] <- C_lake
  i <- i + 1
}
MODFLOW_Mg_cur_irr <- bind_rows(MODFLOW_Mg)
save(MODFLOW_Mg_cur_irr, file="data-raw/MODFLOW_Mg_cur_irr.Rda")

# Wells off
MODFLOW_Mg <- list()
i <- 1
scenario <- "wells_off"
this_scenario <- CSLSdata::MODFLOW %>% filter(.data$scenario == !!scenario)
for (sim in unique(this_scenario$sim)) {
  if (sim %% 20 == 0) {
    message(sprintf("Starting sim %s", sim))
  }
  C_lake     <- calculate_C_lake(C_vals,
                                 scenario = scenario,
                                 sim = sim,
                                 dt = "day")
  MODFLOW_Mg[[i]] <- C_lake
  i <- i + 1
}
MODFLOW_Mg_wells_off <- bind_rows(MODFLOW_Mg)
save(MODFLOW_Mg_wells_off, file="data-raw/MODFLOW_Mg_wells_off.Rda")

load("data-raw/MODFLOW_Mg_no_irr.Rda")
load("data-raw/MODFLOW_Mg_cur_irr.Rda")
load("data-raw/MODFLOW_Mg_wells_off.Rda")

# 4. Evaluate Metrics ----------------------------------------------------------
calculate_Mg_metrics <- function(df) {
  MODFLOW_Mg_metrics <- df %>%
                        group_by(.data$lake, .data$scenario, .data$sim) %>%
                        summarise(median = round(median(.data$C_lake),2),
                                  q10 = quantile(.data$C_lake,
                                                 probs = (1 - 0.1),
                                                 type = 6,
                                                 na.rm = TRUE),
                                  q90 = quantile(.data$C_lake,
                                                 probs = (1 - 0.9),
                                                 type = 6,
                                                 na.rm = TRUE),
                                  min = round(min(.data$C_lake),2),
                                  max = round(max(.data$C_lake),2),
                                  .groups = "drop") %>%
                       mutate(q10 = round(as.numeric(.data$q10),2),
                              q90 = round(as.numeric(.data$q90),2)) %>%
                       select(.data$lake, .data$scenario, .data$sim,
                              .data$min, .data$q90, .data$median, .data$q10,
                              .data$max) %>%
                       melt(id.vars = c("lake", "scenario", "sim")) %>%
                       mutate(metric = "solute_budget")
  return(MODFLOW_Mg_metrics)
}
MODFLOW_Mg_metrics_no_irr <- calculate_Mg_metrics(MODFLOW_Mg_no_irr)
MODFLOW_Mg_metrics_cur_irr <- calculate_Mg_metrics(MODFLOW_Mg_cur_irr)
MODFLOW_Mg_metrics_wells_off <- calculate_Mg_metrics(MODFLOW_Mg_wells_off)

MODFLOW_Mg_metrics <- bind_rows(MODFLOW_Mg_metrics_no_irr,
                             MODFLOW_Mg_metrics_cur_irr)
MODFLOW_Mg_metrics <- bind_rows(MODFLOW_Mg_metrics,
                                MODFLOW_Mg_metrics_wells_off)

MODFLOW_Mg1 <- MODFLOW_Mg_no_irr %>% filter(.data$sim %in% c(1,196,277))
MODFLOW_Mg2 <- MODFLOW_Mg_cur_irr %>% filter(.data$sim %in% c(1,196,277))
MODFLOW_Mg  <- bind_rows(MODFLOW_Mg1, MODFLOW_Mg2)
MODFLOW_Mg  <- bind_rows(MODFLOW_Mg, MODFLOW_Mg_wells_off)

# SAVE: Write out
# usethis::use_data(MODFLOW_Mg, MODFLOW_Mg_metrics, overwrite = TRUE, compress = "xz")
usethis::use_data(MODFLOW_Mg_metrics, overwrite = TRUE, compress = "xz")

