# CSLSfluxes

CSLSfluxes is an R package containing code and results generated by the
Wisconsin Department of Natural Resources related to lake water and solute
budgets in the Central Sands Lakes Study (CSLS).

The key dataset included here contains summary metrics of the solute budget for
each 38-year scenario (`data/MODFLOW_Mg_metrics`). 

The full datasets of monthly lake solute estimates are quite large and a bit 
unwieldy to include here, but can be recreated via the script
`data-raw/calculate_MODFLOW_Mg.R`.

Details on data and code can also be found in `inst/CSLSfluxes_1.0.0.pdf`.


## Installation

Several data files in this package are hosted via Git Large File Storage (LFS).
You may need to install LFS in order for some files to download properly via
git.

To install or explore this R package, you may:

  1. Use devtools to install the package. This will allow you to load the 
  cleaned Rda data files and explore vignettes.
  ```
  devtools::install_github("WDNR-Water-Use/CSLSfluxes", build_vignettes=T)
  ```

  2. Fork from github, clone or download ZIP, then open the CSLSfluxes.prj file in 
  R. This will allow you to also explore the raw data files and cleaning scripts 
  included in data-raw/.

