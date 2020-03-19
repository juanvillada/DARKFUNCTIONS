.onAttach <- function(libname, pkgname) {
  # Startup message
  packageStartupMessage("Auxiliary functions loaded")

  # Create vectors of environmental variables
  CAT.VARS <- c("Sample", "Layer", "Ecosystem", "Vegetation")
  NUM.VARS <- c("SoilWater", "pH", "SOM")
  PROXY.VARS <- c("SnowDepth", "Elevation")
  FLUX.VARS <- c("CH4", "CO2", "N2O")

  assign("CAT.VARS", CAT.VARS, .GlobalEnv)
  assign("NUM.VARS", NUM.VARS, .GlobalEnv)
  assign("PROXY.VARS", PROXY.VARS, .GlobalEnv)
  assign("FLUX.VARS", FLUX.VARS, .GlobalEnv)
}

.onLoad <- function(libname, pkgname) {
  # Create vectors of environmental variables
  CAT.VARS <- c("Sample", "Layer", "Ecosystem", "Vegetation")
  NUM.VARS <- c("SoilWater", "pH", "SOM")
  PROXY.VARS <- c("SnowDepth", "Elevation")
  FLUX.VARS <- c("CH4", "CO2", "N2O")

  assign("CAT.VARS", CAT.VARS, .GlobalEnv)
  assign("NUM.VARS", NUM.VARS, .GlobalEnv)
  assign("PROXY.VARS", PROXY.VARS, .GlobalEnv)
  assign("FLUX.VARS", FLUX.VARS, .GlobalEnv)
}
