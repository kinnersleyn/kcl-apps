# Builds all five apps into docs/ for GitHub Pages.
# Run from the repository root (the folder that contains KCL-Interims, KCL-Hazards, ...):
#   setwd("path/to/kcl-apps"); source("build_all.R")
library(shinylive)

apps <- c(Interims = "KCL-Interims",
          Hazards  = "KCL-Hazards",
          Sim_H0   = "KCL-Sim_H0",
          T12E     = "KCL-T12E",
          Bayesian = "KCL-Bayesian")

for (nm in names(apps)) {
  dest <- file.path("docs", nm)
  unlink(dest, recursive = TRUE)
  shinylive::export(file.path(apps[[nm]], "app"), dest)
}
message("Done. Commit and push, then wait a minute or two for GitHub Pages to update.")
