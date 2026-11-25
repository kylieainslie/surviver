## Simulated ground-truth datasets for time-varying VE: ve_sim_waning and
## ve_sim_no_waning.
##
## Generated with the HERMES individual-based model, waning-VE extension:
##   https://github.com/kylieainslie/ibm_hermes  (branch waning-vaccine-effectiveness)
##
## Requires a local clone of ibm_hermes. Point HERMES_PATH at it if it is not
## at ~/Documents/ibm_hermes, e.g. in ~/.Renviron:
##   HERMES_PATH=/path/to/ibm_hermes
##
## The HERMES commit used is stored in each dataset (`hermes_commit`).
## Runtime: about 20 seconds per scenario.

hermes_path <- normalizePath(Sys.getenv("HERMES_PATH", "~/Documents/ibm_hermes"),
                             mustWork = TRUE)

hermes_commit <- system2("git", c("-C", shQuote(hermes_path), "rev-parse", "HEAD"),
                         stdout = TRUE)
hermes_dirty  <- length(system2("git", c("-C", shQuote(hermes_path), "status",
                                         "--porcelain", "--", "lib"),
                                stdout = TRUE)) > 0
if (hermes_dirty) {
  warning("ibm_hermes/lib has uncommitted changes; the recorded commit (",
          hermes_commit, ") will not exactly describe the code used.")
}

# HERMES scripts source each other with paths relative to the repo root, so
# load them from there, into their own environment
hermes <- new.env()
withr::with_dir(hermes_path, {
  sys.source("lib/ibm_core.R", envir = hermes)
  sys.source("lib/ibm_waning_scenario.R", envir = hermes)
})

make_scenario <- function(ve_floor) {
  sc <- hermes$simulate_waning_scenario(ve_floor = ve_floor)
  # drop HERMES's per-snapshot `time` column (constant within each table)
  list(
    demographic   = sc$demographic_data[setdiff(names(sc$demographic_data), "time")],
    vaccination   = sc$vaccination_data[setdiff(names(sc$vaccination_data), "time")],
    incidence     = sc$incidence_data_reduced[setdiff(names(sc$incidence_data_reduced), "time")],
    params        = sc$params[setdiff(names(sc$params), "output_dir")],
    study_end     = sc$study_end,
    hermes_commit = hermes_commit
  )
}

ve_sim_waning    <- make_scenario(ve_floor = 0.1)
ve_sim_no_waning <- make_scenario(ve_floor = 0.8)   # = peak VE, so no waning

usethis::use_data(ve_sim_waning, ve_sim_no_waning, overwrite = TRUE)
