# Lab notebook 2026-11-26 #11: HERMES vaccinates infected and recovered
# people, and vaccination resets their state to V (susceptible with vaccine
# protection). How often does this happen, and how many cases among
# vaccinated people does it produce?
# Run from the package root; needs a local ibm_hermes clone (HERMES_PATH);
# about 1 minute.
devtools::load_all(quiet = TRUE)
hermes_path <- normalizePath(Sys.getenv("HERMES_PATH", "~/Documents/ibm_hermes"), mustWork = TRUE)
hermes <- new.env()
withr::with_dir(hermes_path, {
  sys.source("lib/ibm_core.R", envir = hermes)
  sys.source("lib/ibm_waning_scenario.R", envir = hermes)
})

for (nm in c("ve_sim_waning", "ve_sim_no_waning")) {
  sim <- get(nm)
  od <- file.path(tempdir(), nm)
  sc <- hermes$simulate_waning_scenario(ve_floor = sim$params$ve_floor, output_dir = od)
  stopifnot(identical(sc$incidence_data_reduced$time_of_reporting, sim$incidence$time_of_reporting))

  # state of each person just before the round in which they were vaccinated
  prev <- readRDS(file.path(od, "pop_data.rds"))
  state_at_vax <- character(0)
  for (i in seq_len(16)) {
    cur <- readRDS(file.path(od, paste0("pop_data_intervention", i, ".rds")))
    newly <- which(!is.na(cur$time_of_vaccination) & is.na(prev$time_of_vaccination))
    state_at_vax[rownames(prev)[newly]] <- prev$health[newly]
    prev <- cur
  }
  fin <- prev
  was_R <- names(state_at_vax)[state_at_vax == "R"]
  reinf <- !is.na(fin[was_R, "time_of_infection"]) &
    fin[was_R, "time_of_infection"] > fin[was_R, "time_of_vaccination"]

  inc <- sim$incidence
  tv  <- fin[inc$ind_id, "time_of_vaccination"]
  vcase <- !is.na(tv) & inc$time_of_reporting > tv & inc$time_of_reporting <= sim$study_end

  cat("\n==", nm, "\n")
  print(table(state_at_vaccination = state_at_vax))
  cat("Recovered when vaccinated:", length(was_R),
      "| reinfected after vaccination:", sum(reinf), "\n")
  cat("Cases among vaccinated (reported after vaccination, by day", sim$study_end, "):",
      sum(vcase), "| from recovered-then-vaccinated people:",
      sum(inc$ind_id[vcase] %in% was_R), "\n")
}
