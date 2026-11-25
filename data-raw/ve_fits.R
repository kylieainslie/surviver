## Precomputed fits of the standard analysis to ve_sim_waning and
## ve_sim_no_waning, stored in inst/extdata/ve_fits.rds for the "Get started"
## vignette (which cannot run Stan during R CMD check or on GitHub Actions).
##
## Rerun after changing the model, the datasets or the standard analysis:
##   source("data-raw/ve_fits.R")     # from the package root; ~2 minutes
## Needs cmdstanr and CmdStan.

devtools::load_all(quiet = TRUE)

censor_window <- 7

run_standard_analysis <- function(sim) {
  pop   <- ve_population(sim)
  inc   <- daily_incidence(pop)
  prev  <- daily_prevalence(inc, duration_infection = 7, reporting_delay = 2)
  foi   <- force_of_infection(prev, age_group_sizes(pop), contact_matrix = polymod_nl)
  early <- reported_soon_after_vaccination(pop, window = censor_window)
  cp    <- counting_process(pop, censor_at_vaccination = early)
  sd    <- make_stan_data(expand_person_days(cp), foi)
  fit   <- fit_ve(sd, seed = 42, refresh = 0, show_messages = FALSE)

  diag <- fit$diagnostic_summary()
  list(
    ve       = extract_ve(fit),
    expected = expected_cases(fit, sd, by = c("day", "vax_status")),
    diagnostics = list(
      max_rhat    = max(fit$summary(c("log_beta", "sigma_gamma", "VE"))$rhat, na.rm = TRUE),
      divergences = sum(diag$num_divergent),
      log_beta    = fit$summary("log_beta")$median,
      sigma_gamma = fit$summary("sigma_gamma")$median
    ),
    settings = list(
      model = "noncentered", contact_matrix = "polymod_nl", theta = 0,
      duration_infection = 7, reporting_delay = 2,
      censor_window = censor_window, n_censored = length(early),
      chains = 4, iter_warmup = 1000, iter_sampling = 1000, seed = 42
    )
  )
}

ve_fits <- list(
  waning    = run_standard_analysis(ve_sim_waning),
  no_waning = run_standard_analysis(ve_sim_no_waning),
  created   = Sys.Date(),
  versions  = list(
    surviver = as.character(utils::packageVersion("surviver")),
    cmdstan  = as.character(cmdstanr::cmdstan_version()),
    hermes_commit = ve_sim_waning$hermes_commit
  )
)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
saveRDS(ve_fits, "inst/extdata/ve_fits.rds", compress = "xz")
