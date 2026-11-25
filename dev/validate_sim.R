# Validate fit_ve() against the known VE(t) in both HERMES scenarios.
# Run from the package root. Takes about 1-2 minutes; needs cmdstanr + CmdStan.
devtools::load_all(quiet = TRUE)
out_dir <- file.path(tempdir(), "surviver-validation")
dir.create(out_dir, showWarnings = FALSE)

cat("CmdStan:", cmdstanr::cmdstan_path(), as.character(cmdstanr::cmdstan_version()), "\n")

results <- list()
for (nm in c("ve_sim_waning", "ve_sim_no_waning")) {
  sim  <- get(nm)
  pop  <- ve_population(sim)
  inc  <- daily_incidence(pop)
  prev <- daily_prevalence(inc)
  foi  <- force_of_infection(prev, age_group_sizes(pop))
  pd   <- expand_person_days(counting_process(pop))
  sd   <- make_stan_data(pd, foi)
  x <- sd$cells
  cat(nm, ": events on floored-FoI days:", sum(x$Y[x$Ibar <= 1e-9]),
      "| first day Ibar>floor:", min(x$day[x$Ibar > 1e-9]),
      "| first event day:", min(x$day[x$Y > 0]), "\n")

  t0  <- Sys.time()
  fit <- fit_ve(sd, seed = 42, refresh = 0, show_messages = FALSE)
  cat(nm, ": fit took", format(Sys.time() - t0), "\n")
  print(fit$diagnostic_summary())
  print(fit$summary(c("log_beta", "sigma_gamma")))

  ve <- extract_ve(fit)
  ve$truth    <- true_ve(ve$days_since_vax, sim$params)
  ve$scenario <- nm
  results[[nm]] <- ve
  fit$save_object(file.path(out_dir, paste0("fit_", nm, ".rds")))
}
ve_all <- do.call(rbind, results)
saveRDS(ve_all, file.path(out_dir, "ve_all.rds"))

for (nm in names(results)) {
  v <- results[[nm]]
  v$band <- cut(v$days_since_vax, c(0, 7, 14, 30, 45, 60, 90, Inf), right = FALSE)
  cat("\n==", nm, "\n")
  print(aggregate(cbind(truth, median, lower, upper) ~ band, data = v,
                  FUN = function(z) round(mean(z), 2)))
  cat("truth inside 90% CrI:", round(mean(v$truth >= v$lower & v$truth <= v$upper), 2),
      "| max rhat:", round(max(v$rhat), 3), "| min ess_bulk:", round(min(v$ess_bulk)), "\n")
  print(plot_ve(v, truth = get(nm)$params) + ggplot2::ggtitle(nm))
}
cat("\nOutputs in", out_dir, "\n")
