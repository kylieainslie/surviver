# Lab notebook entry 8 (2026-11-25). Run from the package root; needs cmdstanr and a local ibm_hermes clone (HERMES_PATH).
# Does removing unreported (but immune) infections from the risk set fix the
# late downward bias in VE(t)? Same model, same events (reported cases at
# report day), same FoI; only change: people infected but never reported
# are censored at their true infection time.
devtools::load_all(quiet = TRUE)
out_dir <- Sys.getenv("SURVIVER_CHECKS_OUT", file.path(tempdir(), "surviver-checks"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

hermes_path <- normalizePath(Sys.getenv("HERMES_PATH", "~/Documents/ibm_hermes"), mustWork = TRUE)
hermes <- new.env()
withr::with_dir(hermes_path, {
  sys.source("lib/ibm_core.R", envir = hermes)
  sys.source("lib/ibm_waning_scenario.R", envir = hermes)
})

results <- list()
for (nm in c("ve_sim_waning", "ve_sim_no_waning")) {
  sim <- get(nm)

  # regenerate with HERMES to recover the full individual-level state
  run_dir <- file.path(tempdir(), nm)
  sc <- hermes$simulate_waning_scenario(ve_floor = sim$params$ve_floor,
                                        output_dir = run_dir)
  full <- readRDS(file.path(run_dir, "pop_data_intervention16.rds"))
  full$ind_id <- rownames(full)

  # the regenerated run must match the packaged data exactly
  stopifnot(
    identical(sc$vaccination_data$ind_id, sim$vaccination$ind_id),
    identical(sc$incidence_data_reduced$ind_id, sim$incidence$ind_id),
    identical(sc$incidence_data_reduced$time_of_reporting, sim$incidence$time_of_reporting)
  )

  pop <- ve_population(sim)
  se  <- attr(pop, "study_end")
  pop$time_of_infection <- full$time_of_infection[match(pop$ind_id, full$ind_id)]
  unrep <- !is.na(pop$time_of_infection) & pop$time_of_infection <= se &
    !(pop$infected & pop$time_of_reporting <= se)
  cat("\n==", nm, "\n")
  cat("ever infected by study end:", sum(!is.na(pop$time_of_infection) & pop$time_of_infection <= se),
      "| reported by study end:", sum(pop$infected & pop$time_of_reporting <= se),
      "| unreported (censored at infection):", sum(unrep), "\n")
  cat("unreported by vaccination status at infection: vaccinated",
      sum(unrep & pop$vaccinated & pop$time_of_vaccination < pop$time_of_infection),
      "| unvaccinated", sum(unrep & !(pop$vaccinated & pop$time_of_vaccination < pop$time_of_infection)), "\n")

  inc  <- daily_incidence(pop)
  foi  <- force_of_infection(daily_prevalence(inc), age_group_sizes(pop))
  cp   <- counting_process(pop)

  # censor unreported infections at their infection time
  t_inf <- pop$time_of_infection[match(cp$ind_id, pop$ind_id)]
  is_unrep <- unrep[match(cp$ind_id, pop$ind_id)]
  cut_here <- is_unrep & t_inf < cp$tstop
  cp_true <- cp[!(cut_here & cp$tstart >= t_inf), ]
  t_inf2 <- pop$time_of_infection[match(cp_true$ind_id, pop$ind_id)]
  trunc  <- unrep[match(cp_true$ind_id, pop$ind_id)] & t_inf2 < cp_true$tstop
  cp_true$tstop[trunc] <- t_inf2[trunc]
  cp_true$event[trunc] <- 0
  cp_true <- cp_true[cp_true$tstop > cp_true$tstart, ]
  stopifnot(sum(cp_true$event) == sum(cp$event))

  for (version in c("reported_only", "true_risk_set")) {
    this_cp <- if (version == "reported_only") cp else cp_true
    sd  <- make_stan_data(expand_person_days(this_cp), foi)
    fit <- fit_ve(sd, seed = 42, refresh = 0, show_messages = FALSE)
    d <- fit$diagnostic_summary()
    cat(version, ": person-days", sum(exp(sd$data$log_persdays)),
        "| divergences", sum(d$num_divergent),
        "| max rhat", round(max(fit$summary("VE")$rhat, na.rm = TRUE), 3), "\n")
    ve <- extract_ve(fit)
    ve$truth    <- true_ve(ve$days_since_vax, sim$params)
    ve$scenario <- nm
    ve$version  <- version
    results[[paste(nm, version)]] <- ve
  }
}
ve_all <- do.call(rbind, results)
saveRDS(ve_all, file.path(out_dir, "ve_true_risk_set.rds"))

cat("\nMean VE by band (median [90% CrI]); truth = true VE\n")
ve_all$band <- cut(ve_all$days_since_vax, c(0, 7, 14, 30, 45, 60, 90, 120, Inf), right = FALSE)
for (nm in unique(ve_all$scenario)) {
  cat("\n==", nm, "\n")
  x <- ve_all[ve_all$scenario == nm, ]
  tab <- aggregate(cbind(truth, median, lower, upper) ~ band + version, data = x,
                   FUN = function(z) round(mean(z), 2))
  print(reshape(tab, idvar = "band", timevar = "version", direction = "wide"))
  for (v in unique(x$version)) {
    y <- x[x$version == v & x$days_since_vax >= 90, ]
    cat(v, ": days 90+ mean error (median - truth):", round(mean(y$median - y$truth), 3),
        "| truth in 90% CrI (all days):",
        round(mean(x$truth[x$version == v] >= x$lower[x$version == v] &
                     x$truth[x$version == v] <= x$upper[x$version == v]), 2), "\n")
  }
}

p <- ggplot2::ggplot(ve_all, ggplot2::aes(days_since_vax, median, colour = version, fill = version)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_line(ggplot2::aes(y = truth), colour = "black", linetype = "dashed", linewidth = 0.6) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
  ggplot2::facet_wrap(~ scenario) +
  ggplot2::coord_cartesian(ylim = c(-0.6, 1)) +
  ggplot2::labs(x = "Days since vaccination", y = "VE", colour = NULL, fill = NULL,
                title = "Reported-only risk set vs true risk set (black dashed = true VE)") +
  ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "true_risk_set.png"), p, width = 10, height = 4.5, dpi = 150)
