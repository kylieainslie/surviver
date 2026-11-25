# Lab notebook 2026-11-26 #12: censoring people infected before vaccination.
# Compares no censoring, an oracle rule (true infection status from HERMES)
# and observable window rules (reported within k days of vaccination).
# Run from the package root; needs cmdstanr and a local ibm_hermes clone
# (HERMES_PATH); about 10 minutes.
devtools::load_all(quiet = TRUE)
out_dir <- Sys.getenv("SURVIVER_CHECKS_OUT", file.path(tempdir(), "surviver-checks"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

hermes_path <- normalizePath(Sys.getenv("HERMES_PATH", "~/Documents/ibm_hermes"), mustWork = TRUE)
hermes <- new.env()
withr::with_dir(hermes_path, {
  sys.source("lib/ibm_core.R", envir = hermes)
  sys.source("lib/ibm_waning_scenario.R", envir = hermes)
})

breaks <- c(0, 7, 14, 30, 45, 60, 90, 120, Inf)
labels <- c("1-6", "7-13", "14-29", "30-44", "45-59", "60-89", "90-119", "120+")

stan_res <- list(); glm_res <- list(); rule_res <- list()
for (nm in c("ve_sim_waning", "ve_sim_no_waning")) {
  sim <- get(nm)
  od  <- file.path(tempdir(), nm)
  sc  <- hermes$simulate_waning_scenario(ve_floor = sim$params$ve_floor, output_dir = od)
  stopifnot(identical(sc$incidence_data_reduced$time_of_reporting, sim$incidence$time_of_reporting))

  # true state just before each person's vaccination round
  prev <- readRDS(file.path(od, "pop_data.rds"))
  state_at_vax <- character(0)
  for (i in seq_len(16)) {
    cur <- readRDS(file.path(od, paste0("pop_data_intervention", i, ".rds")))
    newly <- which(!is.na(cur$time_of_vaccination) & is.na(prev$time_of_vaccination))
    state_at_vax[rownames(prev)[newly]] <- prev$health[newly]
    prev <- cur
  }
  infected_before <- names(state_at_vax)[state_at_vax %in% c("I", "R")]

  pop  <- ve_population(sim)
  se   <- sim$study_end
  foi  <- force_of_infection(daily_prevalence(daily_incidence(pop)), age_group_sizes(pop))
  # people whose censoring would change anything: vaccinated, reported after
  # vaccination and by study end
  affected <- as.character(pop$ind_id[pop$vaccinated & pop$infected &
    pop$time_of_reporting > pop$time_of_vaccination & pop$time_of_reporting <= se])

  rules <- list(
    none      = character(0),
    oracle    = infected_before,
    window_7  = reported_soon_after_vaccination(pop, 7),
    window_14 = reported_soon_after_vaccination(pop, 14)
  )

  for (r in names(rules)) {
    ids <- intersect(rules[[r]], affected)
    rule_res[[paste(nm, r)]] <- data.frame(
      scenario = nm, rule = r, censored_cases = length(ids),
      truly_infected_before = sum(ids %in% infected_before),
      infected_before_missed = sum(!(intersect(infected_before, affected) %in% ids)))

    cp <- counting_process(pop, censor_at_vaccination = rules[[r]])
    sd <- make_stan_data(expand_person_days(cp), foi)

    fit <- fit_ve(sd, seed = 42, refresh = 0, show_messages = FALSE)
    ve  <- extract_ve(fit)
    ve$truth <- true_ve(ve$days_since_vax, sim$params)
    ve$scenario <- nm; ve$rule <- r
    stan_res[[paste(nm, r)]] <- ve

    cells <- sd$cells
    cells$band <- factor(ifelse(cells$vax_status == 1,
      as.character(cut(cells$tvax_bin - 1L, breaks, labels, right = FALSE)), "unvax"),
      levels = c("unvax", labels))
    cd <- cells[cells$day %in% unique(cells$day[cells$Y > 0]), ]
    cd$fday <- factor(cd$day)
    g  <- glm(Y ~ fday + band, offset = log(Ibar) + log(n_persdays), family = poisson, data = cd)
    co <- summary(g)$coefficients
    rows <- paste0("band", labels)
    vc <- cells[cells$vax_status == 1, ]
    vc$true <- true_ve(vc$tvax_bin - 1L, sim$params)
    truth <- tapply(vc$true * vc$n_persdays, vc$band, sum) / tapply(vc$n_persdays, vc$band, sum)
    glm_res[[paste(nm, r)]] <- data.frame(
      scenario = nm, rule = r, band = labels,
      events = as.integer(tapply(vc$Y, vc$band, sum)[labels]),
      truth = truth[labels],
      est = ifelse(rows %in% rownames(co), 1 - exp(co[match(rows, rownames(co)), 1]), NA))
  }
}
stan_all <- do.call(rbind, stan_res); glm_all <- do.call(rbind, glm_res)
rule_all <- do.call(rbind, rule_res)
saveRDS(list(stan = stan_all, glm = glm_all, rules = rule_all),
        file.path(out_dir, "censor_prior_infection.rds"))

cat("\n== Who each rule censors (among vaccinated people reported after vaccination)\n")
print(rule_all, row.names = FALSE)

cat("\n== Stan: mean error (median - truth) by range of days since vaccination, and coverage\n")
for (k in names(stan_res)) {
  v <- stan_res[[k]]
  v$band <- cut(v$days_since_vax, breaks, labels, right = FALSE)
  err <- tapply(v$median - v$truth, v$band, mean)
  cat(sprintf("%-28s %s | coverage %.2f\n", k,
              paste(sprintf("%s:%+.2f", labels, err[labels]), collapse = " "),
              mean(v$truth >= v$lower & v$truth <= v$upper)))
}

cat("\n== GLM with per-day baseline (model B): VE by range, truth in brackets\n")
for (k in names(glm_res)) {
  x <- glm_res[[k]]
  cat(sprintf("%-28s %s\n", k, paste(sprintf("%s:%.2f[%.2f](n=%d)", x$band, x$est, x$truth, x$events), collapse = " ")))
}

p <- ggplot2::ggplot(stan_all, ggplot2::aes(days_since_vax, median, colour = rule)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_line(ggplot2::aes(y = truth), colour = "black", linetype = "dashed", linewidth = 0.6) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
  ggplot2::facet_wrap(~ scenario) + ggplot2::coord_cartesian(xlim = c(0, 60), ylim = c(-0.4, 1)) +
  ggplot2::labs(x = "Days since vaccination", y = "VE (posterior median)", colour = "Censoring",
                title = "Censoring people infected before vaccination (first 60 days; black dashed = true VE)") +
  ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "censor_prior_infection.png"), p, width = 10, height = 4.5, dpi = 150)
cat("\nOutputs in", out_dir, "\n")
