# Lab notebook entry 10 (2026-11-25): FoI built from the HERMES contact structure
# instead of POLYMOD. Run from the package root; needs cmdstanr; ~5 minutes.
devtools::load_all(quiet = TRUE)
out_dir <- Sys.getenv("SURVIVER_CHECKS_OUT", file.path(tempdir(), "surviver-checks"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

breaks <- c(0, 7, 14, 30, 45, 60, 90, 120, Inf)
labels <- c("1-6", "7-13", "14-29", "30-44", "45-59", "60-89", "90-119", "120+")

stan_res <- list(); glm_res <- list(); fit_res <- list()
for (nm in c("ve_sim_waning", "ve_sim_no_waning")) {
  sim  <- get(nm)
  pop  <- ve_population(sim)
  inc  <- daily_incidence(pop)
  prev <- daily_prevalence(inc)
  n    <- age_group_sizes(pop)
  pd   <- expand_person_days(counting_process(pop))

  fois <- list(
    POLYMOD = force_of_infection(prev, n, polymod_nl),
    HERMES  = force_of_infection(prev, n, contact_matrix_hermes(sim$demographic, sim$params)),
    overall = force_of_infection(prev, n, age_structured = FALSE)
  )

  for (f in names(fois)) {
    sd <- make_stan_data(pd, fois[[f]])
    cells <- sd$cells

    # 1. how well does the FoI offset track the unvaccinated hazard?
    u <- cells[cells$vax_status == 0, ]
    m0 <- glm(Y ~ 1, offset = log(Ibar) + log(n_persdays), family = poisson, data = u)
    m_age <- glm(Y ~ age_group, offset = log(Ibar) + log(n_persdays), family = poisson, data = u)
    fit_res[[paste(nm, f)]] <- data.frame(
      scenario = nm, foi = f,
      deviance_const = round(deviance(m0)), deviance_age = round(deviance(m_age)),
      age_rr_range = paste(round(range(exp(c(0, coef(m_age)[-1]))), 2), collapse = "-"))

    # 2. Stan recovery
    fit <- fit_ve(sd, seed = 42, refresh = 0, show_messages = FALSE)
    ve <- extract_ve(fit)
    ve$truth <- true_ve(ve$days_since_vax, sim$params)
    ve$scenario <- nm; ve$foi <- f
    ve$divergent <- sum(fit$diagnostic_summary()$num_divergent)
    stan_res[[paste(nm, f)]] <- ve

    # 3. piecewise GLMs: constant baseline (A) vs day baseline (B) with this FoI
    cells$band <- factor(ifelse(cells$vax_status == 1,
      as.character(cut(cells$tvax_bin - 1L, breaks, labels, right = FALSE)), "unvax"),
      levels = c("unvax", labels))
    cd <- cells[cells$day %in% unique(cells$day[cells$Y > 0]), ]
    cd$fday <- factor(cd$day)
    for (m in c("A", "B")) {
      g <- if (m == "A") {
        glm(Y ~ band, offset = log(Ibar) + log(n_persdays), family = poisson, data = cells)
      } else {
        glm(Y ~ fday + band, offset = log(Ibar) + log(n_persdays), family = poisson, data = cd)
      }
      co <- summary(g)$coefficients[paste0("band", labels), , drop = FALSE]
      vc <- cells[cells$vax_status == 1, ]
      vc$true <- true_ve(vc$tvax_bin - 1L, sim$params)
      truth <- tapply(vc$true * vc$n_persdays, vc$band, sum) / tapply(vc$n_persdays, vc$band, sum)
      ev <- tapply(vc$Y, vc$band, sum)
      glm_res[[paste(nm, f, m)]] <- data.frame(
        scenario = nm, foi = f, model = m, band = labels,
        events = as.integer(ev[labels]), truth = truth[labels],
        est = 1 - exp(co[, 1]), se = co[, 2])
    }
  }
}
stan_all <- do.call(rbind, stan_res); glm_all <- do.call(rbind, glm_res)
saveRDS(list(stan = stan_all, glm = glm_all, fit = do.call(rbind, fit_res)),
        file.path(out_dir, "hermes_contacts.rds"))

cat("\n== 1. FoI fit to unvaccinated cases (Poisson deviance; lower = better)\n")
print(do.call(rbind, fit_res), row.names = FALSE)

cat("\n== 2. Stan recovery\n")
for (k in names(stan_res)) {
  v <- stan_res[[k]]; late <- v$days_since_vax >= 90; mid <- v$days_since_vax %in% 14:89
  cat(sprintf("%-28s coverage %.2f | MAE days 1-120 %.3f | error 14-89 %+.3f | error 90+ %+.3f | mean CrI width %.2f | max rhat %.3f | div %d\n",
      k, mean(v$truth >= v$lower & v$truth <= v$upper),
      mean(abs(v$median - v$truth)[v$days_since_vax <= 120]),
      mean((v$median - v$truth)[mid]), mean((v$median - v$truth)[late]),
      mean((v$upper - v$lower)[v$days_since_vax <= 90]), max(v$rhat), v$divergent[1]))
}

cat("\n== 3. Piecewise GLM: case-weighted error from day 90 (and mean SE of log HR, days 14-89)\n")
for (nm in unique(glm_all$scenario)) for (f in unique(glm_all$foi)) for (m in c("A", "B")) {
  x <- glm_all[glm_all$scenario == nm & glm_all$foi == f & glm_all$model == m, ]
  late <- x$band %in% c("90-119", "120+") & x$events > 0
  cat(sprintf("%-18s %-8s %s: late error %+.2f | mean SE 14-89 %.3f\n", nm, f, m,
              weighted.mean(x$est[late] - x$truth[late], x$events[late]),
              mean(x$se[x$band %in% c("14-29", "30-44", "45-59", "60-89")])))
}

p <- ggplot2::ggplot(stan_all, ggplot2::aes(days_since_vax, median, colour = foi, fill = foi)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.1, colour = NA) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_line(ggplot2::aes(y = truth), colour = "black", linetype = "dashed", linewidth = 0.6) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
  ggplot2::facet_wrap(~ scenario) + ggplot2::coord_cartesian(ylim = c(-0.6, 1)) +
  ggplot2::labs(x = "Days since vaccination", y = "VE", colour = "FoI", fill = "FoI",
                title = "VE(t) with POLYMOD, HERMES and overall-prevalence FoI (black dashed = true VE)") +
  ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "hermes_contacts.png"), p, width = 10, height = 4.5, dpi = 150)
cat("\nOutputs in", out_dir, "\n")
