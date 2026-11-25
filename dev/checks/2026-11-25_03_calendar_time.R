# Lab notebook 2026-11-25 #3. Run from the package root; base R glm only, runs in seconds.
# Is the late VE drop caused by calendar time leaking into time since
# vaccination? Compare piecewise-constant VE from Poisson GLMs on the same
# cells as the Stan model, with and without per-day baselines.
devtools::load_all(quiet = TRUE)
out_dir <- Sys.getenv("SURVIVER_CHECKS_OUT", file.path(tempdir(), "surviver-checks"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

breaks <- c(0, 7, 14, 30, 45, 60, 90, 120, Inf)
labels <- c("1-6", "7-13", "14-29", "30-44", "45-59", "60-89", "90-119", "120+")

all_res <- list()
for (nm in c("ve_sim_waning", "ve_sim_no_waning")) {
  sim  <- get(nm)
  pop  <- ve_population(sim)
  inc  <- daily_incidence(pop)
  foi  <- force_of_infection(daily_prevalence(inc), age_group_sizes(pop))
  cells <- make_stan_data(expand_person_days(counting_process(pop)), foi)$cells

  cells$dsv  <- cells$tvax_bin - 1L
  cells$band <- factor(ifelse(cells$vax_status == 1,
                              as.character(cut(cells$dsv, breaks, labels, right = FALSE)),
                              "unvax"),
                       levels = c("unvax", labels))
  cells$fday <- factor(cells$day)
  # days with no events at all carry no information once each day has its
  # own baseline; drop them to avoid infinite day coefficients
  event_days <- unique(cells$day[cells$Y > 0])
  cd <- cells[cells$day %in% event_days, ]

  fits <- list(
    "A: constant baseline + FoI offset" =
      glm(Y ~ band, offset = log(Ibar) + log(n_persdays), family = poisson, data = cells),
    "B: day baseline + FoI offset" =
      glm(Y ~ fday + band, offset = log(Ibar) + log(n_persdays), family = poisson, data = cd),
    "C: day baseline + age, no FoI (Cox-like)" =
      glm(Y ~ fday + age_group + band, offset = log(n_persdays), family = poisson, data = cd)
  )

  # truth per band: vaccinated person-day weighted mean of true VE
  vc <- cells[cells$vax_status == 1, ]
  vc$true <- true_ve(vc$dsv, sim$params)
  truth <- tapply(vc$true * vc$n_persdays, vc$band, sum) /
    tapply(vc$n_persdays, vc$band, sum)
  events <- tapply(vc$Y, vc$band, sum)

  res <- data.frame(band = labels, events = as.integer(events[labels]),
                    truth = round(truth[labels], 2))
  for (f in names(fits)) {
    co <- summary(fits[[f]])$coefficients
    rows <- paste0("band", labels)
    est <- co[rows, "Estimate"]; se <- co[rows, "Std. Error"]
    res[[f]] <- sprintf("%.2f (%.2f, %.2f)", 1 - exp(est),
                        1 - exp(est + 1.645 * se), 1 - exp(est - 1.645 * se))
    res[[paste0("pt_", substr(f, 1, 1))]] <- 1 - exp(est)
    res[[paste0("lo_", substr(f, 1, 1))]] <- 1 - exp(est + 1.645 * se)
    res[[paste0("hi_", substr(f, 1, 1))]] <- 1 - exp(est - 1.645 * se)
  }
  cat("\n==", nm, "  [VE (90% CI) per band]\n")
  print(res[, c("band", "events", "truth", names(fits))], row.names = FALSE)
  late <- res$band %in% c("90-119", "120+")
  w <- res$events[late]
  cat("Late (90+) event-weighted error, estimate - truth:  A",
      round(weighted.mean(res$pt_A[late] - res$truth[late], w), 2),
      "| B", round(weighted.mean(res$pt_B[late] - res$truth[late], w), 2),
      "| C", round(weighted.mean(res$pt_C[late] - res$truth[late], w), 2), "\n")
  res$scenario <- nm
  all_res[[nm]] <- res
}
saveRDS(all_res, file.path(out_dir, "calendar_day.rds"))

# figure: band estimates by model, waning scenario
w <- all_res[["ve_sim_waning"]]
long <- do.call(rbind, lapply(c("A", "B", "C"), function(m) data.frame(
  band = factor(w$band, levels = labels), truth = w$truth, model = m,
  est = w[[paste0("pt_", m)]], lo = w[[paste0("lo_", m)]], hi = w[[paste0("hi_", m)]])))
long$model <- factor(long$model, labels = c("A: constant baseline + FoI",
                                            "B: day baseline + FoI",
                                            "C: day baseline + age (Cox-like)"))
p <- ggplot2::ggplot(long, ggplot2::aes(band, est, colour = model)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
  ggplot2::geom_pointrange(ggplot2::aes(ymin = pmax(lo, -1), ymax = hi),
                           position = ggplot2::position_dodge(width = 0.6)) +
  ggplot2::geom_point(ggplot2::aes(band, truth), colour = "black", shape = 4, size = 3,
                      inherit.aes = FALSE, data = w |> transform(band = factor(band, levels = labels))) +
  ggplot2::coord_cartesian(ylim = c(-1, 1)) +
  ggplot2::labs(x = "Days since vaccination", y = "VE (90% CI)", colour = NULL,
                title = "Waning scenario: effect of a per-day baseline (x = true VE)") +
  ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "calendar_time.png"), p, width = 9, height = 4.5, dpi = 150)
cat("Figure:", file.path(out_dir, "calendar_time.png"), "\n")
