# Lab notebook 2026-11-25: describe the packaged simulated datasets
# (ve_sim_waning, ve_sim_no_waning). Run from the package root; seconds.
devtools::load_all(quiet = TRUE)
out_dir <- Sys.getenv("SURVIVER_CHECKS_OUT", file.path(tempdir(), "surviver-checks"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

labels <- c(ve_sim_waning = "Waning", ve_sim_no_waning = "No waning")
sims <- list(ve_sim_waning = ve_sim_waning, ve_sim_no_waning = ve_sim_no_waning)

# true VE
ve_df <- do.call(rbind, lapply(names(sims), function(s) data.frame(
  scenario = labels[[s]], t = 0:180, ve = true_ve(0:180, sims[[s]]$params))))
p1 <- ggplot2::ggplot(ve_df, ggplot2::aes(t, ve, colour = scenario)) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::labs(x = "Days since vaccination", y = "True VE", colour = NULL) +
  ggplot2::theme_minimal() + ggplot2::theme(legend.position = "top")
ggplot2::ggsave(file.path(out_dir, "true_ve.png"), p1, width = 6, height = 3.5, dpi = 150)

# weekly reported cases and cumulative vaccination coverage
epi <- do.call(rbind, lapply(names(sims), function(s) {
  wk <- floor(sims[[s]]$incidence$time_of_reporting / 7) * 7
  tab <- as.data.frame(table(week = wk), stringsAsFactors = FALSE)
  data.frame(scenario = labels[[s]], week = as.numeric(tab$week), cases = tab$Freq)
}))
cov <- do.call(rbind, lapply(names(sims), function(s) {
  v <- sims[[s]]$vaccination$time_of_vaccination
  data.frame(scenario = labels[[s]], day = 0:180,
             coverage = sapply(0:180, function(d) sum(v <= d)) / nrow(sims[[s]]$demographic))
}))
k <- max(epi$cases) / 0.5
p2 <- ggplot2::ggplot() +
  ggplot2::geom_col(data = epi, ggplot2::aes(week + 3.5, cases), fill = "grey60", width = 6.5) +
  ggplot2::geom_step(data = cov, ggplot2::aes(day, coverage * k), colour = "#1f78b4", linewidth = 0.9) +
  ggplot2::facet_wrap(~ scenario, ncol = 1) +
  ggplot2::scale_y_continuous("Reported cases per week",
    sec.axis = ggplot2::sec_axis(~ . / k, name = "Vaccination coverage",
                                 labels = function(x) paste0(100 * x, "%"))) +
  ggplot2::labs(x = "Study day") + ggplot2::theme_minimal()
ggplot2::ggsave(file.path(out_dir, "epidemic.png"), p2, width = 7, height = 5, dpi = 150)

# vaccinated cases by days since vaccination (reported after vaccination, by day 180)
for (s in names(sims)) {
  pop <- ve_population(sims[[s]])
  v <- pop[pop$vaccinated & pop$infected & pop$time_of_reporting > pop$time_of_vaccination &
             pop$time_of_reporting <= sims[[s]]$study_end, ]
  cat(labels[[s]], ": vaccinated cases by days since vaccination\n")
  print(table(cut(v$time_of_reporting - v$time_of_vaccination,
                  c(0, 7, 14, 30, 45, 60, 90, Inf), right = FALSE)))
}
cat("Figures in", out_dir, "\n")
