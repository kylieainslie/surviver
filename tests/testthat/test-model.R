prev_toy <- function() {
  # two age groups, 3 days; group "a" has 2 unvaccinated infectious on day 1,
  # group "b" has 1 vaccinated infectious on day 2
  grid <- expand.grid(day = 0:3, age_group = c("a", "b"), vax_status = 0:1,
                      stringsAsFactors = FALSE)
  grid$age_group <- factor(grid$age_group, levels = c("a", "b"))
  grid$n_infectious <- 0L
  grid$n_infectious[grid$day == 1 & grid$age_group == "a" & grid$vax_status == 0] <- 2L
  grid$n_infectious[grid$day == 2 & grid$age_group == "b" & grid$vax_status == 1] <- 1L
  grid
}
cm_toy <- matrix(c(1, 2, 3, 4), 2, dimnames = list(c("a", "b"), c("a", "b")))
n_toy  <- c(a = 10, b = 20)

test_that("force_of_infection matches the formula with a one-day lag", {
  foi <- force_of_infection(prev_toy(), n_toy, cm_toy, theta = 0)
  get <- function(d, g) foi$Ibar[foi$day == d & foi$age_group == g]
  # day 2 uses day-1 prevalence: I_a = 2
  expect_equal(get(2, "a"), cm_toy["a", "a"] * 2 / 10)
  expect_equal(get(2, "b"), cm_toy["a", "b"] * 2 / 20)
  # day 3 uses day-2 prevalence: I_b = 1 (vaccinated, theta = 0)
  expect_equal(get(3, "a"), cm_toy["b", "a"] * 1 / 10)
  # day 1 uses day 0: nothing infectious, floored
  expect_equal(get(1, "a"), 1e-10)
})

test_that("theta down-weights vaccinated infectious", {
  foi <- force_of_infection(prev_toy(), n_toy, cm_toy, theta = 1)
  expect_equal(foi$Ibar[foi$day == 3], c(1e-10, 1e-10))
  half <- force_of_infection(prev_toy(), n_toy, cm_toy, theta = 0.5)
  expect_equal(half$Ibar[half$day == 3 & half$age_group == "a"],
               cm_toy["b", "a"] * 0.5 / 10)
})

test_that("non-age-structured FoI is total prevalence over population", {
  foi <- force_of_infection(prev_toy(), n_toy, age_structured = FALSE)
  expect_equal(foi$Ibar[foi$day == 2], rep(2 / 30, 2))
})

test_that("force_of_infection checks age-group names", {
  expect_error(force_of_infection(prev_toy(), c(x = 1, y = 2), cm_toy),
               "named by the age groups")
})

test_that("make_stan_data aggregates person-days without losing any", {
  pop <- ve_population(ve_sim_waning)
  coh <- sample_cohort(pop, 300, seed = 2)
  pd  <- expand_person_days(counting_process(coh))
  foi <- force_of_infection(daily_prevalence(daily_incidence(pop)),
                            age_group_sizes(pop))
  sd  <- make_stan_data(pd, foi)
  expect_equal(sum(sd$data$Y), sum(pd$event))
  expect_equal(sum(exp(sd$data$log_persdays)), nrow(pd))
  expect_equal(sd$data$N, length(sd$data$Y))
  expect_true(all(sd$data$tvax_bin[sd$data$vax == 0] == 1L))
  expect_true(all(sd$data$tvax_bin[sd$data$vax == 1] >= 2L))
  expect_lte(max(sd$data$tvax_bin), sd$data$T_vax)

  capped <- make_stan_data(pd, foi, max_days_since_vax = 30)
  expect_equal(capped$data$T_vax, 31L)
  expect_equal(max(capped$data$tvax_bin), 31L)
})

test_that("true_ve ramps, peaks and wanes to the floor", {
  p <- ve_sim_waning$params
  expect_equal(true_ve(0, p), 0)
  expect_equal(true_ve(p$ve_ramp_days, p), p$vaccine_effectiveness)
  expect_equal(true_ve(p$ve_ramp_days + p$ve_half_life, p),
               p$ve_floor + (p$vaccine_effectiveness - p$ve_floor) / 2)
  expect_equal(true_ve(1e4, p), p$ve_floor)
  flat <- ve_sim_no_waning$params
  expect_equal(true_ve(c(30, 100), flat), c(0.8, 0.8))
})

test_that("simulate_counts replaces Y reproducibly", {
  d <- list(N = 3, T_vax = 2, log_Ibar = log(c(1, 1, 1)),
            log_persdays = log(c(100, 100, 100)), vax = c(0, 1, 1),
            tvax_bin = c(1, 2, 2), Y = c(0, 0, 0))
  a <- simulate_counts(d, log_beta = 0, gamma = c(0, -1), seed = 1)
  expect_identical(a, simulate_counts(d, 0, c(0, -1), seed = 1))
  expect_type(a$Y, "integer")
  expect_error(simulate_counts(d, 0, gamma = 0), "length T_vax")
})

test_that("bundled Stan models exist", {
  expect_true(file.exists(ve_stan_file("noncentered")))
  expect_true(file.exists(ve_stan_file("centered")))
})

test_that("plot_ve builds with and without colour and truth", {
  skip_if_not_installed("ggplot2")
  ve <- data.frame(days_since_vax = rep(1:5, 2), median = 0.5, lower = 0.3,
                   upper = 0.7, model = rep(c("a", "b"), each = 5))
  expect_s3_class(plot_ve(ve[1:5, ]), "ggplot")
  p <- plot_ve(ve, truth = ve_sim_waning$params, colour = "model", max_days = 4)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})
