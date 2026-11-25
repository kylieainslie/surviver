# Toy population: 4 people, 10-day study
#   A: never vaccinated, never reported
#   B: vaccinated day 3, reported day 6
#   C: reported day 4, vaccinated day 7 (after infection)
#   D: vaccinated day 2, reported day 12 (after study end -> censored)
toy <- list(
  demographic = data.frame(ind_id = c("A", "B", "C", "D"),
                           age = c(3, 30, 70, 10),
                           sex = c("F", "M", "F", "M")),
  vaccination = data.frame(ind_id = c("B", "C", "D"),
                           time_of_vaccination = c(3, 7, 2)),
  incidence   = data.frame(ind_id = c("B", "C", "D"),
                           time_of_reporting = c(6, 4, 12)),
  study_end   = 10
)

test_that("ve_population joins one row per person and keeps study_end", {
  pop <- ve_population(toy)
  expect_equal(nrow(pop), 4)
  expect_equal(attr(pop, "study_end"), 10)
  expect_equal(pop$vaccinated, c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(pop$infected,   c(FALSE, TRUE, TRUE, TRUE))
})

test_that("ve_population reads the HERMES directory format", {
  dir <- withr::local_tempdir()
  saveRDS(cbind(toy$demographic, time = 0), file.path(dir, "demographic_data.rds"))
  saveRDS(cbind(toy$vaccination, time = 10), file.path(dir, "vaccination_data.rds"))
  saveRDS(cbind(toy$incidence, time = 10), file.path(dir, "incidence_data_reduced.rds"))
  expect_message(pop <- ve_population(dir), "last reporting day")
  expect_equal(attr(pop, "study_end"), 12)
  expect_false("time" %in% names(pop))
})

test_that("ve_population rejects duplicated ids", {
  bad <- toy
  bad$vaccination <- rbind(bad$vaccination, bad$vaccination[1, ])
  expect_error(ve_population(bad), "Duplicated")
})

test_that("daily_incidence counts cases and at-risk correctly", {
  inc <- daily_incidence(ve_population(toy))
  expect_equal(max(inc$day), 10)
  # D's report on day 12 is after study end and not counted
  expect_equal(sum(inc$n_cases), 2)
  # B reported day 6, vaccinated day 3 -> vaccinated case, age 18-49
  b <- inc[inc$day == 6 & inc$age_group == "18-49" & inc$vax_status == 1, ]
  expect_equal(b$n_cases, 1)
  # C reported day 4 before vaccination -> unvaccinated case
  c4 <- inc[inc$day == 4 & inc$age_group == "65+", ]
  expect_equal(c4$n_cases[c4$vax_status == 0], 1)
  # everyone at risk on day 1; B no longer at risk on day 7
  expect_equal(sum(inc$n_at_risk[inc$day == 1]), 4)
  expect_equal(sum(inc$n_at_risk[inc$day == 7 & inc$age_group == "18-49"]), 0)
  # vaccinated status on day t requires vaccination before t
  d <- inc[inc$age_group == "5-17", ]
  expect_equal(d$n_at_risk[d$day == 2 & d$vax_status == 1], 0)
  expect_equal(d$n_at_risk[d$day == 3 & d$vax_status == 1], 1)
})

test_that("daily_prevalence spreads each case over the infectious period", {
  inc  <- daily_incidence(ve_population(toy))
  prev <- daily_prevalence(inc, duration_infection = 7, reporting_delay = 2)
  b <- prev[prev$age_group == "18-49" & prev$vax_status == 1, ]
  # reported day 6 -> infectious days 4..10
  expect_equal(b$day[b$n_infectious == 1], 4:10)
  expect_equal(sum(prev$n_infectious), 14)
})

test_that("counting_process splits at vaccination and censors at study end", {
  cp <- counting_process(ve_population(toy))
  b <- cp[cp$ind_id == "B", ]
  expect_equal(b$tstart, c(0, 3))
  expect_equal(b$tstop, c(3, 6))
  expect_equal(b$vax_status, c(0L, 1L))
  expect_equal(b$event, c(0, 1))
  # C reported before vaccination: one unvaccinated interval ending in event
  c_ <- cp[cp$ind_id == "C", ]
  expect_equal(nrow(c_), 1)
  expect_equal(c_$vax_status, 0L)
  # D reported after study end: censored at 10
  d <- cp[cp$ind_id == "D", ]
  expect_equal(max(d$tstop), 10)
  expect_equal(sum(d$event), 0)
  # A never vaccinated: vax_status 0, not NA
  expect_equal(cp$vax_status[cp$ind_id == "A"], 0L)
})

test_that("expand_person_days gives one row per day with a single event", {
  pop <- ve_population(toy)
  pd  <- expand_person_days(counting_process(pop), daily_incidence(pop))
  expect_equal(nrow(pd), 10 + 6 + 4 + 10)
  expect_equal(sum(pd$event), 2)
  b <- pd[pd$ind_id == "B", ]
  expect_equal(b$day, 1:6)
  expect_equal(b$event, c(0, 0, 0, 0, 0, 1))
  expect_equal(b$vax_status, c(0, 0, 0, 1, 1, 1))
  expect_equal(b$days_since_vax, c(NA, NA, NA, 1, 2, 3))
  expect_false(anyNA(pd$incidence_rate))
})

test_that("pipeline is consistent on simulated data", {
  pop <- ve_population(ve_sim_waning)
  se  <- attr(pop, "study_end")
  coh <- sample_cohort(pop, 500, seed = 1)
  expect_equal(nrow(coh), 500)
  expect_identical(coh, sample_cohort(pop, 500, seed = 1))
  pd <- expand_person_days(counting_process(coh), daily_incidence(pop))
  expect_equal(sum(pd$event), sum(coh$infected & coh$time_of_reporting <= se))
  followup <- ifelse(coh$infected & coh$time_of_reporting <= se,
                     coh$time_of_reporting, se)
  expect_equal(nrow(pd), sum(followup))
})

test_that("missing study_end gives an informative error", {
  pop <- ve_population(toy)
  attr(pop, "study_end") <- NULL
  expect_error(daily_incidence(pop), "study_end")
})
