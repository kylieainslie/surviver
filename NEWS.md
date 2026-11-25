# surviver (development version)

* Initial package, consolidating the Group D analysis from the Lorentz
  Center workshop (June 2026).
* Data preparation: `ve_population()`, `sample_cohort()`,
  `daily_incidence()`, `daily_prevalence()`, `counting_process()`,
  `expand_person_days()`.
* Exposure: `force_of_infection()` (age-structured or overall),
  `age_group_sizes()`, `contact_matrix_hermes()`.
* Model: `make_stan_data()`, `fit_ve()` with non-centred (default) and
  centred Stan models, `extract_ve()`, `plot_ve()`.
* Validation helpers: `true_ve()`, `simulate_counts()`.
* Datasets: `ve_sim_waning`, `ve_sim_no_waning` (HERMES simulations with
  known VE(t)) and `polymod_nl`.
* Lab notebook article recording development and validation checks.
