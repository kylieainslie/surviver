# Changelog

## surviver 0.0.0.9000

- Initial package, consolidating the Group D analysis from the Lorentz
  Center workshop (June 2026).
- Data preparation:
  [`ve_population()`](https://kylieainslie.github.io/surviver/reference/ve_population.md),
  [`sample_cohort()`](https://kylieainslie.github.io/surviver/reference/sample_cohort.md),
  [`daily_incidence()`](https://kylieainslie.github.io/surviver/reference/daily_incidence.md),
  [`daily_prevalence()`](https://kylieainslie.github.io/surviver/reference/daily_prevalence.md),
  [`counting_process()`](https://kylieainslie.github.io/surviver/reference/counting_process.md),
  [`expand_person_days()`](https://kylieainslie.github.io/surviver/reference/expand_person_days.md).
- Exposure:
  [`force_of_infection()`](https://kylieainslie.github.io/surviver/reference/force_of_infection.md)
  (age-structured or overall),
  [`age_group_sizes()`](https://kylieainslie.github.io/surviver/reference/age_group_sizes.md),
  [`contact_matrix_hermes()`](https://kylieainslie.github.io/surviver/reference/contact_matrix_hermes.md).
- Model:
  [`make_stan_data()`](https://kylieainslie.github.io/surviver/reference/make_stan_data.md),
  [`fit_ve()`](https://kylieainslie.github.io/surviver/reference/fit_ve.md)
  with non-centred (default) and centred Stan models,
  [`extract_ve()`](https://kylieainslie.github.io/surviver/reference/extract_ve.md),
  [`plot_ve()`](https://kylieainslie.github.io/surviver/reference/plot_ve.md).
- `counting_process(censor_at_vaccination = )` and
  [`reported_soon_after_vaccination()`](https://kylieainslie.github.io/surviver/reference/reported_soon_after_vaccination.md)
  to censor people infected before vaccination but reported after it.
- Validation helpers:
  [`true_ve()`](https://kylieainslie.github.io/surviver/reference/true_ve.md),
  [`simulate_counts()`](https://kylieainslie.github.io/surviver/reference/simulate_counts.md).
- Datasets: `ve_sim_waning`, `ve_sim_no_waning` (HERMES simulations with
  known VE(t)) and `polymod_nl`.
- Lab notebook article recording development and validation checks.
