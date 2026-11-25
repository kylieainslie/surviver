# Package index

## Prepare data

Build person-level, daily and counting-process data from vaccination and
case records.

- [`ve_population()`](https://kylieainslie.github.io/surviver/reference/ve_population.md)
  : Build a person-level table from VE input data
- [`sample_cohort()`](https://kylieainslie.github.io/surviver/reference/sample_cohort.md)
  : Randomly sample a cohort
- [`add_age_group()`](https://kylieainslie.github.io/surviver/reference/add_age_group.md)
  : Add an age-group column
- [`default_age_breaks`](https://kylieainslie.github.io/surviver/reference/default_age_breaks.md)
  [`default_age_labels`](https://kylieainslie.github.io/surviver/reference/default_age_breaks.md)
  : Default age groups
- [`daily_incidence()`](https://kylieainslie.github.io/surviver/reference/daily_incidence.md)
  : Daily cases and population at risk
- [`counting_process()`](https://kylieainslie.github.io/surviver/reference/counting_process.md)
  : Counting-process data with time-varying vaccination
- [`reported_soon_after_vaccination()`](https://kylieainslie.github.io/surviver/reference/reported_soon_after_vaccination.md)
  : People reported soon after vaccination
- [`expand_person_days()`](https://kylieainslie.github.io/surviver/reference/expand_person_days.md)
  : Expand counting-process data to person-days

## Exposure

Daily infectious pressure (force of infection) by age group.

- [`daily_prevalence()`](https://kylieainslie.github.io/surviver/reference/daily_prevalence.md)
  : Daily number infectious
- [`age_group_sizes()`](https://kylieainslie.github.io/surviver/reference/age_group_sizes.md)
  : Population size by age group
- [`force_of_infection()`](https://kylieainslie.github.io/surviver/reference/force_of_infection.md)
  : Daily force-of-infection offset
- [`contact_matrix_hermes()`](https://kylieainslie.github.io/surviver/reference/contact_matrix_hermes.md)
  : Contact matrix implied by a HERMES population

## Fit the VE model

Poisson approximation to a Cox model with a random walk on the log
hazard ratio over days since vaccination.

- [`make_stan_data()`](https://kylieainslie.github.io/surviver/reference/make_stan_data.md)
  : Build Stan input for the time-varying VE model
- [`ve_stan_file()`](https://kylieainslie.github.io/surviver/reference/ve_stan_file.md)
  : Path to a bundled Stan model
- [`fit_ve()`](https://kylieainslie.github.io/surviver/reference/fit_ve.md)
  : Fit the time-varying VE model
- [`extract_ve()`](https://kylieainslie.github.io/surviver/reference/extract_ve.md)
  : Summarise VE by day since vaccination
- [`plot_ve()`](https://kylieainslie.github.io/surviver/reference/plot_ve.md)
  : Plot estimated VE by day since vaccination
- [`expected_cases()`](https://kylieainslie.github.io/surviver/reference/expected_cases.md)
  : Observed and model-expected case counts

## Validation

Helpers for checking the method against known truth.

- [`true_ve()`](https://kylieainslie.github.io/surviver/reference/true_ve.md)
  : True VE curve used in the HERMES simulations
- [`simulate_counts()`](https://kylieainslie.github.io/surviver/reference/simulate_counts.md)
  : Simulate event counts from the VE model

## Data

- [`ve_sim_waning`](https://kylieainslie.github.io/surviver/reference/ve_sim.md)
  [`ve_sim_no_waning`](https://kylieainslie.github.io/surviver/reference/ve_sim.md)
  : Simulated ground-truth data for time-varying vaccine effectiveness
- [`polymod_nl`](https://kylieainslie.github.io/surviver/reference/polymod_nl.md)
  : POLYMOD contact matrix for the Netherlands
