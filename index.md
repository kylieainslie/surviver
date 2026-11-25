# surviver

**surviver** estimates vaccine effectiveness as a function of time since
vaccination, VE(t), from individual-level data on vaccination and
reported infections.

It combines survival analysis with infectious disease modelling. Each
person’s hazard of infection is driven by the force of infection they
are exposed to, built from daily prevalence and an age-structured
contact matrix, and vaccination modifies that hazard by a log hazard
ratio γ(t) that varies smoothly with time since vaccination:

``` math
\log h = \log\beta + \log \bar I(t_c, a) + \gamma(t)\,V,
\qquad \mathrm{VE}(t) = 1 - e^{\gamma(t)}
```

The model is a Poisson approximation to a Cox model, fitted in Stan with
a random-walk prior on γ(t), so VE can rise, plateau and wane without
assuming a functional form.

The method was developed by Group D at the Lorentz Center workshop
*Connecting Survival Analysis and Infectious Disease Modelling* (Leiden,
June 2026).

> **Status: experimental.** The method is under active development and
> validation. See the [lab
> notebook](https://kylieainslie.github.io/surviver/articles/lab-notebook.html)
> for what has been checked so far and the known open issues.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("kylieainslie/surviver")
```

Fitting models needs [cmdstanr](https://mc-stan.org/cmdstanr/) and
CmdStan, which are not on CRAN:

``` r

install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

## Example

Estimate VE(t) from a simulated dataset where the true VE is known:

``` r

library(surviver)

# person-level data: demographics, vaccination and reported cases
pop <- ve_population(ve_sim_waning)

# force of infection by day and age group, from reported cases
inc  <- daily_incidence(pop)
prev <- daily_prevalence(inc, duration_infection = 7, reporting_delay = 2)
foi  <- force_of_infection(prev, age_group_sizes(pop), contact_matrix = polymod_nl)

# person-days split at vaccination, aggregated for Stan. People reported
# within 7 days of vaccination are censored at vaccination, since they were
# probably infected before it
early <- reported_soon_after_vaccination(pop, window = 7)
pd <- expand_person_days(counting_process(pop, censor_at_vaccination = early))
sd <- make_stan_data(pd, foi)

# fit and plot against the true curve (VE is not estimable for the
# first 7 days after vaccination once early cases are censored)
fit <- fit_ve(sd, seed = 42)
ve  <- extract_ve(fit)
plot_ve(ve[ve$days_since_vax > 7, ], truth = ve_sim_waning$params)
```

To use your own data, pass
[`ve_population()`](https://kylieainslie.github.io/surviver/reference/ve_population.md)
a folder containing `demographic_data.rds`, `vaccination_data.rds` and
`incidence_data_reduced.rds`, or a list with `demographic`,
`vaccination` and `incidence` data frames (see
[`?ve_population`](https://kylieainslie.github.io/surviver/reference/ve_population.md)).

## Main functions

| Step | Functions |
|----|----|
| Prepare data | [`ve_population()`](https://kylieainslie.github.io/surviver/reference/ve_population.md), [`sample_cohort()`](https://kylieainslie.github.io/surviver/reference/sample_cohort.md), [`add_age_group()`](https://kylieainslie.github.io/surviver/reference/add_age_group.md) |
| Exposure | [`daily_incidence()`](https://kylieainslie.github.io/surviver/reference/daily_incidence.md), [`daily_prevalence()`](https://kylieainslie.github.io/surviver/reference/daily_prevalence.md), [`force_of_infection()`](https://kylieainslie.github.io/surviver/reference/force_of_infection.md), [`age_group_sizes()`](https://kylieainslie.github.io/surviver/reference/age_group_sizes.md) |
| Survival data | [`counting_process()`](https://kylieainslie.github.io/surviver/reference/counting_process.md), [`reported_soon_after_vaccination()`](https://kylieainslie.github.io/surviver/reference/reported_soon_after_vaccination.md), [`expand_person_days()`](https://kylieainslie.github.io/surviver/reference/expand_person_days.md), [`make_stan_data()`](https://kylieainslie.github.io/surviver/reference/make_stan_data.md) |
| Fit and summarise | [`fit_ve()`](https://kylieainslie.github.io/surviver/reference/fit_ve.md), [`extract_ve()`](https://kylieainslie.github.io/surviver/reference/extract_ve.md), [`plot_ve()`](https://kylieainslie.github.io/surviver/reference/plot_ve.md) |
| Validation | [`true_ve()`](https://kylieainslie.github.io/surviver/reference/true_ve.md), [`simulate_counts()`](https://kylieainslie.github.io/surviver/reference/simulate_counts.md), [`contact_matrix_hermes()`](https://kylieainslie.github.io/surviver/reference/contact_matrix_hermes.md) |

Two Stan models are bundled
([`ve_stan_file()`](https://kylieainslie.github.io/surviver/reference/ve_stan_file.md)):
`"noncentered"` (default) and `"centered"` parameterisations of the
random walk.

## Data

| Dataset | Description |
|----|----|
| `ve_sim_waning` | Simulated outbreak, 10,000 people, 180 days. VE rises to 0.8 over 14 days, then wanes (half-life 45 days) towards 0.1 |
| `ve_sim_no_waning` | As above, but VE stays at 0.8 after the ramp |
| `polymod_nl` | POLYMOD contact matrix for the Netherlands, five age groups |

The simulated datasets were generated with the
[HERMES](https://github.com/kylieainslie/ibm_hermes/tree/waning-vaccine-effectiveness)
individual-based model (household, school, workplace and community
transmission), extended with waning vaccine protection. The script is
`data-raw/ve_test_scenarios.R`; it needs a local clone of HERMES (set
`HERMES_PATH`), and each dataset records the HERMES commit it came from.

## Development

- [Handover
  notes](https://kylieainslie.github.io/surviver/articles/handover.html):
  where everything lives, setup, how to regenerate or adapt the
  simulations, and what to do next.
- `vignettes/articles/lab-notebook.Rmd`: dated record of how the method
  and data were developed and every validation check, with results as
  observed.
- `dev/checks/`: scripts that reproduce each lab-notebook entry.
- `data-raw/`: scripts that build the package datasets.

## Licence

MIT © Kylie Ainslie
