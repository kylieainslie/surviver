# surviver

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**surviver** estimates vaccine effectiveness as a function of time since
vaccination, VE(t), from individual-level data on vaccination and reported
infections.

It combines survival analysis with infectious disease modelling. Each
person's hazard of infection is driven by the force of infection they are
exposed to, built from daily prevalence and an age-structured contact
matrix, and vaccination modifies that hazard by a log hazard ratio γ(t)
that varies smoothly with time since vaccination:

$$
\log h = \log\beta + \log \bar I(t_c, a) + \gamma(t)\,V,
\qquad \mathrm{VE}(t) = 1 - e^{\gamma(t)}
$$

The model is a Poisson approximation to a Cox model, fitted in Stan with a
random-walk prior on γ(t), so VE can rise, plateau and wane without assuming
a functional form.

The method was developed by Group D at the Lorentz Center workshop
*Connecting Survival Analysis and Infectious Disease Modelling* (Leiden, June
2026).

> **Status: experimental.** The method is under active development and
> validation. See the [lab notebook](https://kylieainslie.github.io/surviver/articles/lab-notebook.html) for
> what has been checked so far and the known open issues.

## Installation

```r
# install.packages("remotes")
remotes::install_github("kylieainslie/surviver")
```

Fitting models needs [cmdstanr](https://mc-stan.org/cmdstanr/) and CmdStan,
which are not on CRAN:

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

## Example

Estimate VE(t) from a simulated dataset where the true VE is known:

```r
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

![Estimated VE (blue, with 90% credible interval) against the true VE (red dashed) for the two simulated scenarios.](man/figures/README-ve-fit.png)

See the [Get started](https://kylieainslie.github.io/surviver/articles/surviver.html)
vignette for the epidemic dynamics, the model fit to cases, and this
comparison in detail.

To use your own data, pass `ve_population()` a folder containing
`demographic_data.rds`, `vaccination_data.rds` and
`incidence_data_reduced.rds`, or a list with `demographic`, `vaccination`
and `incidence` data frames (see `?ve_population`).

## Main functions

| Step | Functions |
|---|---|
| Prepare data | `ve_population()`, `sample_cohort()`, `add_age_group()` |
| Exposure | `daily_incidence()`, `daily_prevalence()`, `force_of_infection()`, `age_group_sizes()` |
| Survival data | `counting_process()`, `reported_soon_after_vaccination()`, `expand_person_days()`, `make_stan_data()` |
| Fit and summarise | `fit_ve()`, `extract_ve()`, `plot_ve()` |
| Validation | `true_ve()`, `simulate_counts()`, `contact_matrix_hermes()` |

Two Stan models are bundled (`ve_stan_file()`): `"noncentered"` (default)
and `"centered"` parameterisations of the random walk.

## Data

| Dataset | Description |
|---|---|
| `ve_sim_waning` | Simulated outbreak, 10,000 people, 180 days. VE rises to 0.8 over 14 days, then wanes (half-life 45 days) towards 0.1 |
| `ve_sim_no_waning` | As above, but VE stays at 0.8 after the ramp |
| `polymod_nl` | POLYMOD contact matrix for the Netherlands, five age groups |

The simulated datasets were generated with **HERMES** (Health, Epidemic and
Economic R-based Microsimulation Engine for Individual-based Simulations),
an individual-based model of transmission through households, schools,
workplaces and the community, developed by Lander Willem (University of
Antwerp). See the [HERMES repository](https://github.com/lwillem/ibm_hermes) and its
[user manual](https://github.com/lwillem/ibm_hermes/blob/lorentz_center_workshop/doc/user_manual_HERMES.pdf). For this package HERMES was extended with vaccine
protection that wanes with time since vaccination
([`waning-vaccine-effectiveness` branch](https://github.com/kylieainslie/ibm_hermes/tree/waning-vaccine-effectiveness) of our fork). The script is
`data-raw/ve_test_scenarios.R`; it needs a local clone of HERMES (set
`HERMES_PATH`), and each dataset records the HERMES commit it came from.

## Development

- [Handover notes](https://kylieainslie.github.io/surviver/articles/handover.html):
  where everything lives, setup, how to regenerate or adapt the simulations,
  and what to do next.
- `vignettes/articles/lab-notebook.Rmd`: dated record of how the method and
  data were developed and every validation check, with results as observed.
- `dev/checks/`: scripts that reproduce each lab-notebook entry.
- `data-raw/`: scripts that build the package datasets.

## Acknowledgements

The method was developed by Group D at the Lorentz Center workshop
*Connecting Survival Analysis and Infectious Disease Modelling* (Leiden,
June 2026): Fanny Bergström, David Price, Nathalie Dean, Hein Putter and
Kylie Ainslie. The simulated data use HERMES by Lander Willem
(<https://github.com/lwillem/ibm_hermes>).

## Licence

MIT © Kylie Ainslie
