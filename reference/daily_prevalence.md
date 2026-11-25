# Daily number infectious

Approximates the number of people infectious on each day from reported
cases, assuming each case was infectious for `duration_infection` days
starting `reporting_delay` days before they were reported.

## Usage

``` r
daily_prevalence(incidence, duration_infection = 7, reporting_delay = 2)
```

## Arguments

- incidence:

  Output of
  [`daily_incidence()`](https://kylieainslie.github.io/surviver/reference/daily_incidence.md).

- duration_infection:

  Days infectious.

- reporting_delay:

  Days from start of infectiousness to reporting.

## Value

A tibble with one row per `day` x `age_group` x `vax_status`
(vaccination status at reporting) and column `n_infectious`, from day 0
(or earlier, if a case was infectious before day 0) to the last study
day.
