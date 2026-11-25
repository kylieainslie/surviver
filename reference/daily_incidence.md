# Daily cases and population at risk

Counts reported cases and the population at risk on each study day by
age group and vaccination status.

## Usage

``` r
daily_incidence(
  pop,
  study_end = attr(pop, "study_end"),
  age_breaks = default_age_breaks,
  age_labels = default_age_labels
)
```

## Arguments

- pop:

  Person-level data from
  [`ve_population()`](https://kylieainslie.github.io/surviver/reference/ve_population.md).

- study_end:

  Last day of follow-up.

- age_breaks, age_labels:

  Passed to [`cut()`](https://rdrr.io/r/base/cut.html) (left-closed
  intervals).

## Value

A tibble with one row per `day` x `age_group` x `vax_status`:
`n_at_risk`, `n_cases` and `incidence_rate` (`n_cases / n_at_risk`).

## Details

A person is at risk on day `t` if they have not been reported before day
`t`. They count as vaccinated on day `t` if vaccinated before day `t`,
and a case counts as vaccinated if vaccination preceded reporting.
