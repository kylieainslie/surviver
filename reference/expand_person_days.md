# Expand counting-process data to person-days

Turns each `(tstart, tstop]` interval into one row per day, with
`event = 1` only on the day of reporting. Optionally attaches the daily
incidence rate for the person's age group and vaccination status.

## Usage

``` r
expand_person_days(cp, incidence = NULL)
```

## Arguments

- cp:

  Output of
  [`counting_process()`](https://kylieainslie.github.io/surviver/reference/counting_process.md).

- incidence:

  Optional output of
  [`daily_incidence()`](https://kylieainslie.github.io/surviver/reference/daily_incidence.md).

## Value

A tibble with one row per person-day: `ind_id`, `age`, `age_group`,
`sex`, `day`, `tstart`, `tstop` (= `day`), `event`, `vax_status`,
`time_of_vaccination`, `days_since_vax` (`NA` when unvaccinated) and, if
`incidence` is given, `incidence_rate`.
