# Counting-process data with time-varying vaccination

Splits each person's follow-up at vaccination, giving `(tstart, tstop]`
intervals with `vax_status` 0 before and 1 after vaccination. Follow-up
ends at reporting (`event = 1`) or `study_end` (censored); cases
reported after `study_end` are censored at `study_end`.

## Usage

``` r
counting_process(
  cohort,
  study_end = attr(cohort, "study_end"),
  censor_at_vaccination = NULL,
  age_breaks = default_age_breaks,
  age_labels = default_age_labels
)
```

## Arguments

- cohort:

  Person-level data from
  [`ve_population()`](https://kylieainslie.github.io/surviver/reference/ve_population.md)
  or
  [`sample_cohort()`](https://kylieainslie.github.io/surviver/reference/sample_cohort.md).

- study_end:

  Last day of follow-up.

- censor_at_vaccination:

  Optional vector of `ind_id`s to censor at their vaccination date.
  Ignored for people who are unvaccinated or whose follow-up ends before
  vaccination.

- age_breaks, age_labels:

  Passed to [`cut()`](https://rdrr.io/r/base/cut.html) (left-closed
  intervals).

## Value

A data frame from
[`survival::tmerge()`](https://rdrr.io/pkg/survival/man/tmerge.html)
with `ind_id`, `age`, `age_group`, `sex`, `time_of_vaccination`,
`tstart`, `tstop`, `event`, `vax_status`.

## Details

People reported before their vaccination date are already handled: their
follow-up ends at the report, so they never count as vaccinated. People
infected before vaccination but reported after it would otherwise count
as vaccinated cases. `censor_at_vaccination` lets these be censored at
their vaccination date instead (follow-up ends unvaccinated, no event).
Which people to censor is the analyst's choice; see
[`reported_soon_after_vaccination()`](https://kylieainslie.github.io/surviver/reference/reported_soon_after_vaccination.md)
for a rule based only on observed dates.
