# Build a person-level table from VE input data

Joins demographic, vaccination and reported-case data into one row per
person.

## Usage

``` r
ve_population(x, study_end = NULL)
```

## Arguments

- x:

  Either a directory containing `demographic_data.rds`,
  `vaccination_data.rds` and `incidence_data_reduced.rds` (the HERMES
  export format, optionally with `scenario_params.rds`), or a list with
  elements `demographic`, `vaccination` and `incidence` such as
  [ve_sim_waning](https://kylieainslie.github.io/surviver/reference/ve_sim.md).

- study_end:

  Last day of follow-up. Defaults to `x$study_end` (or the value in
  `scenario_params.rds`); if neither exists, the last whole reporting
  day, with a message.

## Value

A data frame with one row per person: the demographic columns plus
`time_of_vaccination`, `time_of_reporting` (`NA` if never vaccinated /
reported), and logical `vaccinated` and `infected`. The study end is
stored in the `"study_end"` attribute.

## Examples

``` r
pop <- ve_population(ve_sim_waning)
attr(pop, "study_end")
#> [1] 180
```
