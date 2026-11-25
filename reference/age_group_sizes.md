# Population size by age group

Population size by age group

## Usage

``` r
age_group_sizes(
  pop,
  age_breaks = default_age_breaks,
  age_labels = default_age_labels
)
```

## Arguments

- pop:

  Person-level data from
  [`ve_population()`](https://kylieainslie.github.io/surviver/reference/ve_population.md).

- age_breaks, age_labels:

  Passed to [`cut()`](https://rdrr.io/r/base/cut.html) (left-closed
  intervals).

## Value

A named integer vector of group sizes.
