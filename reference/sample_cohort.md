# Randomly sample a cohort

Randomly sample a cohort

## Usage

``` r
sample_cohort(pop, n, seed = NULL)
```

## Arguments

- pop:

  Person-level data from
  [`ve_population()`](https://kylieainslie.github.io/surviver/reference/ve_population.md).

- n:

  Cohort size. If `n >= nrow(pop)` the whole population is returned.

- seed:

  Optional random seed, applied locally.

## Value

A subset of `pop` (attributes preserved).
