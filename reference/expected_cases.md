# Observed and model-expected case counts

Sums observed cases and the posterior distribution of expected cases
(\\\mu\\ from the fitted model) over groups of Stan cells, e.g. by
calendar day and vaccination status. Comparing the two shows how well
the model fits the case data.

## Usage

``` r
expected_cases(
  fit,
  stan_data,
  by = c("day", "vax_status"),
  probs = c(0.05, 0.95)
)
```

## Arguments

- fit:

  A fit from
  [`fit_ve()`](https://kylieainslie.github.io/surviver/reference/fit_ve.md).

- stan_data:

  The output of
  [`make_stan_data()`](https://kylieainslie.github.io/surviver/reference/make_stan_data.md)
  used for the fit (with its `cells`).

- by:

  Columns of `stan_data$cells` to group by.

- probs:

  Lower and upper quantiles for the interval.

## Value

A tibble with the `by` columns, `observed`, and the posterior `mean`,
`lower` and `upper` of the expected count.
