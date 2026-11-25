# Summarise VE by day since vaccination

Summarise VE by day since vaccination

## Usage

``` r
extract_ve(fit, probs = c(0.05, 0.95))
```

## Arguments

- fit:

  A fit from
  [`fit_ve()`](https://kylieainslie.github.io/surviver/reference/fit_ve.md).

- probs:

  Lower and upper quantiles for the credible interval.

## Value

A tibble with `days_since_vax`, `median`, `lower`, `upper` (VE scale),
plus `rhat` and `ess_bulk`.
