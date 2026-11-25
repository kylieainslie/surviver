# Simulate event counts from the VE model

Replaces the outcome in Stan input with Poisson counts generated under
known parameters, keeping the real covariate structure. Used to check
that
[`fit_ve()`](https://kylieainslie.github.io/surviver/reference/fit_ve.md)
recovers known values.

## Usage

``` r
simulate_counts(stan_data, log_beta, gamma, seed = NULL)
```

## Arguments

- stan_data:

  Output of
  [`make_stan_data()`](https://kylieainslie.github.io/surviver/reference/make_stan_data.md),
  or its `data` element.

- log_beta:

  True log transmission scaling.

- gamma:

  True log hazard ratio per bin (length `T_vax`; element 1 is the
  unvaccinated bin and is ignored).

- seed:

  Optional random seed, applied locally.

## Value

`stan_data` (same form as supplied) with `Y` replaced.
