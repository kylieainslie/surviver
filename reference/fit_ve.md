# Fit the time-varying VE model

Fits a Poisson approximation to a Cox model in which the log hazard
ratio for vaccination, \\\gamma(t)\\, follows a random walk over days
since vaccination: \$\$\log \mu = \log\beta + \log \bar I +
\log(\text{person-days}) + \gamma(t)\\\text{vax}.\$\$ VE is \\1 -
\exp(\gamma(t))\\. Requires cmdstanr and CmdStan.

## Usage

``` r
fit_ve(
  stan_data,
  model = c("noncentered", "centered"),
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  max_treedepth = 10,
  seed = NULL,
  ...
)
```

## Arguments

- stan_data:

  Output of
  [`make_stan_data()`](https://kylieainslie.github.io/surviver/reference/make_stan_data.md),
  or its `data` element.

- model:

  Stan model; see
  [`ve_stan_file()`](https://kylieainslie.github.io/surviver/reference/ve_stan_file.md).

- chains, iter_warmup, iter_sampling, adapt_delta, max_treedepth, seed:

  Passed to `$sample()`.

- ...:

  Further arguments to `$sample()`.

## Value

A `CmdStanMCMC` fit.
