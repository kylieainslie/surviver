# Build Stan input for the time-varying VE model

Aggregates person-days into cells of age group x vaccination status x
day x day since vaccination, and attaches the force-of-infection offset.

## Usage

``` r
make_stan_data(person_days, foi, max_days_since_vax = NULL)
```

## Arguments

- person_days:

  Output of
  [`expand_person_days()`](https://kylieainslie.github.io/surviver/reference/expand_person_days.md).

- foi:

  Output of
  [`force_of_infection()`](https://kylieainslie.github.io/surviver/reference/force_of_infection.md).

- max_days_since_vax:

  Last day since vaccination with its own bin. Defaults to the longest
  observed.

## Value

A list with `data` (the list passed to Stan) and `cells` (the aggregated
data frame, for inspection).

## Details

Days since vaccination are binned by day: bin 1 is unvaccinated, bin
`t + 1` is day `t` after vaccination. Days beyond `max_days_since_vax`
share the last bin.
