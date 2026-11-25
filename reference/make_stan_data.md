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

A list with two elements:

**`cells`**: a tibble with one row per cell, i.e. per combination of age
group, vaccination status, study day and day since vaccination that has
any person-time:

- `age_group`:

  Age group (factor).

- `vax_status`:

  1 if the person-days are after vaccination, 0 if before vaccination or
  never vaccinated.

- `day`:

  Study (calendar) day.

- `tvax_bin`:

  Day-since-vaccination bin: 1 for unvaccinated person-days, `t + 1` for
  day `t` after vaccination.

- `Y`:

  Number of reported cases (events) in the cell.

- `n_persdays`:

  Person-days at risk in the cell.

- `Ibar`:

  Force-of-infection offset for that day and age group, from `foi`.

**`data`**: the list passed to Stan, built from `cells` (same row
order):

- `N`:

  Number of cells.

- `K`:

  Number of age groups (declared in the Stan models but not used).

- `T_vax`:

  Number of day-since-vaccination bins, `max_days_since_vax + 1`.

- `Y`:

  `cells$Y`.

- `log_Ibar`, `log_persdays`:

  Logs of `cells$Ibar` and `cells$n_persdays`, used as offsets.

- `vax`, `tvax_bin`:

  `cells$vax_status` and `cells$tvax_bin`.

## Details

Days since vaccination are binned by day: bin 1 is unvaccinated, bin
`t + 1` is day `t` after vaccination. Days beyond `max_days_since_vax`
share the last bin.
