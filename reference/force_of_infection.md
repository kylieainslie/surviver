# Daily force-of-infection offset

Computes the infectious pressure \\\bar I(t, a)\\ on each study day for
each age group, used as an offset in the VE model.

## Usage

``` r
force_of_infection(
  prevalence,
  pop_size,
  contact_matrix = surviver::polymod_nl,
  theta = 0,
  age_structured = TRUE,
  lag = 1L,
  study_end = max(prevalence$day),
  floor = 1e-10
)
```

## Arguments

- prevalence:

  Output of
  [`daily_prevalence()`](https://kylieainslie.github.io/surviver/reference/daily_prevalence.md).

- pop_size:

  Named vector of population size by age group, e.g. from
  [`age_group_sizes()`](https://kylieainslie.github.io/surviver/reference/age_group_sizes.md).

- contact_matrix:

  Square matrix of contacts with rows and columns in the order of
  `names(pop_size)`; `contact_matrix[a', a]` is the mean daily contacts
  a person in `a'` has with group `a`. Ignored when
  `age_structured = FALSE`.

- theta:

  Relative reduction in infectiousness of vaccinated cases (0 = as
  infectious as unvaccinated, 1 = not infectious). The HERMES
  simulations assume vaccination does not affect infectiousness, i.e.
  `theta = 0`.

- age_structured:

  Use the contact matrix (`TRUE`) or overall prevalence (`FALSE`).

- lag:

  Days between prevalence and the infections it causes.

- study_end:

  Last day to compute; defaults to the last prevalence day.

- floor:

  Lower bound applied so the log offset is finite.

## Value

A tibble with `day` (1 to `study_end`), `age_group` and `Ibar`.

## Details

With age structure, for recipient age group \\a\\ on day \\t\\: \$\$\bar
I(t, a) = \sum\_{a'} \frac{C\_{a'a}}{N_a} I^{\mathrm{eff}}\_{a'}(t -
\ell)\$\$ where \\I^{\mathrm{eff}}\_{a'} = I^u\_{a'} + (1 - \theta)
I^v\_{a'}\\ is the number infectious in group \\a'\\ (unvaccinated plus
vaccinated, the latter down-weighted by \\\theta\\) and \\\ell\\ is
`lag`.

Without age structure every group gets \\\bar I(t) = \sum\_{a'}
I^{\mathrm{eff}}\_{a'}(t - \ell) / N\\.
