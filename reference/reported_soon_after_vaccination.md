# People reported soon after vaccination

A rule, based only on observed dates, for identifying vaccinated people
who may have been infected before vaccination: those reported within
`window` days after their vaccination date. Infection precedes reporting
by the incubation period plus the reporting delay, so a case reported
shortly after vaccination may have been infected before it. Pass the
result to `counting_process(censor_at_vaccination = )`.

## Usage

``` r
reported_soon_after_vaccination(pop, window)
```

## Arguments

- pop:

  Person-level data from
  [`ve_population()`](https://kylieainslie.github.io/surviver/reference/ve_population.md).

- window:

  Days after vaccination within which a report triggers censoring
  (reports on days `1` to `window` after vaccination).

## Value

Character vector of `ind_id`s.

## Details

Censoring also removes some genuine early vaccine failures, so VE is not
estimable for the first `window` days after vaccination. Longer windows
catch more prior infections but discard more genuine cases.

## Examples

``` r
pop <- ve_population(ve_sim_waning)
ids <- reported_soon_after_vaccination(pop, window = 7)
length(ids)
#> [1] 25
cp <- counting_process(pop, censor_at_vaccination = ids)
```
