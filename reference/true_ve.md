# True VE curve used in the HERMES simulations

Linear rise from 0 to `vaccine_effectiveness` over `ve_ramp_days`, then
exponential decay with half-life `ve_half_life` towards `ve_floor`.
Matches `compute_time_varying_ve()` in the [HERMES waning-VE
extension](https://github.com/kylieainslie/ibm_hermes/tree/waning-vaccine-effectiveness).

## Usage

``` r
true_ve(t, params)
```

## Arguments

- t:

  Days since vaccination.

- params:

  A list with `vaccine_effectiveness`, `ve_ramp_days`, `ve_half_life`,
  `ve_floor`, e.g. `ve_sim_waning$params`.

## Value

Numeric vector of VE values.

## Examples

``` r
true_ve(c(7, 14, 59), ve_sim_waning$params)
#> [1] 0.40 0.80 0.45
```
