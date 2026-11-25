# Plot estimated VE by day since vaccination

Plot estimated VE by day since vaccination

## Usage

``` r
plot_ve(ve, truth = NULL, colour = NULL, max_days = NULL)
```

## Arguments

- ve:

  Output of
  [`extract_ve()`](https://kylieainslie.github.io/surviver/reference/extract_ve.md).
  Several fits can be compared by binding their outputs with a `model`
  (or other) column named in `colour`.

- truth:

  Optional HERMES parameter list (e.g. `ve_sim_waning$params`); the true
  VE(t) is drawn as a dashed line.

- colour:

  Optional column of `ve` to colour by.

- max_days:

  Upper limit of the x axis.

## Value

A ggplot object.
