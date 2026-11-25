# Add an age-group column

Add an age-group column

## Usage

``` r
add_age_group(
  data,
  age_breaks = default_age_breaks,
  age_labels = default_age_labels
)
```

## Arguments

- data:

  Data frame with an `age` column.

- age_breaks, age_labels:

  Passed to [`cut()`](https://rdrr.io/r/base/cut.html) (left-closed
  intervals).

## Value

`data` with a factor column `age_group`.
