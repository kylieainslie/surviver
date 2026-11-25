# POLYMOD contact matrix for the Netherlands

Mean daily contacts by age group from the POLYMOD survey (Netherlands),
symmetrised, aggregated to surviver's default age groups
([default_age_labels](https://kylieainslie.github.io/surviver/reference/default_age_breaks.md)).

## Usage

``` r
polymod_nl
```

## Format

A 5 x 5 numeric matrix. Rows are participant age groups and columns are
contact age groups, so `polymod_nl[i, j]` is the mean number of daily
contacts a person in group `i` has with people in group `j`.

## Source

Mossong J, et al. (2008) Social contacts and mixing patterns relevant to
the spread of infectious diseases. *PLoS Med* 5(3):e74. Computed with
[`socialmixr::contact_matrix()`](https://epiforecasts.io/socialmixr/reference/contact_matrix.html);
see `data-raw/polymod_nl.R`.
