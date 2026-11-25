# Contact matrix implied by a HERMES population

Computes the age-by-age contact matrix that
[HERMES](https://github.com/lwillem/ibm_hermes) effectively uses for
transmission, from the household, school-class and workplace structure
of a simulated population.

## Usage

``` r
contact_matrix_hermes(
  demographic,
  params,
  age_breaks = default_age_breaks,
  age_labels = default_age_labels
)
```

## Arguments

- demographic:

  Data frame with `age`, `hh_id`, `classroom_id` and `workplace_id`
  (`NA` when not in school or work), e.g. `ve_sim_waning$demographic`.

- params:

  HERMES parameter list with `contact_prob_household`,
  `contact_prob_school`, `contact_prob_workplace` and
  `num_contacts_community_day`, e.g. `ve_sim_waning$params`.

- age_breaks, age_labels:

  Passed to [`cut()`](https://rdrr.io/r/base/cut.html) (left-closed
  intervals).

## Value

A square numeric matrix with dimnames `age_group` (rows) and
`contact_age_group` (columns).

## Details

Each day, HERMES gives every infectious person \\i\\ and susceptible
person \\j\\ a contact weight \$\$w\_{ij} = c\_\text{community} +
c\_\text{household} 1\[\text{same household}\] + c\_\text{school}
1\[\text{same class}\] + c\_\text{workplace} 1\[\text{same
workplace}\],\$\$ with \\c\_\text{community} = 1 - \exp(-k/N)\\ for `k`
community contacts per day, and infection probability
`transmission_prob` \\\times w\_{ij}\\. This function returns \\C\_{a'a}
= \frac{1}{N\_{a'}} \sum\_{i \in a'} \sum\_{j \in a, j \ne i} w\_{ij}\\,
the mean daily weighted contacts a person in group \\a'\\ has with group
\\a\\, in the same orientation as
[polymod_nl](https://kylieainslie.github.io/surviver/reference/polymod_nl.md).
It is reciprocal: \\N\_{a'} C\_{a'a} = N_a C\_{aa'}\\.

Averaging over each age group ignores clustering: a household member of
a case is at much higher risk than the age-group average suggests.

## Examples

``` r
cm <- contact_matrix_hermes(ve_sim_waning$demographic, ve_sim_waning$params)
round(cm, 2)
#>          contact_age_group
#> age_group  0-4  5-17 18-49 50-64  65+
#>     0-4   5.42  0.91  2.91  0.77 0.79
#>     5-17  0.28 10.42  2.85  0.83 0.79
#>     18-49 0.38  1.24  3.29  1.14 0.81
#>     50-64 0.21  0.73  2.29  1.67 0.88
#>     65+   0.21  0.68  1.60  0.87 1.40
```
