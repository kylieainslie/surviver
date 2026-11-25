## POLYMOD contact matrix for the Netherlands, aggregated to surviver's
## default age groups. Originally computed in the Lorentz workshop repository
## (analysis/group-D/prepare_stan_data.R); reproduced exactly with
## socialmixr 0.6.x.
##
## socialmixr >= 0.6.0 warns that the automatic country-population lookup
## used here is deprecated; a future version will need `survey_pop` supplied.

age_labels <- c("0-4", "5-17", "18-49", "50-64", "65+")  # = default_age_labels

polymod_nl <- socialmixr::contact_matrix(
  socialmixr::polymod,
  countries             = "Netherlands",
  age.limits            = c(0, 5, 18, 50, 65),
  symmetric             = TRUE,
  estimated.contact.age = "mean",
  missing.contact.age   = "remove"
)$matrix

dimnames(polymod_nl) <- list(
  age_group         = age_labels,
  contact_age_group = age_labels
)

usethis::use_data(polymod_nl, overwrite = TRUE)
