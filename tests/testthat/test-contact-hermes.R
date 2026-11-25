test_that("contact_matrix_hermes matches a person-by-person calculation", {
  d <- data.frame(
    age          = c(2, 30, 32, 10, 12, 70, 45),
    hh_id        = c(1, 1, 1, 2, 2, 3, 4),
    classroom_id = c("k1", NA, NA, "c1", "c1", NA, NA),
    workplace_id = c(NA, 10, 11, NA, NA, NA, 10)
  )
  p <- list(contact_prob_household = 0.9, contact_prob_school = 0.5,
            contact_prob_workplace = 0.1, num_contacts_community_day = 4)
  cm <- contact_matrix_hermes(d, p)

  same <- function(x) {
    o <- outer(as.character(x), as.character(x), "==")
    o[is.na(o)] <- FALSE
    o
  }
  w <- (1 - exp(-4 / nrow(d))) + 0.9 * same(d$hh_id) + 0.5 * same(d$classroom_id) +
    0.1 * same(d$workplace_id)
  diag(w) <- 0
  g <- as.integer(add_age_group(d)$age_group)
  expected <- matrix(NA_real_, 5, 5)
  for (ap in 1:5) for (a in 1:5) {
    expected[ap, a] <- if (any(g == ap)) sum(w[g == ap, g == a]) / sum(g == ap) else NaN
  }
  expect_equal(unname(cm), expected)
  expect_equal(rownames(cm), default_age_labels)
})

test_that("contact_matrix_hermes is reciprocal on simulated data", {
  cm <- contact_matrix_hermes(ve_sim_waning$demographic, ve_sim_waning$params)
  n  <- age_group_sizes(ve_population(ve_sim_waning))
  expect_equal(unname(n * cm), unname(t(n * cm)))
  expect_true(all(cm > 0))
})
