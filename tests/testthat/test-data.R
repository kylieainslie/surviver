test_that("simulated VE datasets have the expected structure", {
  for (d in list(ve_sim_waning, ve_sim_no_waning)) {
    expect_named(d, c("demographic", "vaccination", "incidence", "params",
                      "study_end", "hermes_commit"))
    expect_equal(nrow(d$demographic), d$params$pop_size)
    expect_false(anyDuplicated(d$demographic$ind_id) > 0)
    expect_true(all(d$vaccination$ind_id %in% d$demographic$ind_id))
    expect_true(all(d$incidence$ind_id %in% d$demographic$ind_id))
    expect_true(all(d$vaccination$time_of_vaccination <= d$study_end))
  }
})

test_that("scenarios differ only in the VE floor", {
  p_w  <- ve_sim_waning$params
  p_nw <- ve_sim_no_waning$params
  expect_equal(p_w$ve_floor, 0.1)
  expect_equal(p_nw$ve_floor, p_nw$vaccine_effectiveness)
  p_w$ve_floor <- p_nw$ve_floor <- NULL
  expect_identical(p_w, p_nw)
  expect_identical(ve_sim_waning$demographic, ve_sim_no_waning$demographic)
})
