test_that("as_direct_model.data.frame with default settings an object with the correct classes", { # nolint: line_length_linter.
  prep_obs <- as_direct_model(sim_obs)
  expect_s3_class(prep_obs, "data.frame")
  expect_s3_class(prep_obs, "epidist_direct_model")
})

test_that("as_direct_model.data.frame errors when passed incorrect inputs", { # nolint: line_length_linter.
  expect_error(as_direct_model(list()))
  expect_error(as_direct_model(sim_obs[, 1]))
  expect_error({
    sim_obs$case <- paste("case_", seq_len(nrow(sim_obs)))
    as_direct_model(sim_obs)
  })
})

# Make this data available for other tests
prep_obs <- as_direct_model(sim_obs)
family_lognormal <- epidist_family(prep_obs, family = brms::lognormal())

test_that("is_direct_model returns TRUE for correct input", { # nolint: line_length_linter.
  expect_true(is_direct_model(prep_obs))
  expect_true({
    x <- list()
    class(x) <- "epidist_direct_model"
    is_direct_model(x)
  })
})

test_that("is_direct_model returns FALSE for incorrect input", { # nolint: line_length_linter.
  expect_false(is_direct_model(list()))
  expect_false({
    x <- list()
    class(x) <- "epidist_direct_model_extension"
    is_direct_model(x)
  })
})

test_that("epidist_validate.epidist_direct_model doesn't produce an error for correct input", { # nolint: line_length_linter.
  expect_no_error(epidist_validate(prep_obs))
})

test_that("epidist_validate.epidist_direct_model returns FALSE for incorrect input", { # nolint: line_length_linter.
  expect_error(epidist_validate(list()))
  expect_error(epidist_validate(prep_obs[, 1]))
  expect_error({
    x <- list()
    class(x) <- "epidist_direct_model"
    epidist_validate(x)
  })
})
