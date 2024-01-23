# Tests for get_wb_mod()

mod <- fungible::simFA(Seed = 42)
set.seed(42)

test_that("Errors are thrown when `mod` isn't a simFA() object.", {
  expect_error(get_wb_mod(mod = "a"))
  expect_error(get_wb_mod(mod = list(a = 1, b = 2)))
}
)

test_that("Errors are thrown when `n` is not valid.", {
  expect_error(get_wb_mod(mod, n = NA))
  expect_error(get_wb_mod(mod, n = -1))
  expect_error(get_wb_mod(mod, n = "a"))
})

test_that("Errors are thrown when `values` is not valid.", {
  expect_error(get_wb_mod(mod, values = NA))
  expect_error(get_wb_mod(mod, values = -1))
  expect_error(get_wb_mod(mod, values = 1))
})

test_that("Errors are thrown when `lower` is not valid.", {
  expect_error(get_wb_mod(mod, lower = NA))
  expect_error(get_wb_mod(mod, lower = -.1))
  expect_error(get_wb_mod(mod, lower = 0))
})

test_that("Errors are thrown when `upper` is not valid.", {
  expect_error(get_wb_mod(mod, upper = NA))
})

test_that("Warning is generated when `upper` is too large.", {
  expect_warning(get_wb_mod(mod, upper = 1))
})

