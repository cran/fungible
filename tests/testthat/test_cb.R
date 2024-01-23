# Tests for cb()

library(noisemaker)

mod <- fungible::simFA(Seed = 42)
set.seed(42)
Sigma <- cb(mod, target_rmsea = 0.05)
Omega <- mod$Rpop

test_that("Errors are thrown when invalid target RMSEA values are given", {
  expect_error(cb(mod, target_rmsea = "a"))
  expect_error(cb(mod, target_rmsea = -.01))
  expect_error(cb(mod, target_rmsea = 1.01))
  expect_error(cb(mod, target_rmsea = NULL))
}
)

test_that("Errors are thrown when invalid mod values are given", {
  expect_error(cb(mod = "a", target_rmsea = 0.05))
  expect_error(cb(mod = list(a = 1, b = 2, c = 3), target_rmsea = 0.05))
  expect_error(cb(mod = NULL, target_rmsea = 0.05))
})

test_that("Function output has the expected dimension and type", {
  Omega <- mod$Rpop
  expect_equal(dim(Omega), dim(Sigma))
  expect_false(any(eigen(Sigma)$values < 0))
  expect_false(any(diag(Sigma) != 1))
  expect_false(any(abs(Sigma) > 1))
})

test_that("Error is thrown if Sigma is indefinite", {
  expect_error(cb(mod, target_rmsea = .5))
})

test_that("RMSEA value is in the ballpark of the target RMSEA value", {
  expect_true(abs(rmsea(Sigma, Omega, k = 3) - 0.05) < 0.01)
})
