# Test rmsea()

library(noisemaker)

mod <- fungible::simFA(ModelError = list(ModelError = TRUE),
                       Seed = 42)
set.seed(42)
X <- fungible::rcor(5)
Y <- fungible::rcor(6)
Z <- Y
Z[1,2] <- 1

test_that("RMSEA value agrees with the RMSEA value from `simFA()`", {
  expect_equal(rmsea(mod$RpopME, mod$Rpop, k = 3),
               mod$ModelErrorFitStats$RMSEA_theta)
})

test_that("The function throws an error if the arguments aren't matrices or if they have different dimensions.", {
  expect_error(rmsea("a", "b", k = 3))
  expect_error(rmsea(X, Y, k = 3))
  expect_error(rmsea(Y, Z, k = 3))
  expect_error(rmsea(Z, Y, k = 3))
})

test_that("The function throws an error if k is not a non-negative integer", {
  expect_error(rmsea(X, X, k = 1.2))
  expect_error(rmsea(X, X, k = "a"))
  expect_error(rmsea(X, X, k = -0.1))
})
