# Test cfi()

library(noisemaker)
mod <- fungible::simFA(ModelError = list(ModelError = TRUE),
                       Seed = 42)
set.seed(42)
X <- fungible::rcor(5)
Y <- fungible::rcor(6)
Z <- Y
Z[1,2] <- 1

test_that("CFI value agrees with the CFI value from `simFA()`", {
  expect_equal(cfi(mod$RpopME, mod$Rpop),
               mod$ModelErrorFitStats$CFI_theta)
})

test_that("The function throws an error if the arguments aren't matrices or if they have different dimensions.", {
  expect_error(cfi("a", "b"))
  expect_error(cfi(X, Y))
  expect_error(cfi(Y, Z))
  expect_error(cfi(Z, Y))
})
