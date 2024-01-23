# Tests for wb()

mod <- fungible::simFA(Seed = 42)
Omega <- mod$Rpop

set.seed(42)
Sigma <- wb(mod, target_rmsea = 0.05)$Sigma

# Load an indefinite matrix (BadRBY) and create a non-symmetric matrix
data(BadRBY, package = "fungible")
nonsymmetric_matrix <- fungible::rcor(5)
nonsymmetric_matrix[1,2] <- .1

test_that("Errors are thrown when invalid target RMSEA values are given", {
  expect_error(wb(mod = mod, target_rmsea = "a"))
  expect_error(wb(mod = mod, target_rmsea = -.01))
  expect_error(wb(mod = mod, target_rmsea = 1.01))
  expect_error(wb(mod = mod, target_rmsea = NULL))
}
)

test_that("Errors are thrown when invalid mod values are given", {
  expect_error(wb(mod = "a", target_rmsea = 0.05))
  expect_error(wb(mod = list("a" = 1, "b" = 2, "c" = 3), target_rmsea = 0.05))
  expect_error(wb(mod = NULL, target_rmsea = 0.05))
})

test_that("Function output has the expected dimension and type", {
  expect_equal(dim(Omega), dim(Sigma))
  expect_false(any(eigen(Sigma)$values < 0))
  expect_false(any(diag(Sigma) != 1))
  expect_false(any(abs(Sigma) > 1))
})

test_that("Function works when wb_mod is specified.", {
  wb_mod <- get_wb_mod(mod)
  expect_error(wb(mod, target_rmsea = 0.05, wb_mod = -0.1))
  expect_lte(abs(rmsea(wb(mod, target_rmsea = 0.05, wb_mod = wb_mod)$Sigma,
                       Omega, k = ncol(mod$loadings)) - 0.05), 0.01)
})

test_that("Function works when target_rmsea value is large", {
  expect_lte(abs(rmsea(wb(mod, target_rmsea = 0.1)$Sigma,
                       Omega, k = ncol(mod$loadings)) - .1),
             0.03)
})

test_that("An error is thrown if RMSEA is too large.", {
  expect_error(wb(mod, target_rmsea = .35, adjust_target = FALSE))
})
