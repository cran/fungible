# Test noisemaker()

library(noisemaker)
library(fungible)

mod <- fungible::simFA(Seed = 42)
set.seed(42)

test_that("Errors are thrown when invalid target RMSEA or CFI values are given", {
  expect_error(noisemaker(mod, method = "TKL", target_rmsea = "a"))
  expect_error(noisemaker(mod, method = "CB", target_rmsea = "a"))
  expect_error(noisemaker(mod, method = "WB", target_rmsea = "a"))

  expect_error(noisemaker(mod, method = "TKL", target_rmsea = -.01))
  expect_error(noisemaker(mod, method = "CB",  target_rmsea = -.01))
  expect_error(noisemaker(mod, method = "WB",  target_rmsea = -.01))

  expect_error(noisemaker(mod, method = "TKL", target_rmsea = 1.01))
  expect_error(noisemaker(mod, method = "CB",  target_rmsea = 1.01))
  expect_error(noisemaker(mod, method = "WB",  target_rmsea = 1.01))

  expect_error(noisemaker(mod, method = "TKL", target_rmsea = NULL))
  expect_error(noisemaker(mod, method = "CB",  target_rmsea = NULL))
  expect_error(noisemaker(mod, method = "WB",  target_rmsea = NULL))

  expect_error(noisemaker(mod, method = "TKL", target_cfi = 1.01))
  expect_error(noisemaker(mod, method = "CB",  target_cfi = 0.95))
  expect_error(noisemaker(mod, method = "WB",  target_cfi = 0.95))

  expect_error(noisemaker(mod, method = "TKL", target_cfi = "a"))
  expect_error(noisemaker(mod, method = "CB",  target_cfi = "a"))
  expect_error(noisemaker(mod, method = "WB",  target_cfi = "a"))
}
)

test_that("An error is thrown when an invalid method is given", {
  expect_error(noisemaker(mod, method = "AB"))
})

test_that("An error is thrown if an invalid `simFA()` object is given", {
  mod2 <- list(a = 1, b = 2, c = 3)
  expect_error(noisemaker(mod2, method = "WB"))
})

test_that("A warning is given if a `simFA()` model is provided that includes model error", {
  mod3 <- simFA(ModelError = list(ModelError = TRUE))
  expect_warning(noisemaker(mod3, method = "WB"))
  expect_warning(noisemaker(mod3, method = "CB"))
  expect_warning(noisemaker(mod3, method = "TKL"))
})
