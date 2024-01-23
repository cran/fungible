# Test tkl()

library(noisemaker)
library(fungible)

mod <- fungible::simFA()
set.seed(42)
Sigma <- tkl(mod, target_rmsea = 0.05, target_cfi = 0.95)

test_that("Errors are thrown when invalid target RMSEA or CFI values are given", {
  expect_error(tkl(mod, target_rmsea = "a"))
  expect_error(tkl(mod, target_rmsea = -.01))
  expect_error(tkl(mod, target_rmsea = 1.01))
  expect_error(tkl(mod, target_rmsea = NULL))
  expect_error(tkl(mod, target_cfi = NULL))
  expect_error(tkl(mod, target_cfi = 1.01))
  expect_error(tkl(mod, target_cfi = NULL, target_rmsea = NULL))
}
)

test_that("Errors are thrown when invalid tkl_ctrl arguments are given", {
  expect_error(tkl(mod, target_rmsea = .05,
                   tkl_ctrl = list(weights = c("a", 2))))
  expect_error(tkl(mod, target_rmsea = .05,
                   tkl_ctrl = list(v_start = 2)))
  expect_error(tkl(mod, target_rmsea = .05,
                   tkl_ctrl = list(eps_start = 2)))
  expect_error(tkl(mod, target_rmsea = .05,
                   tkl_ctrl = list(NMinorFac = -1)))
  expect_error(tkl(mod, target_rmsea = .05,
                   tkl_ctrl = list(optim_type = "ag")))
  expect_error(tkl(mod, target_rmsea = .05,
                   tkl_ctrl = list(penalty = -1)))
  expect_error(tkl(mod, target_rmsea = .05,
                   tkl_ctrl = list(WmaxLoading = -1)))
  expect_error(tkl(mod, target_rmsea = .05,
                   tkl_ctrl = list(NWmaxLoading = -1)))
})

test_that("Errors are thrown when invalid mod values are given", {
  expect_error(tkl(mod = "a", target_rmsea = 0.05, mod = "a"))
  expect_error(tkl(mod = list(a = 1, b = "test", c = 3), target_rmsea = 0.05))
  expect_error(tkl(mod = NULL, target_rmsea = 0.05))
})

test_that("Function works when only RMSEA is specified and `ga()` is used.", {
  sol <- tkl(mod = mod, target_rmsea = 0.05, tkl_ctrl = list(optim_type = "ga"))
  expect_equal(sol$Rpop, mod$Rpop)
  expect_equal(dim(sol$RpopME), dim(mod$Rpop))
  expect_true(abs(sol$rmsea - 0.05) < 0.01)
})

test_that("Function works when only CFI is specified and `ga()` is used.", {
  sol <- tkl(mod = mod, target_cfi = 0.95, tkl_ctrl = list(optim_type = "ga"))
  expect_equal(sol$Rpop, mod$Rpop)
  expect_equal(dim(sol$RpopME), dim(mod$Rpop))
  expect_true(abs(sol$cfi - 0.95) < 0.01)
})

test_that("Function works when CFI and RMSEA are specified and `ga()` is used.", {
  sol <- tkl(mod = mod, target_cfi = 0.95, target_rmsea = 0.05,
             tkl_ctrl = list(optim_type = "ga"))
  expect_equal(sol$Rpop, mod$Rpop)
  expect_equal(dim(sol$RpopME), dim(mod$Rpop))
  expect_true(abs(sol$cfi - 0.95) < 0.05)
  expect_true(abs(sol$rmsea - 0.05) < 0.05)
})

test_that("Function works when only RMSEA is specified and `optim()` is used.", {
  sol <- tkl(mod = mod, target_rmsea = 0.05, tkl_ctrl = list(optim_type = "optim"))
  expect_equal(sol$Rpop, mod$Rpop)
  expect_equal(dim(sol$RpopME), dim(mod$Rpop))
  expect_true(abs(sol$rmsea - 0.05) < 0.01)
})

test_that("Function works when only CFI is specified and `optim()` is used.", {
  sol <- tkl(mod = mod, target_cfi = 0.95,
             tkl_ctrl = list(optim_type = "optim"))
  expect_equal(sol$Rpop, mod$Rpop)
  expect_equal(dim(sol$RpopME), dim(mod$Rpop))
  expect_true(abs(sol$cfi - 0.95) < 0.01)
})

test_that("Function works when CFI and RMSEA are specified and `optim()` is used.", {
  sol <- tkl(mod = mod, target_cfi = 0.95, target_rmsea = 0.05,
             tkl_ctrl = list(optim_type = "optim"))
  expect_equal(sol$Rpop, mod$Rpop)
  expect_equal(dim(sol$RpopME), dim(mod$Rpop))
  expect_true(abs(sol$cfi - 0.95) < 0.05)
  expect_true(abs(sol$rmsea - 0.05) < 0.05)
})

test_that("The bounds on v and eps are working as intended.", {
  sol <- tkl(mod = mod, target_cfi = 0.95, target_rmsea = 0.05,
             tkl_ctrl = list(v_bounds = c(0, .2),
                             eps_bounds = c(0, .05)))
  expect_true(sol$v <= .2)
  expect_true(sol$eps <= .05)
})

test_that("Function works when CFI and RMSEA are specified and `ga()` is used with NWmaxLoading constraints.", {
  sol <- tkl(mod = mod, target_cfi = 0.95, target_rmsea = 0.05,
             tkl_ctrl = list(optim_type = "ga",
                             WmaxLoading = .3,
                             NWmaxLoading = 3))
  expect_equal(sol$Rpop, mod$Rpop)
  expect_equal(dim(sol$RpopME), dim(mod$Rpop))
  expect_true(abs(sol$cfi - 0.95) < 0.05)
  expect_true(abs(sol$rmsea - 0.05) < 0.05)
  expect_false(any(apply(sol$W, 2, function(x) sum(x >= .3)) > 3))
})
