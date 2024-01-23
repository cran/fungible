## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE---------------------------------------------------
library(fungible)

## ----create-model-------------------------------------------------------------
# Specify the factor model
Lambda <- matrix(c(.5, .5, .5,  0,  0,  0,  0,  0,  0,
                    0,  0,  0, .6, .6, .6,  0,  0,  0,
                    0,  0,  0,  0,  0,  0, .7, .7, .7),
                 ncol = 3, byrow = FALSE)
Phi <- matrix(c( 1, .3, .3,
                .3,  1, .3,
                .3, .3,  1),
              ncol = 3, byrow = TRUE)

mod <- fungible::simFA(
  Model = list(NFac = 3,
               NItemPerFac = 3,
               Model = "oblique"),
  Loadings = list(FacPattern = Lambda),
  Phi = list(PhiType = "user",
             UserPhi = Phi),
  Seed = 42
)

## ----factor-pattern-----------------------------------------------------------
# factor-pattern matrix
mod$loadings

## ----factor-correlation-------------------------------------------------------
# factor correlation matrix
mod$Phi

## -----------------------------------------------------------------------------
mod$Rpop

## ----TKL-rmsea----------------------------------------------------------------
set.seed(42)
TKL_m1 <- noisemaker(mod, method = "TKL", target_rmsea = 0.05)

TKL_m1

## ----TKL-rmsea-and-cfi--------------------------------------------------------
TKL_m2 <- noisemaker(mod, method = "TKL", 
                       target_rmsea = 0.05, 
                       target_cfi = 0.95)

TKL_m2

## ----TKL-rmsea-cfi-weights----------------------------------------------------
TKL_m3 <- noisemaker(mod, method = "TKL", 
                       target_rmsea = 0.05, 
                       target_cfi = 0.95,
                     tkl_ctrl = list(weights = c(rmsea = 4, cfi = 1)))

TKL_m3

## -----------------------------------------------------------------------------
TKL_m4 <- noisemaker(mod, method = "TKL", 
                       target_rmsea = 0.05, 
                       target_cfi = 0.95,
                     tkl_ctrl = list(weights = c(rmsea = 1, cfi = 4)))

TKL_m4

## ----TKL----------------------------------------------------------------------
TKL_m5 <-  noisemaker(mod, method = "TKL", 
                      target_rmsea = 0.05, 
                      target_cfi = 0.95,
                      tkl_ctrl = list(WmaxLoading = 0.3,
                                      NWmaxLoading = 2))

TKL_m5

## ----bounds-on-v-and-eps------------------------------------------------------
TKL_m6 <-  noisemaker(mod, method = "TKL", 
                      target_rmsea = 0.05, 
                      target_cfi = 0.95,
                      tkl_ctrl = list(v_bounds = c(0, .2),
                                      eps_bounds = c(0, 1)))

TKL_m6

## ----cb-example---------------------------------------------------------------
CB_m1 <- noisemaker(mod, method = "CB", target_rmsea = 0.05)
CB_m1

## ----cb-error, error = TRUE---------------------------------------------------
CB_m1 <- noisemaker(mod, method = "CB", target_rmsea = 0.5)

## ----get-wb-mod---------------------------------------------------------------
wb_mod <- get_wb_mod(
  mod,           # simFA() model specification
  n = 50,        # Number of matrices to simulate at each target RMSEA value
  values = 10,   # Number of target RMSEA values to test
  lower = 0.01,  # 'lower' and 'upper' are the endpoints of the RMSEA sequence
  upper = 0.095  
)

summary(wb_mod)

## ----wb-example---------------------------------------------------------------
noisemaker(mod, method = "WB", target_rmsea = 0.05, wb_mod = wb_mod)

