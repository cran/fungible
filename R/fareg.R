#' Regularized Factor Analysis 
#'
#' This function applies the regularized factoring method to extract an unrotated factor structure matrix.
#'
#' @param R (Matrix) A correlation matrix to be analyzed.
#' @param numFactors (Integer) The number of factors to extract. Default: numFactors = 1.
#' @param facMethod (Character) "rls" for regularized least squares estimation or 
#'         "rml" for regularized maximum likelihood estimation. Default: facMethod = "rls".
#' @return The main output is the matrix of unrotated factor loadings.
#' \itemize{
#'   \item \strong{loadings}: (Matrix) A matrix of unrotated factor loadings.
#'   \item \strong{h2}: (Vector) A vector of estimated communality values.
#'   \item \strong{L}: (Numeric) Value of the estimated penality parameter.
#'   \item \strong{Heywood} (Logical) TRUE if a Heywood case is detected 
#'                (this should never happen).
#'   }
#'
#' @author
#'  Niels G. Waller (nwaller@umn.edu)
#'  
#' @family Factor Analysis Routines
#' 
#' @references Jung, S. & Takane, Y.  (2008).  Regularized common factor analysis.  
#' New trends in psychometrics, 141-149.  
#'
#' @examples
#' 
#' data("HW")
#'
#' # load first HW data set
#' 
#' RHW <- cor(x = HW$HW6)
#' 
#' # Compute principal axis factor analysis
#' fapaOut <- faMain(R = RHW, 
#'                  numFactors = 3, 
#'                  facMethod = "fapa", 
#'                  rotate = "oblimin",
#'                  faControl = list(treatHeywood = FALSE))
#'
#'
#' fapaOut$faFit$Heywood
#' round(fapaOut$h2, 2)
#'
#'  # Conduct a regularized factor analysis
#' regOut <- fareg(R = RHW, 
#'                numFactors = 3,
#'                facMethod = "rls")
#' regOut$L
#' regOut$Heywood
#'
#'
#' # rotate regularized loadings and align with 
#' # population structure
#' regOutRot <- faMain(urLoadings = regOut$loadings,
#'                    rotate = "oblimin")
#'
#' # ALign
#' FHW  <- faAlign(HW$popLoadings, fapaOut$loadings)$F2
#' Freg <- faAlign(HW$popLoadings, regOutRot$loadings)$F2
#'
#' AllSolutions <- round(cbind(HW$popLoadings, Freg, FHW),2) 
#' colnames(AllSolutions) <- c("F1", "F2", "F3", "Fr1", "Fr2", "Fr3", 
#'                            "Fhw1", "Fhw2", "Fhw3")
#' AllSolutions
#' 
#' 
#' rmsdHW <- rmsd(HW$popLoadings, FHW, 
#'               IncludeDiag = FALSE, 
#'               Symmetric = FALSE)
#'
#' rmsdReg <- rmsd(HW$popLoadings, Freg, 
#'                IncludeDiag = FALSE, 
#'                Symmetric = FALSE)
#'
#' cat("\nrmsd HW =  ", round(rmsdHW,3),
#'     "\nrmsd reg = ", round(rmsdReg,3))
#'
#'
#'
#' @export

#--------------------------------------------------------#
fareg <- function(R, 
                  numFactors = 1,
                  facMethod = "rls"){
  
  
  # Date: 04/19/2019
  # Arguments:
  #
  # R  (numeric matrix) Input correlation matrix.
  # numfactors (scaler) Number of factors to extract, Default = 1
  # facMethod (Character) "rls" for regulaized least squares estimation or
  #                       "rml" for regularixed maximum likelihood estimation.
  #Values
  # loadings,
  # h2,
  # L,
  # Heywood
  
  if(!isSymmetric(R)) stop("\nInput data must be a correlation matrix")

  if( !all.equal( as.vector(diag(R)), rep(1, dim(R)[1]) ) ){
     stop("\nInput data must be a correlation matrix")
  }  
  
  PlausibleMethods <- c("rls", "rml")
  if(facMethod %in% PlausibleMethods == FALSE){
    stop("\nInvalid argument for facMethod")
  }
  p <- nrow(R)
  
  
  communality <- "SMC"  
  ## Check to ensure that R is positive definite 
  eigenVal <- eigen(R)$value
  if ( min(eigenVal) <= 1E-6 ) {
    warning("Inverting the correlation matrix for SMC communality estimates requires a positive-definite matrix.")
    communality <- "maxr"  
  } # END if ( min(eigenVal) <= 1E-8 )
  
  
  if(communality == "SMC"){
     Psi <- diag(diag(solve(R))^-1)
  }
  if(communality == "maxr"){
    R0 <- R
    diag(R0) <- 0
    h2temp <- apply(abs(R0),2, max)
    Psi <- diag(1 - h2temp)
  }  
  
    #regularization parameter
    Lstart = 1

    # function to minimize
  FncRidgeULS <- function(L){
      Ri <- R - (L * Psi)
      UDU <- eigen(Ri)
      U <- UDU$vectors
      D <- UDU$values
      
      if( numFactors == 1){
        F <- U[,1] * sqrt(D[1])
      }
      
      if( numFactors > 1){
        F <- U[,1:numFactors] %*% diag(sqrt(D[1:numFactors]))
      }
      
      
      Rhat <- F %*% t(F)
    
      
      #return fnc
      rmsd(Ri, Rhat, 
           Symmetric = TRUE,
           IncludeDiag = FALSE)
  }
  
  # function to minimize
  FncRidgeMLE <- function(L){
    Ri <- R - (L * Psi)
    UDU <- eigen(Ri)
    U <- UDU$vectors
    D <- UDU$values
    
    if( numFactors == 1){
      F <- U[,1] * sqrt(D[1])
    }
    
    if( numFactors > 1){
      F <- U[,1:numFactors] %*% diag(sqrt(D[1:numFactors]))
    }
    
    
    Rhat <- F %*% t(F)
    diag(Rhat) <- 1
    RhatInv <- solve(Rhat)
    log(det(Rhat)) + sum(diag(R %*% RhatInv )) - log(det(R)) - p
  }
  
  
  # function to produce loadings from optimal L
  faSolution <- function(L){
    Ri <- R - (L * Psi)
    UDU <- eigen(Ri)
    U <- UDU$vectors
    D <- UDU$values
    
    # return unrotated factor loading matrix
    if( numFactors == 1){
      F <- U[,1] * sqrt(D[1])
    }
    
    if( numFactors > 1){
      F <- U[,1:numFactors] %*% diag(sqrt(D[1:numFactors]))
    }

    F
  }
  
  maxPsi <- max(diag(Psi))
  Lup <- 1/maxPsi + .001
 
  
  # optimize function
  if(facMethod == "rls"){
     out <- optimize(f = FncRidgeULS, lower = .001, upper = Lup)
  }
  if(facMethod == "rml"){
    out <- optimize(f = FncRidgeMLE, lower = 0, upper = Lup)
  }
  
  
  Lopt <- out$minimum
  
  # set convergence status
  converged <- FALSE
  if(Lopt > 0 & Lopt < Lup) converged <- TRUE
  
  # Compute regularized loadings
  loadings <- faSolution(Lopt)
  
  if(numFactors == 1){
    loadings <- as.matrix(loadings)
    h2 <- as.vector(loadings^2)
  }  
  
  if( numFactors > 1){
      h2 <- apply(loadings^2, 1, sum)
  }    
  
  
  Heywood <- FALSE
  # By design, Heywood should never be TRUE
  if(max(h2 > 1)) Heywood <- TRUE
  
  list(loadings = loadings,
       h2 = h2,
       L = Lopt,
       converged = converged,
       Heywood = Heywood)
}

