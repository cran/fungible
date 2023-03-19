#' Complete a Partially Specified Correlation Matrix by the Method of Differential Evolution
#'
#' This function completes a partially specified
#' correlation matrix by the method of differential evolution.
#'
#' @param  Rna  (matrix) An n x n incomplete correlation matrix.  Missing entries must
#'     be specified by \code{NA} values.
#' @param NMatrices (integer) \code{CompleteRDEV} will complete \code{NMatrices}
#' correlation matrices.
#' @param  MaxDet (logical) If MaxDet = TRUE then the correlation matrix will
#'     be completed with entries that maximize the determinant of R.
#' @param  MaxIter (integer) The maximum number of iterations
#'       (i.e., generations) allowed. Default \code{MaxIter = 200}.
#' @param  delta (numeric > 0) A number that controls the convergence
#'     accuracy of the differential evolution algorithm. Default \code{delta = 1E-8}.
#' @param  PRINT (logical) When PRINT = TRUE the algorithm convergence status is printed.
#'     Default  \code{PRINT = FALSE}.
#' @param  Seed (integer) Initial random number seed. Default (\code{Seed = NULL}).
#'
#' @return \code{CompleteRdev} returns the following objects:
#' \itemize{
#'    \item  \strong{R} (matrix) A PSD completed correlation matrix.
#'    \item  \strong{converged}: (logical) a logical that indicates the convergence status of the optimizaton.
#'    \item \strong{iter} (integer) The number of cycles needed to reach converged solution.
#' }
#'
#'
#' @references Ardia, D., Boudt, K., Carl, P., Mullen, K.M., Peterson, B.G. (2011) Differential
#' Evolution with DEoptim. An Application to Non-Convex Portfolio Optimization.
#' URL The R Journal, 3(1), 27-34.
#' URL https://journal.r-project.org/archive/2011-1/RJournal_2011-1_Ardia~et~al.pdf.
#'
#' @references Georgescu, D. I., Higham, N. J., and Peters, G. W.  (2018).  Explicit
#' solutions to correlation matrix completion problems, with
#' an application to risk management and insurance.  Royal Society Open
#' Science, 5(3), 172348.
#'
#' @references Mauro, R.  (1990).  Understanding L.O.V.E. (left out variables
#' error): a method for estimating the effects of omitted variables.
#' Psychological Bulletin, 108(2), 314-329.
#
#' @references Mishra, S. K.  (2007).  Completing correlation matrices
#' of arbitrary order by differential evolution method of global optimization:
#' a Fortran program.  Available at SSRN 968373.
#'
#' @references Mullen, K.M, Ardia, D., Gil, D., Windover, D., Cline, J. (2011). DEoptim: An
#' R Package for Global Optimization by Differential Evolution. Journal of Statistical Software, 40(6), 1-26. URL http://www.jstatsoft.org/v40/i06/.
#'
#' @references Price, K.V., Storn, R.M., Lampinen J.A. (2005) Differential Evolution - A Practical Approach to Global Optimization. Berlin Heidelberg:
#' Springer-Verlag. ISBN 3540209506.
#'
#' @references  Zhang, J. and Sanderson, A. (2009) Adaptive Differential
#' Evolution Springer-Verlag. ISBN 978-3-642-01526-7
#'
#' @author Niels G. Waller
#'
#' @examples
#' ## Example 1: Generate random 4 x 4 Correlation matrices.
#'   Rmiss <- matrix(NA, nrow = 4, ncol = 4)
#'   diag(Rmiss) <- 1
#'
#'   out <- CompleteRdev(Rna = Rmiss,
#'                       NMatrices = 4,
#'                       PRINT = TRUE,
#'                       Seed = 1)
#'
#'   print( round( out$R[[1]] , 3) )
#'
#'\dontrun{
#' # Example 2: Complete a partially specified R matrix.
#' # Example from Georgescu, D. I., Higham, N. J., and
#' #              Peters, G. W.  (2018).
#'
#' Rmiss <- matrix(
#'      c( 1,  .25, .6,  .55, .65,  0,  .4,   .6,  .2,  .3,
#'        .25, 1,    0,   0,   0,   0,  NA,   NA,  NA,  NA,
#'        .6,  0,   1,   .75, .75,  0,  NA,   NA,  NA,  NA,
#'        .55, 0,   .75, 1,   .5,   0,  NA,   NA,  NA,  NA,
#'        .65, 0,   .75,  .5, 1,    0,  NA,   NA,  NA,  NA,
#'         0,  0,    0,   0,   0,  1,   NA,   NA,  NA,  NA,
#'         .4, NA,   NA,  NA,  NA,  NA, 1,   .25, .25,  .5,
#'         .6, NA,   NA,  NA,  NA,  NA, .25,  1,  .25,  0,
#'         .2, NA,   NA,  NA,  NA,  NA, .25,  .25, 1,   0,
#'         .3, NA,   NA,  NA,  NA,  NA, .5,    0,   0,  1), 10,10)
#'
#' # Complete Rmiss with values that maximize
#' # the matrix determinant (this is the MLE solution)
#'  set.seed(123)
#'  out <- CompleteRdev(Rna = Rmiss,
#'                      MaxDet = TRUE,
#'                      MaxIter = 1000,
#'                      delta = 1E-8,
#'                      PRINT = FALSE)
#'
#' cat("\nConverged = ", out$converged,"\n")
#' print( round(out$R, 3))
#' print( det(out$R))
#' print( eigen(out$R)$values, digits = 5)
#'}
#'
#' @import DEoptim
#' @export
CompleteRdev <- function(Rna,
                         NMatrices = 1,
                         MaxDet = FALSE,
                         MaxIter = 200,
                         delta = 1E-8,
                         PRINT = FALSE,
                         Seed = NULL){


  # MaxDet is unique
  if(MaxDet) NMatrices = 1
  p <- nrow(Rna)

  #
  Rlist <- list(matrix(NA, p,p))
  ConvergeVec <- iterVec <- rep(NA, NMatrices)


  ## Generate random seed if not supplied
  if( is.null(Seed) ){
    Seed <- as.integer((as.double(Sys.time()) * 1000 + Sys.getpid()) %% 2^31)
  } # END if (is.null(Seed))
  set.seed(Seed)

  # Imat = Identity matrix
  Imat <- diag(nrow(Rna))
  # Grab rij from upper triangle of Rna
  UpRna <- Rna[upper.tri(Rna, diag=FALSE)]
  # Count the number of missing rij values
  num_miss <- sum(is.na(UpRna))


  # After optimization insert final values in R matrix ---
  InsertFinalValues <- function(Rna, X){
     #Step 1 replace missing r's with final values
      UpRna <- Rna[upper.tri(Rna, diag=FALSE)]
      UpRna[is.na(UpRna)]<- X
      Ri <- Imat
      Ri[upper.tri(Ri, diag=FALSE)] <- UpRna
      Ri <- Ri + t(Ri) - Imat
      Ri
   }#END InsertFinalValues


  # Constraint: R is PSD
  EvoFunc <- function(x){
     UpRna <- Rna[upper.tri(Rna, diag=FALSE)]
     UpRna[is.na(UpRna)]<- x
     Ri <- Imat
     Ri[upper.tri(Ri, diag=FALSE)] <- UpRna
     Ri <- Ri + t(Ri) - Imat
     L <- eigen(Ri, only.values=TRUE)$values
     Lneg <- L[L < 0]
     if( length(Lneg) > 0){
        #penalty to insure PSD
        ans <- 500 * sum(abs(Lneg))
     }else ans <- 0
     ans
  }#END EvoFunc

  # Complete R by Evolutionary Algorithm
  # find completed R with max(det(R_full))
  MaxDetFunc <- function(x){
    UpRna <- Rna[upper.tri(Rna, diag=FALSE)]
    UpRna[is.na(UpRna)]<- x
    Ri <- Imat
    Ri[upper.tri(Ri, diag=FALSE)] <- UpRna
    Ri <- Ri + t(Ri) - Imat
    L <- eigen(Ri, only.values=TRUE)$values
    Lneg <- L[L < 0]
    # constraint
    if( length(Lneg) > 0){
      ans <-  -log(abs(det(Ri))) + 500*(sum(abs(Lneg)))
    }else ans <- -log(abs(det(Ri)))
    ans
  }#END MaxDetFunc


 # Run Optimization --------
 if(MaxDet){
   #find RFull that maximizes det(RFull)
   # Remember: DEoptim minimizes a function
          dout <<- DEoptim::DEoptim(fn = MaxDetFunc,
                           lower = rep(-1,num_miss),
                           upper = rep(1,num_miss),
                           control = list(trace = FALSE,
                                          itermax = MaxIter))
    # assess convergence
          converged = FALSE
          Delta <-   abs(dout$member$bestval[MaxIter] - dout$member$bestval[MaxIter-5])
          if(Delta <= delta)  converged <- TRUE

    # Save output
          Rlist <- InsertFinalValues(Rna, dout$optim$bestmem)
          ConvergeVec <- converged
          iterVec <- dout$optim$iter

          # Print convergence status
          if(PRINT){
            cat("\nConverged = ", converged, "\n")
          }
    }#END MaxDet
  else{
   # else find NMatrices feasible solutions

    k = 1
    while( k <= NMatrices){

       dout <- tryCatch(DEoptim::DEoptim(fn = EvoFunc,
                               lower = rep(-1,num_miss),
                               upper = rep(1,num_miss),
                               control = list(trace = FALSE,
                                              VTR = 0,
                                              itermax = MaxIter)),
                        error = " ") # do nothing on error
     # assess convergence
       converged = FALSE
       if(round(dout$optim$bestval,8) == 0) converged <- TRUE


       if(converged){
         if(PRINT) cat(" ", k)
         Rlist[[k]] <- InsertFinalValues(Rna, dout$optim$bestmem)
         ConvergeVec[k] <- converged
         iterVec[k] <- dout$optim$iter
         k  <- k + 1
       } #END if converged

    } #END while loop

  }#END if MaxDet



  # RETURN LIST
  list(R = Rlist,
       converged = ConvergeVec,
       iter = iterVec)
} #END CompleteRdev


