# Niels Waller
# June 29, 2023
#'
#' Generate random R matrices with user-defined bounds on the correlation 
#' coefficients via differential evolution (DE).
#' 
#' @description
#' \code{Rbounds} can generate uniformly sampled correlation matrices with 
#' user-defined bounds on the correlation coefficients via differential 
#' evolution (DE). Unconstrained \eqn{R} matrices (i.e., with no constraints placed 
#' on the \eqn{r_{ij}}) computed from 12 or fewer variables can be generated relatively 
#' quickly on a personal computer.  Larger matrices may require 
#' very long execution times. \code{Rbounds} can 
#' generate larger matrices when the correlations are tightly 
#' bounded (e.g., \eqn{0 < r_{ij} < .5} for all \eqn{i \neq j}). To generate 
#' uniformly sampled \eqn{R} matrices, users should leave  
#' \code{NPopFactor} and \code{crAdaption} at 
#' their default values.
#' 
#'
#' @param Nvar  (integer) The order of the generated correlation matrices.  
#' @param NMatrices (integer) Generate \code{NMatrices}
#' correlation matrices.
#' @param Minr  (numeric > -1 and < Maxr)  The lower bound for all  \eqn{r_{ij}} in 
#' the generated R matrices.  Default \code{Minr = -1}.   
#' @param Maxr (numeric > Minr and <= 1). The upper bound for all \eqn{r_{ij}} in the 
#' generated \eqn{R} matrices.  Default \code{Maxr = 1}. 
#' @param MinEig (numeric). Minimum size of the last eigenvalue of R. Default 
#' \code{MinEig = 0}. By setting \code{MinEig} to a value slightly greater than 
#' 0 (e.g., 1E-3), all generated matrices will be positive definite.     
#' @param MaxIter (integer) The maximum number of iterations
#'       (i.e., generations) for the DE optimizer. Default \code{MaxIter = 200}.
#' @param NPopFactor (numeric > 0). If \eqn{R} is an \eqn{n \times n} matrix, then each generation
#'  will contain \code{NPopFactor} \eqn{\times n(n-1)/2}  members.   Default \code{NPOP = 10}.
#' @param crAdaption (numeric (0,1]). Controls the speed of the crossover adaption.  
#'  This parameter is called `c' in the  DEoptim.control help page.  
#'  Default \code{crAdaption = 0}.       
#' @param delta (numeric > 0) A number that controls the convergence. See the DEoptim.control
#'     accuracy of the differential evolution algorithm. Default \code{delta = 1E-8}.
#' @param PRINT (logical) When PRINT = TRUE the algorithm convergence status is printed.
#'     Default  \code{PRINT = FALSE}.
#' @param Seed (integer) Initial random number seed. Default (\code{Seed = NULL}).
#'
#' @return \code{Rbounds} returns the following objects:
#' \itemize{
#'    \item  \strong{R} (matrix) A list of generated correlation matrices.
#'    \item  \strong{converged}: (logical) a logical that indicates the 
#'    convergence status of the optimization for each matrix.
#'    \item \strong{iter} (integer) The number of cycles needed to reach a 
#'    converged solution for each matrix.
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
#' @references Mishra, S. K.  (2007).  Completing correlation matrices
#' of arbitrary order by differential evolution method of global optimization:
#' a Fortran program.  Available at SSRN 968373.
#'
#' @references Mullen, K.M, Ardia, D., Gil, D., Windover, D., Cline, J. (2011). DEoptim: An
#' R Package for Global Optimization by Differential Evolution. \emph{Journal of Statistical Software, 40}, 1-26. URL http://www.jstatsoft.org/v40/i06/.
#'
#' @references Price, K.V., Storn, R.M., Lampinen J.A. (2005) \emph{Differential Evolution - A Practical Approach 
#' to Global Optimization}. Berlin Heidelberg: Springer-Verlag. ISBN 3540209506.
#'
#' @references  Zhang, J. and Sanderson, A. (2009) \emph{Adaptive Differential
#' Evolution}. Springer-Verlag. ISBN 978-3-642-01526-7
#'
#' @author Niels G. Waller
#'
#' @examples
#' ## Example 1: Generate random 4 x 4 Correlation matrices with all rij >= 0.
#'
#'   out <- Rbounds(Nvar = 4,
#'               NMatrices = 4,
#'               Minr = 0,
#'               Maxr = 1,
#'               PRINT = TRUE,
#'               Seed = 1)
#'                       
#'   # Check convergence status of matrices                     
#'   print( table(out$converged) )                     
#'
#'   print( round( out$R[[1]] , 3) )
#'
#'
#' @import DEoptim
#' @export
Rbounds <- function(Nvar = 3,
                    NMatrices = 1,
                    Minr = -1,
                    Maxr = 1,
                    MinEig = 0,
                    MaxIter = 200,
                    NPopFactor = 10,
                    crAdaption = 0,
                    delta = 1E-8,
                    PRINT = FALSE,
                    Seed = NULL){

# initialize a matrix with all NA off-diagonals
  Rna <- matrix(NA, nrow = Nvar, ncol = Nvar)
  diag(Rna) <- 1
  
  p <- nrow(Rna)
  
  #library(parallel)
  Ptype = "none"
  #if(Parallel) Ptype = "parallel"

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
  # Count the number of missing rij values
  num_miss <- p*(p-1)/2


  # After optimization insert final values into R matrix ---
  InsertFinalValues <- function(Rna, X){
      Ri <- Imat
      Ri[upper.tri(Ri, diag=FALSE)] <- X
      Ri <- Ri + t(Ri) - Imat
      Ri
   }#END InsertFinalValues


  # ---- Fit Function ----
  EvoFunc <- function(x){
     Ri <- Imat
     Ri[upper.tri(Ri, diag = FALSE)] <- x
     Ri <- Ri + t(Ri) - Imat
     L <- eigen(Ri, 
                symmetric = TRUE, 
                only.values=TRUE)$values
     
     fit <- 0
     
     # Constrain min eigenvalue >= MinEig
     Lneg <- L[ L < MinEig ]
     
     if( length(Lneg) > 0){
        # penalty to insure PSD
        fit <- 100 *  sum(abs(Lneg))
     }
     
     fit
     
  }#END EvoFunc

 # Run Optimization --------
 # find NMatrices feasible solutions
  
    k = 1
    while( k <= NMatrices){
     
      # Providing (-1,1) start values  speeds up convergence
     initPop <-matrix(runif(NPopFactor*(num_miss^2),  min = Minr, max = Maxr), 
                        nrow=NPopFactor*num_miss, ncol=num_miss)

     suppressWarnings(dout <- tryCatch(
                      DEoptim::DEoptim(fn = EvoFunc,
                               lower = rep(Minr, num_miss),
                               upper = rep(Maxr, num_miss),
                               control = list(trace = FALSE,
                                              VTR = 0,      # Value to be reached
                                              strategy = 2, # 2=default
                                              bs = FALSE,   #TRUE kills program
                                              initialpop = initPop,
                                              NP =  NPopFactor * num_miss,    # NPopFactor * n(n-1)/2,
                                              itermax = MaxIter,
                                              CR = .50, #crossover probability (.5 Default)
                                              F = 0.8, # differential weighting factor, Default = 0.8
                                              storepopfreq = 1,
                                              # higher values of c speed up convergence
                                              # when c is specified then the dist of rij is wrong!
                                              c = crAdaption, #controls speed of crossover Default = 0
                                              # I have yet to find an example in which parallel
                                              # speeds up convergence
                                              parallelType = "none")),
                        error = " "))# do nothing on error
     # assess convergence
       converged = FALSE
       if(round(dout$optim$bestval,8) <= delta) converged <- TRUE



       if(converged){
         if(PRINT) cat(" ", k)
         Rlist[[k]] <- InsertFinalValues(Rna, dout$optim$bestmem)
         ConvergeVec[k] <- converged
         iterVec[k] <- dout$optim$iter
         k  <- k + 1
       } #END if converged

    } #END while loop


  # RETURN LIST
  list(R = Rlist,
       converged = ConvergeVec,
       iter = iterVec)
       #dout = dout)
} #END Rbounds


