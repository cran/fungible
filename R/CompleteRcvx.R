#' Complete a Partially Specified Correlation Matrix by Convex Optimization
#'
#' This function completes a partially specified
#' correlation matrix by the method of convex optimization.  The completed
#' matrix will maximize the log(det(R)) over the space of PSD R matrices.
#'
#' @param  Rna  (matrix) An n x n incomplete correlation matrix.  Missing entries must
#'     be specified by \code{NA} values.
#' @param Check_Convexity (logical)  If \code{Check_Convexity= FALSE} the
#' program will not check the convexity of the objective function.  Since the
#' convexity of the R completion problem is known to be true, setting this argument
#' to FALSE can decrease computation time.
#' @param PRINT (logical)  If \code{PRINT = TRUE} then the program will print the
#' convergence status of the final solution.
#'
#' @return The \code{CompleteCvxR} function returns the following objects.
#' \itemize{
#'    \item  \strong{R} (matrix) A PSD completed correlation matrix.
#'    \item  \strong{converged}: (Logical) a logical that indicates the convergence status of the optimization.
#'    \item  \strong{max_delta} The maximum absolute difference between the
#'    known elements in the partially specified R matrix and the estimated matrix.
#'    \item \strong{convergence_status} (list) A list containing additional information about the convergence status of the solution.
#' }
#'
#' @references Georgescu, D. I., Higham, N. J., and Peters, G. W.  (2018).  Explicit
#' solutions to correlation matrix completion problems, with
#' an application to risk management and insurance.  Royal Society Open
#' Science, 5(3), 172348.
#'
#' Olvera Astivia, O. L. (2021). A Note on the general solution to completing
#' partially specified correlation matrices. Measurement: Interdisciplinary
#' Research and Perspectives, 19(2), 115--123.
#'
#' @author Niels G. Waller
#'
#'@examples
#'\dontrun{
#'   Rmiss <- matrix(
#'     c( 1,  .25, .6,  .55, .65,  0,  .4,   .6,  .2,  .3,
#'        .25, 1,    0,   0,   0,   0,  NA,   NA,  NA,  NA,
#'        .6,  0,   1,   .75, .75,  0,  NA,   NA,  NA,  NA,
#'        .55, 0,   .75, 1,   .5,   0,  NA,   NA,  NA,  NA,
#'        .65, 0,   .75,  .5, 1,    0,  NA,   NA,  NA,  NA,
#'        0,  0,    0,   0,   0,  1,   NA,   NA,  NA,  NA,
#'        .4, NA,   NA,  NA,  NA,  NA, 1,   .25, .25,  .5,
#'        .6, NA,   NA,  NA,  NA,  NA, .25,  1,  .25,  0,
#'        .2, NA,   NA,  NA,  NA,  NA, .25,  .25, 1,   0,
#'        .3, NA,   NA,  NA,  NA,  NA, .5,    0,   0,  1), 10,10)
#'
#'   out <- CompleteRcvx(Rna = Rmiss,
#'                       Check_Convexity = FALSE,
#'                       PRINT = FALSE)
#'
#'   round(out$R, 3)
#'}
#' @importFrom CVXR Variable Maximize Problem solve
#' @export
#'
CompleteRcvx <- function(Rna, 
                        Check_Convexity = TRUE, 
						PRINT = TRUE){

    Nvar <- nrow(Rna)

    R <- CVXR::Variable(Nvar, Nvar, PSD=TRUE)
    constraints <- list(R[!is.na(Rna)] == Rna[!is.na(Rna)])

    objective <-  CVXR::Maximize(CVXR::log_det(R))

    program <- CVXR::Problem(objective, constraints)

    soln <- CVXR::solve(program,
                  ignore_dcp = Check_Convexity)

    Rcomplete <- as.matrix(soln$getValue(R))

    #Assess solution accuracy
    max_delta <- max(abs(Rcomplete[!is.na(Rna)] - Rna[!is.na(Rna)]))

    # clean up matrix
    diag(Rcomplete) <- 1
    Rcomplete[!is.na(Rna)] <- Rna[!is.na(Rna)]

    converged <- FALSE

    # The status of the solution.
    # Can be "optimal", "optimal_inaccurate",
    # "infeasible", "infeasible_inaccurate",
    # "unbounded", "unbounded_inaccurate", or
    # "solver_error".
    if(soln$status == "optimal") converged <- TRUE

    if(PRINT){
        cat("\nConvergence status = ", converged,
            "\nmax_delta = ", max_delta, "\n\n")
    }

    list (R = Rcomplete,
          converged = converged,
          max_delta = max_delta,
          convergence_status = soln$status )
}#END CompleteCvxR






