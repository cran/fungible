#' Generate Correlation Matrices with Specified Eigenvalues
#' 
#' rGivens generates correlation matrices with user-specified eigenvalues via a
#' series of Givens rotations by methods described in Bendel & Mickey (1978)
#' and Davis & Higham (2000).
#' 
#' 
#' @param eigs A vector of eigenvalues that must sum to the order of the
#' desired correlation matrix. A fatal error will occur if sum(eigs) !=
#' length(eigs).
#' @param Seed Either a user supplied seed for the random number generator or
#' `NULL' for a function generated seed. Default Seed = `NULL'.
#' @return \item{R}{A correlation matrix with desired spectrum.}
#' \item{Frob}{The Frobenius norm of the difference between the initial and
#' final matrices with the desired spectrum.} \item{convergence}{(Logical) TRUE
#' if rGivens converged to a feasible solution, otherwise FALSE.}
#' @references Bendel, R. B. & Mickey, M. R. (1978). Population correlation
#' matrices for sampling experiments, Commun. Statist. Simulation Comput., B7,
#' pp. 163-182.
#' 
#' Davies, P. I, & Higham,N. J. (2000). Numerically stable generation of
#' correlation matrices and their factors, BIT, 40 (2000), pp. 640-651.
#' @keywords datagen
#' @export
#' @examples
#' 
#' 
#' ## Example
#' ## Generate a correlation matrix with user-specified eigenvalues
#' 
#' out <- rGivens(c(2.5, 1, 1, .3, .2), Seed = 123)
#' 
#' #> eigen(out$R)$values
#' #[1] 2.5 1.0 1.0 0.3 0.2
#' 
#' print(out)
#' #$R
#' #           [,1]       [,2]        [,3]        [,4]       [,5]
#' #[1,]  1.0000000 -0.1104098 -0.24512327  0.46497370  0.2392817
#' #[2,] -0.1104098  1.0000000  0.33564370 -0.46640155 -0.7645915
#' #[3,] -0.2451233  0.3356437  1.00000000 -0.02935466 -0.2024926
#' #[4,]  0.4649737 -0.4664016 -0.02935466  1.00000000  0.6225880
#' #[5,]  0.2392817 -0.7645915 -0.20249261  0.62258797  1.0000000
#' #
#' #$Frob
#' #[1] 2.691613
#' #
#' ##$S0
#' #           [,1]        [,2]        [,3]        [,4]        [,5]
#' #[1,]  1.0349665  0.22537748 -0.46827121 -0.10448336 -0.24730565
#' #[2,]  0.2253775  0.31833805 -0.23208078  0.06591368 -0.14504161
#' #[3,] -0.4682712 -0.23208078  2.28911499  0.05430754  0.06964858
#' #[4,] -0.1044834  0.06591368  0.05430754  0.94884439 -0.14439623
#' #[5,] -0.2473056 -0.14504161  0.06964858 -0.14439623  0.40873606
#' #
#' #$convergence
#' #[1] TRUE
#' 
#' 
rGivens <- function(eigs, Seed = NULL){
# rGivens: Generate R matrices with user-specified eigenvalues via
# Givens rotations based on theory described in:
  
#  References:
#   [1] R. B. Bendel and M. R. Mickey, Population correlation matrices
#      for sampling experiments, Commun. Statist. Simulation Comput.,
#      B7 (1978), pp. 163-182.
#   [2] P. I. Davies and N. J. Higham, Numerically stable generation of
#       correlation matrices and their factors, BIT, 40 (2000), pp. 640-651.
#


  eps<-1e-12
  n <- length(eigs)
  
  

   if(abs(sum(eigs)-n)/n > 100*eps){
     stop("Elements of eigs must sum to n,\n")
   }


  A = diag(eigs);
## generate a random orthogonal matrix
  if (is.null(Seed)) 
    Seed <- as.integer((as.double(Sys.time()) * 1000 + Sys.getpid())%%2^31)
  set.seed(Seed)
  
  M <- matrix(rnorm(n * n), nrow = n, ncol = n)
  Q <- qr.Q(qr(M))
  S0 <- A <- Q %*% A %*% t(Q);  # Not exploiting symmetry here.

  
  a = diag(A);
  y = which(a < 1);
  z = which(a > 1);
  

while( length(y) > 0 && length(z) > 0){
  
    i = y[ceiling(runif(1)*length(y))];
    j = z[ceiling(runif(1)*length(z))];
  if( i > j){
     temp = i; 
     i = j; 
     j = temp
  }
  
  alpha = sqrt(A[i,j]^2 - (a[i]-1)*(a[j]-1));
  
  t<-c(0,0)
  t[1] = (A[i,j] + sign(A[i,j])*alpha)/(a[j]-1);
  
  t[2] = (a[i]-1)/((a[j]-1)*t[1]);
  
  t = t[ceiling(runif(1)*2)];  # Choose randomly from the two roots.
  c = 1/sqrt(1 + t^2);
  s = t*c;
  
  A[, c(i,j)] =  A[, c(i,j)] %*% matrix(c(c, s, -s, c),2,2, byrow=TRUE);
  A[c(i,j),] = matrix(c(c, -s, s,c), 2, 2, byrow=TRUE) %*% A[c(i,j), ];
  
  # Ensure (i,i) element is exactly 1.
  A[i,i] = 1;
  
  a = diag(A);
  y = which(a < 1);
  z = which(a > 1);
  
} ## End while loop
  
  # Setlast diagonal to 1:
  diag(A) = 1
  
  A = (A + t(A))/2; # Force symmetry.
  convergence <- TRUE
  if(max(abs(A))>1) convergence <- FALSE
  
  normF<-function(M) sqrt(sum(M^2))
  Frob = normF(A - S0)  
 
  list(R = A, Frob = Frob, S0 = S0, convergence = convergence)
} # END rGivens
##################################################


