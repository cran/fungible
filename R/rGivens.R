

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


