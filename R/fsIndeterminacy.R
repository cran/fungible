#' Understanding Factor Score Indeterminacy with Finite Dimensional Vector Spaces
#' 
#' 
#' This function illustrates the algebra of factor score indeterminacy 
#' using concepts from finite dimensional vector spaces. Given any factor loading 
#' matrix, factor correlation matrix, and desired sample size, the program will 
#' compute a matrix of observed scores and multiple sets of factors
#' scores.  Each set of (m common and p unique) factors scores will fit the model 
#' perfectly.
#' @param Lambda (Matrix) A p x m matrix of factor loadings.
#' @param Phi (Matrix) An m x m factor correlation matrix.
#' @param N (Integer) The desired sample size. 
#' @param X (Matrix) an optional N x p matrix of observed scores. Note that the observed scores
#' are expected to fit the factor model (as they will if they are generated 
#' from simFA and Population = TRUE is specified). Default (\code{X = NULL}).
#' @param SeedX (Integer) Starting seed for generating the matrix of observed scores, X.
#' @param SeedBasis (Integer) Starting seed for generating a basis for all scores. 
#' @param SeedW  (Integer) Starting seed for generating a weight matrix that is 
#' used to construct those parts of the factor scores that lie outside of span(X).
#' @param SeedT (Integer) Starting seed for generating a rotation matrix that 
#' creates a new set of factor scores from an existing set of scores such that 
#' the new set also perfectly fits the factor model. 
#' @param DoFCorrection (Logical) Degrees of freedom correction.  If DoFCorrection = TRUE
#' then var(x) = 1/(N-1) * t(x) \%*\% x; else var(x) = 1/N * t(x) \%*\% x. 
#' Default (\code{DoFCorrection = TRUE}). 
#' @param Print (Character) If \code{Print = "none"} no summary information
#' will be printed.  If \code{Print = "short"} then basic output for evaluating
#' the factor scores will be printed. If \code{Print = "long"} extended output
#' will be printed. Default (\code{Print = "short"}).
#' @param Digits  (Integer) Sets the number of significant digits to print when
#' printing is requested. 
#' @param Example (Logical) If Example = TRUE the program will
#' execute the orthogonal two factor model  described  in Waller (2021). 
#' @return 
#' \itemize{
#'  \item \strong{"Sigma"}: The p x p model implied covariance matrix.
#'  \item \strong{"X"}:  An N x p data matrix for the observed variables. 
#'  \item \strong{"Fhat"}: An N x (m + p) matrix of regression factor score estimates. 
#'  \item \strong{"Fi"}:  A possible set of common and unique factor scores.
#'  \item \strong{"Fj"}: The set of factor scores that are minimally correlated with Fi. 
#'  \item \strong{"Fk"}: Another set of common and unique factor scores. 
#'  Note that in a 1-factor model, Fk = Fi.
#'  \item \strong{"Fl"}: The set of factor scores that are minimally correlated with Fk. 
#'  Note that in a 1-factor model, Fj = Fl. 
#'  \item \strong{"Ei"}: Residual scores for Fi.
#'  \item \strong{"Ej"}: Residual scores for Fj.
#'  \item \strong{"Ek"}: Residual scores for Fk.
#'  \item \strong{"El"}: Residual scores for Fl. 
#'  \item \strong{"L"}: The factor loading super matrix.
#'  \item \strong{"C"}: The factor correlation super matrix.
#'  \item \strong{"V"}: A (non unique) basis for R^N.
#'  \item \strong{"W"}: Weight matrix for generating  Zi.
#'  \item \strong{"Tmat"}: The orthogonal transformation matrix used to construct Fk  from Fi .
#'  \item \strong{"B}:  The matrix that takes Ei to Ek (Ek = Ei B).
#'  \item \strong{"Bstar"} In an orthogonal factor model, Bstar takes Fi to Fk (Fk = Fi Bstar). 
#'               In an oblique model the program returns Bstar=NULL.
#'  \item \strong{"P"}: The matrix that imposes the proper covariance structure on Ei.
#'  \item \strong{"SeedX"}: Starting seed for X.   
#'  \item \strong{"SeedBasis"}: Starting seed for the basis. 
#'  \item \strong{"SeedW"}: Starting seed for weight matrix W. 
#'  \item \strong{"SeedT"}: Starting seed for rotation matrix T.
#'  \item \strong{"Guttman"}:  Guttman indeterminacy measures for the common and unique factors.
#'  \item \strong{"CovFhat"}: Covariance matrix of estimated factor scores.
#'  }
#' @references Guttman, L.  (1955).  The determinacy of factor score matrices 
#'  with applications for five other problems of common factor theory.  
#'  \emph{British Journal of Statistical Psychology, 8}, 65-82. 
#' @references Ledermann, W.  (1938).  The orthogonal transformation of a factorial matrix 
#'  into itself.  \emph{Psychometrika, 3}, 181-187. 
#' @references Schönemann, P. H. (1971). The minimum average correlation 
#' between equivalent sets of uncorrelated factors. \emph{Psychometrika, 36}, 
#' 21-30. 
#' @references Steiger, J. H. and Schonemann, P. H.  (1978).  In Shye, S. (Ed.), 
#'  \emph{A history of factor indeterminacy} (pp. 136--178). San  Francisco: Jossey-Bass. 
#' @references Waller, N. G. (2021) Understanding factor indeterminacy through the lens of finite 
#' dimensional vector spaces. Manuscript under review.
#'
#' @family Factor Analysis Routines
#'
#' @author Niels G. Waller (nwaller@umn.edu)
#' @examples
#' # ---- Example 1: ----
#' # To run the example in Waller (2021) enter:
#' out1 <- fsIndeterminacy(Example = TRUE)
#' 
#' # ---- Example 1: Extended Version: ----
#' 
#' N <- 10 # number of observations
#' # Generate Lambda: common factor loadings 
#' #          Phi: Common factor correlation matrix
#' 
#' Lambda <- matrix(c(.8,  0,
#'                    .7,  0,
#'                    .6,  0,
#'                     0, .5,
#'                     0, .4,
#'                     0, .3), 6, 2, byrow=TRUE)
#' 
#' out1  <- fsIndeterminacy(Lambda,
#'                          Phi = NULL,    # orthogonal model
#'                          SeedX = 1,     # Seed for X
#'                          SeedBasis = 2, # Seed for Basis
#'                          SeedW = 3,     # Seed for Weight matrix
#'                          SeedT = 5,     # Seed for Transformation matrix
#'                          N = 10,        # Number of subjects
#'                          Print = "long",
#'                          Digits = 3)
#'
#' # Four sets of factor scores
#'   Fi <- out1$Fi
#'   Fj <- out1$Fj
#'   Fk <- out1$Fk
#'   Fl <- out1$Fl
#'
#' # Estimated Factor scores
#'   Fhat <- out1$Fhat
#'
#' # B wipes out Fhat (in an orthogonal model)
#'   B <- out1$B
#'   round( cbind(Fhat[1:5,1:2], (Fhat %*% B)[1:5,1:2]), 3) 
#' 
#' # B takes Ei -> Ek
#'   Ei <- out1$Ei
#'   Ek <- out1$Ek
#'   Ek - (Ei %*% B)
#'
#' # The Transformation Approach
#' # Bstar takes Fi --> Fk
#'   Bstar <- out1$Bstar
#'   round( Fk - Fi %*% Bstar, 3)
#'
#' # Bstar L' = L'
#'   L <- out1$L
#'   round( L %*% t(Bstar), 3)[,1:2]  
#'
#'
#' # ---- Example 3 ----
#' # We choose a different seed for T
#'
#' out2  <- fsIndeterminacy(Lambda , 
#'                         Phi = NULL, 
#'                         X = NULL,
#'                         SeedX = 1,     # Seed for X 
#'                         SeedBasis = 2, #  Seed for Basis
#'                         SeedW = 3,     #  Seed for Weight matrix
#'                         SeedT = 4,     # Seed for Transformation matrix
#'                         N,             
#'                         Print = "long",
#'                         Digits = 3,
#'                         Example = FALSE)
#'  Fi <- out2$Fi
#'  Fj <- out2$Fj
#'  Fk <- out2$Fk
#'  Fl <- out2$Fl
#'  X  <- out2$X
#'  
#' # Notice that all sets of factor scores are model consistent 
#'  round( t( solve(t(Fi) %*% Fi) %*% t(Fi) %*% X) ,3)
#'  round( t( solve(t(Fj) %*% Fj) %*% t(Fj) %*% X) ,3)
#'  round( t( solve(t(Fk) %*% Fk) %*% t(Fk) %*% X) ,3)
#'  round( t( solve(t(Fl) %*% Fl) %*% t(Fl) %*% X) ,3)
#'  
#' # Guttman's Indeterminacy Index
#' round( (1/N * t(Fi) %*% Fj)[1:2,1:2], 3)
#' 
#' @export
fsIndeterminacy <- function(Lambda = NULL, 
                            Phi = NULL, 
                            N = NULL,
                            X = NULL,
                            SeedX = NULL,
                            SeedBasis = NULL,
                            SeedW = NULL,
                            SeedT = 1,
                            DoFCorrection = TRUE,
                            Print = "short",
                            Digits = 3,
                            Example = FALSE){
  
 
 # Lambda: The p x m factor pattern matrix
 # Phi: The  m x m factor correlation matrix
 # N: Sample size  
 # SeedX: Starting seed for the Observed Scores (X)  
 # SeedBasis: Starting seed for the Basis
 # SeedW: Starting seed for the weight matrix: W  
 # SeedT: Starting seed for rotation matrix: T  
 # Print:  (character)  Print = {"none", "short", "long"}
 # Digits: Number of significant digits to print in model output.
 # Example (logical) if Example = TRUE print example from paper
  
  
  # if X is supplied we assume that it has been standardized
  # with sample sdevs s.t. DF = N - 1
  XSupplied = FALSE
  if(!is.null(X)) XSupplied = TRUE
  
  if(!is.null(SeedX) && !is.null(SeedBasis)){
    if(SeedX == SeedBasis) stop("\nFATAL ERROR: SeedX and SeedBasis must be differet numbers\n")
  }
   
  ## generate random seeds if not supplied
  if(is.null(SeedX)) SeedX<- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) 
  if(is.null(SeedBasis)) SeedBasis<- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) 
  if(is.null(SeedW)) SeedW<- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) 
  
  
  if( is.null(Lambda) & isFALSE(Example) ) stop("\n\nFATAL ERROR: Lambda is a required argument")
  if(is.null(N) & isFALSE(Example)) stop("\n\nN must be specified\n") 
  
  
  # --- Example == TRUE ----
  if(Example){
     SeedX = 1
     SeedBasis = 2
     SeedW = 3
     Lambda <- matrix(c(.8,  0,
                        .7,  0,
                        .6,  0,
                         0, .5,
                         0, .4,
                         0, .3), 6, 2, byrow=TRUE)
      Phi <- matrix(c(1,.0,.0,1),2, 2)
      N = 10
  }# END if(Example)
  
   
   m <- ncol(Lambda)   # Number of common factors
   p <- nrow(Lambda)
   # Create Phi if not supplied
   if(is.null(Phi)) Phi <- diag(m)
   
   # Check if N is large enough
   if(N < (m + p + 1)) stop("\n\n N must be > ", m + p + 1, "\n")
   
   # create the factor correlation super matrix
   C <- diag(p + m)
   C[1:m, 1:m] <- Phi
   # calculate commonalities
   h2 <- diag(Lambda %*% Phi %*% t(Lambda))
   # calculate Psi
   Psi <- diag(sqrt(1 - h2))
   # concatenate Lambda and Psi to create L
   L <- cbind(Lambda, Psi)
   
   # calculate  (model implied) R.X  and its inverse
   # and generate X 
   R.X <- L %*% C %*% t(L)
   R.X.inv <- solve(R.X)
   
   # This is the easiest method for generating X
   # X <- MASS::mvrnorm(n = N, 
   #                       mu = rep(0, p), 
   #                       Sigma = R.X, 
   #                       empirical = TRUE)
   # Rescale X so 1/N X'X is a correlation matrix
   # X <- sqrt(N)/sqrt(N-1) * X
   
   # This is the method described in the paper
   I_N <- diag(N)
   J <- matrix(1, nrow = N, ncol = 1)
   P.J <- 1/N * J %*% t(J)
  
   if(is.null(X)){
      set.seed(SeedX)
      S1 <- matrix(rnorm(N * p), N, p)
      # S2 = matrix of random deviation scores
      S2 = (I_N - P.J) %*% S1
      K1 <- chol(1/N * t(S2) %*% S2)
      K2 <- chol(R.X)
      X <-  S2 %*% solve(K1) %*% K2
   } # END if(is.null(X))
   
   # # &&&&&&&&
   # if(XSupplied == FALSE & DoFCorrection == TRUE){
   #    X = sqrt((N - 1)/(N)) * X
   # } #END dof correction    
   # 

   # Calculate Fhat: the LS factor score estimates 
    Fhat <- X %*% R.X.inv %*% L %*% C  
   # Fhat is in the span (column space) of X

   # Calculate X^+ and Projection matrix P.X 
    X.ginv <- solve( t(X) %*% X) %*% t(X)
    P.X <- X %*% X.ginv

   # Create V: a basis for R^N ( N dimensional space)
   #   1. the leading p cols of V = X
   #   2. the next  N - p - 1  cols are random vectors in 
   #      deviation score form that are orthogonal to X
   #   3. the final column = J, vector of all ones 
    set.seed(SeedBasis)  
    Zstar <- matrix(rnorm(N * (N - p - 1)), nrow = N, ncol = (N-p-1))
   # project Zstar into the left null space J and then
   # the left null space of X
    Z <- (I_N - P.X  )  %*% (I_N - P.J) %*% Zstar
   # standardize cols of  Z  (not really necessary)
    Z <- sqrt(N-1)/sqrt(N) * apply(Z, 2, scale)
   # Create V
    V <- cbind(X, Z, J)  
    
    
    # ----Label Matrices ----
    #Add column names to basis: V
    colnames(V) <- c(paste0("x",1:p), paste0("z",1:(N-p-1)),"J")
    row.names(V) <- rep("", N)
    colnames(Lambda) <- paste0("f",1:m)
    row.names(Lambda) <- paste0("x", 1:p)
    colnames(Phi) <- rownames(Phi) <- paste0("f",1:m)
    colnames(X) <- paste0("x", 1:p)
    rownames(X) <- paste0("Subj", 1:N)
    
    
    
   # Generate F 
   # Take any collection of m randomly weighted
   # vectors in Z
   set.seed(SeedW)  
   W <- matrix(runif( (N-p-1)*m, -1,1), nrow = N-p-1, ncol = m)
   
   if(Example){
      W <- matrix(c( 1, 0, 
                     0, 1, 
                     0, 0), 3, 2, byrow=TRUE)
   }
   
   
   Zr <- Z %*% W
   
   
   # X created by program
   # orthogonalize and standardize Zr. Zi is an orthogonal basis
   # for a random m-dimensional subspace of Z
   Zi <- sqrt(N) * svd(Zr)$u
   
   # X supplied by user
   # if standardized X is supplied then we assume that
   # t(x) %*% x = N - 1. This is true for simFA
   if( XSupplied ) Zi <- sqrt(N - 1) * svd(Zr)$u
   
   # Using the Zi basis, create p+m vectors with
   # desired (NPD) Cov matrix: C - CovFhat
   
   CovFhat <- C %*% t(L) %*% R.X.inv %*% L %*% C
  
   # Define K = C - CovFhat 
   # note that K is NPD (with rank m)
   K <- C - CovFhat
   # compute eigenstructure of K
   QDQt <- eigen(K, symmetric = TRUE)
   Q <- QDQt$vectors
   D <- QDQt$values 
   
   # ---- Generate T ----
   # generate a random orthogonal matrix from a (uniform) 
   # Haar distribution
   # T can be interpreted as an orthogonal rotation matrix
   set.seed(SeedT)     
   M <- matrix(rnorm(m * m), nrow = m, ncol = m)
   Tmat <- qr.Q(qr(M))
   
    
   if(Example){
      # rotate basis by 45 deg
      Tmat <- matrix(c( cos(pi/4), -sin(pi/4), 
                     sin(pi/4),  cos(pi/4)), 2, 2, byrow=TRUE)
   } 
   
   # By including T in the construction equation we combine 
   # the (a) construction and (b) transformation approach
   # to generating factor scores.
   
  
   # Calculate Fi and Fj [without T]
   # Ei = that part of F that is in the column space of Z
   
   if(m == 1){  # Nfac = 1
     P <-  as.matrix(Q[,1:m, drop=FALSE] * sqrt(D[1]) )
   }
   if(m >= 2){  # NFac > 1 
   P <-  as.matrix(Q[,1:m, drop=FALSE] %*% diag(sqrt(D[1:m,drop=FALSE]))) 
   }
  
     
   # ---- Construct B ----
   B <- P %*% solve( t(P) %*% P) %*% Tmat %*% t(P)

   # ---- Construct Bstar ----
   V_Fhat <- svd(Fhat)$v[,1:p]
   Bstar <- (V_Fhat %*% t(V_Fhat)  + B)
   
  
   # Bstar only works for orthogonal factor models
   if( round( sum(C - diag(p+m)), 8) > 0 ){
      Bstar = NULL
   }
   
   
   Ei <-  Zi %*% t(P) 
   Ej <- -Ei # to compute the least correlated set of FS with Fk
     
   # Calculate Fk and Fl [with T]
   # we could do this:
   # Ek <-  Zi %*% Tmat %*% t(P)
   # This is an alternative way of executing the 
   # transformation method of Thomson (1935), Ledermann (1938), 
   # and Schonemann (1971)
   Ek <- Ei %*% B
   El <- -Ek # to compute the least correlated set of FS with Fk
   
   
   #---- Degrees of Freedom Correction ----
   if(XSupplied == FALSE & 
      DoFCorrection == TRUE){
     
     s <- sqrt((N - 1)/( (N )) )
     X <- s * X
     Fhat <- s * Fhat
     Ei   <- s * Ei
     Ej   <- s * Ej
     Ek   <- s * Ek
     El   <- s * El
     
   } #END DofCorrection

   
   
   Fi <- Fhat + Ei
   Fj <- Fhat + Ej
   Fk <- Fhat + Ek
   Fl <- Fhat + El
   

   # Gutmman's indeterminacy measure
   rho <- sqrt(diag(CovFhat))
   # Correlation of Fk and Fj where Fk = Fhat + Ek
   Guttman <- 2*rho^2 - 1
   

      #----fnc to Print Section Headings 
      sectionHeading <- function(string){
         stringLength <- nchar(string)
         cat("\n\n")
         cat(c(paste(rep("=", stringLength ),collapse =""),"\n"))
         cat(c(string, "\n"))
         cat(c(paste(rep("=", stringLength ),collapse =""),"\n"))
      } # END sectionHeading
      
      
      # ---- Print short ----
      if(Print == "short"){
         sectionHeading("A Demonstration of Factor Score Indeterminacy")
         
         sectionHeading("Lambda: Factor Pattern Matrix")
         print(round(Lambda, Digits) )
         
         sectionHeading("Phi: Factor Correlation Matrix")
         print(round(Phi, Digits) )
         
         sectionHeading("Observed Scores")
         print(round(X,Digits))
      

         sectionHeading("Two Alternative Sets of Common Factor Scores")
         Fcommon <- cbind(Fi[,1:m], Fj[,1:m])
         rownames(Fcommon) <- paste0("Subj",1:N)
         colnames(Fcommon) <- c(paste0("f", 1:m, "i"), paste0("f", 1:m, "j"))
         print(round( Fcommon, Digits ) )
        
         FcommonT <- cbind(Fk[,1:m], Fl[,1:m])
         rownames(FcommonT) <- paste0("Subj",1:N)
         colnames(FcommonT) <- c(paste0("f", 1:m, "k"), paste0("f", 1:m, "l"))
         sectionHeading("Transformed Sets of Common Factor Scores")
         print(round( FcommonT, Digits ) )
         
         sectionHeading("Corr(Fhat, F) for common factors")
         rvec <- matrix(sqrt(diag(CovFhat[1:m, 1:m, drop = FALSE])),1,m)
         colnames(rvec) <- paste0("f",1:m)
         row.names(rvec) <- ""
         print(round(rvec , Digits) )
         
         sectionHeading("Guttman's Indeterminacy Index: Corr(Fi,Fj)")
         Gttmn <- Guttman[1:m]
         names(Gttmn) <- paste0("f",1:m)
         print(round(Gttmn, Digits))
      }# END print short version
      
      
      # ---- Print long ----
      if(Print == "long"){
      sectionHeading("A Demonstration of Factor Score Indeterminacy")
     
      sectionHeading("Lambda: Factor Pattern Matrix")
      print(round(Lambda, Digits) )
      
      sectionHeading("Phi: Factor Correlation Matrix")
      print(round(Phi, Digits) )
      
      sectionHeading("Observed Scores")
      print(round(X,Digits))
      
      sectionHeading("A basis for R^N")
      print(round(V,Digits))
      
      sectionHeading("W: Weight matrix used to construct a basis for Ei")
      print(round( W, Digits ))
      
      sectionHeading("Zi: An Orthogonal Basis for Ei")
      print(round( Zi, Digits ))
      
      sectionHeading("Two Alternative Sets of Common Factor Scores")
      Fcommon <- cbind(Fi[,1:m], Fj[,1:m])
      rownames(Fcommon) <- paste0("Subj",1:N)
      colnames(Fcommon) <- c(paste0("f", 1:m, "i"), paste0("f", 1:m, "j"))
      print(round( Fcommon, Digits ) )
      
      sectionHeading("T: Basis Rotation Matrix")
      print(round( Tmat, Digits ))
      
      sectionHeading("B: The matrix that takes Ei to Ek (Ek = Ei B)")
      print(round( B, Digits ))
         
      sectionHeading("Transformed Sets of Common Factor Scores")
      FcommonT <- cbind(Fk[,1:m], Fl[,1:m])
      
      rownames(FcommonT) <- paste0("Subj",1:N)
      colnames(FcommonT) <- c(paste0("f", 1:m, "k"), paste0("f", 1:m,"l"))
      print(round( FcommonT, Digits ) )
      
      sectionHeading("Estimated Common Factor Scores")
      Fhatcommon <- Fhat[,1:m, drop=FALSE]
      rownames(Fhatcommon) <- paste0("Subj",1:N)
      colnames(Fhatcommon) <- paste0("f", 1:m, "hat")
      print(round( Fhatcommon, Digits ) )

      
      sectionHeading("Corr(Fhat, F) for common factors")
      rvec <- matrix(sqrt(diag(CovFhat[1:m, 1:m, drop = FALSE])),1,m)
      colnames(rvec) <- paste0("f",1:m)
      row.names(rvec) <- ""
      print(round(rvec , Digits) )
         
      sectionHeading("Guttman's Indeterminacy Index: Corr(Fi,Fj)")
      Gttmn <- Guttman[1:m, drop = FALSE]
      names(Gttmn) <- paste0("f",1:m)
      print(round(Gttmn, Digits))
   }# END print long version
   
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   #   ---- Return ----
   list(Sigma = R.X,  # model implied cov matrix
        X = X,        # Data for observed variables
        Fhat = Fhat,  # L.S. estimated common and unique fac scores
        Fi = Fi,      # A possible set of factor scores
        Fj = Fj,      # Another possible set of factor scores
        Fk = Fk,
        Fl = Fl,
        Ei = Ei,      # Residual scores for Fi
        Ej = Ej,      # Residual scores for Fj
        Ek = Ek,
        El = El,
        L = L,        # factor loading super matrix
        C = C,        # factor correlation super matrix
        V = V,        # a basis for R^N
        W = W,        # weights to generate Zr
        Tmat = Tmat,  # orthogonal transformation matrix
        B = B,        # matrix that takes Ei to Ek
        Bstar = Bstar,# matrix that takes Fi to Fk
        P = P,        # matrix that generates proper cov for Ei
        SeedX = SeedX,   
        SeedBasis = SeedBasis,
        SeedW = SeedW,
        SeedT = SeedT,
        Guttman = Guttman, # Guttman measure of indeterminacy
        CovFhat = CovFhat) # Cov matrix of estimated factor scores
   
 } ##END fsIndeterminacy


