#' Complete a Partially Specified Correlation Matrix by the Method of Alternating Projections
#'
#' This function completes a (possibly) partially specified
#' correlation matrix by a modified alternating projections algorithm.
#
#' @param  Rna  (matrix) An n x n incomplete correlation matrix.  Missing entries must
#'     be specified by NA values. If all off diagonal values are NA then the function
#'     will generate a random correlation matrix.
#' @param NMatrices (integer) \code{CompleteRmap} will complete \code{NMatrices}
#' correlation matrices.
#' @param RBounds (logical) If \code{RBounds = TRUE} then the function will attempt to
#' produce a matrix on the surface of the associated elliptope (i.e., the space of all
#' possible PSD R matrices of a given dimension).
#' When \code{RBounds = FALSE}, during each cycle of the alternating projections
#' algorithm all negative eigenvalues of the provisional R matrix are replaced by
#' (sorted) uniform random numbers between the smallest positive eigenvalue and zero (inclusive) of the indefinite matrix.
#' Default \code{RBounds = FALSE}.
#' @param LB (numeric) The lower bound for the random number generator when generating
#' initial estimates for the missing elements of a partially specified correlation matrix.
#' @param UB (numeric) The upper bound for the random number generator when generating
#' initial estimates for the missing elements of a partially specified correlation matrix. Start values
#' (for missing correlations) are sampled from a uniform distribution with bounds \code{[LB, UB]}.
#' @param delta (numeric) A small number that controls the precision of the estimated solution.
#' Default \code{delta = 1E-16}.
#' @param MinLambda (numeric) A small value greater than or equal to 0  used to replace negative
#' eigenvalues during the modified alternating projections algorithm.
#' @param MaxIter (integer) The maximum number of cycles of the
#' alternating projections algorithm. Default \code{MaxIter = 1000}.
#' @param detSort (logical). If \code{detSort = TRUE} then all results will be sorted
#' according to the sizes of the matrix determinants (det(Ri)). Default \code{detSort = FALSE}
#' @param Parallel (logical). If \code{Parallel = TRUE} parallel processing will be used to
#'  generate the completed correlation matrices. Default:  \code{Parallel = FALSE}.
#' @param ProgressBar (logical).  If \code{Parallel = TRUE} and \code{ProgressBar = TRUE} a progress bar
#' will be printed to screen. Default \code{ProgressBar = FALSE}.
#' @param PrintLevel (integer) The \code{PrintLevel} argument can take one of three values:
#'   \itemize{
#'     \item 0  No output will be printed. Default (PrintLevel = 0).
#'     \item 1  Print \code{Delta} and the minimum eigenvalue of the currently completed correlation matrix.
#'     \item 2  Print convergence history.
#'     }
#' @param Digits (integer) Controls the number of printed significant digits if PrintLevel = 2.
#' @param Seed  (integer) Initial random number seed. If reproducible results are desired then
#' it is necessary to specify  \code{ProgressBar = FALSE}. Default \code{Seed = NULL}.
#'@return
#'\itemize{
#'    \item  \strong{CALL} The function call.
#'    \item \strong{NMatrices} The number of completed R matrices.
#'    \item \strong{Rna} The input partially specified R matrix.
#'    \item  \strong{Ri} A list of the completed R matrices.
#'    \item \strong{RiEigs} A list of eigenvalues for each \code{Ri}.
#'    \item \strong{RiDet}  A list of the determinants for each \code{Ri}.
#'    \item \strong{converged} The convergence status (TRUE/FALSE) for each \code{Ri}.
#' }
#'
#' @references
#' Higham, N. J.  (2002).  Computing the nearest correlation matrix: A problem
#' from finance.  IMA Journal of Numerical Analysis, 22(3), 329--343.
#'
#' Waller, N. G.  (2020).  Generating correlation matrices with specified
#' eigenvalues using the method of alternating projections.
#' The American Statistician, 74(1), 21-28.
#'
#' @author Niels G. Waller
#'
#' @examples
#' \dontrun{
#' Rna4 <- matrix(c( 1,  NA,  .29, .18,
#'                   NA, 1,   .11, .24,
#'                  .29, .11, 1,   .06,
#'                  .18, .24, .06, 1), 4, 4)
#'
#' Out4  <- CompleteRmap(Rna = Rna4,
#'                       NMatrices = 5,
#'                       RBounds = FALSE,
#'                       LB = -1,
#'                       UB = 1,
#'                       delta = 1e-16,
#'                       MinLambda = 0,
#'                       MaxIter = 5000,
#'                       detSort = FALSE,
#'                       ProgressBar = TRUE,
#'                       Parallel = TRUE,
#'                       PrintLevel = 1,
#'                       Digits = 3,
#'                       Seed = 1)
#'
#' summary(Out4,
#'         PrintLevel = 2,
#'         Digits = 5)
#' }
#' @import parallel
#' @import pbmcapply
#' @export
#'
CompleteRmap <- function(Rna,
                         NMatrices = 1,
                         RBounds = FALSE,
                         LB = -1,
                         UB = 1,
                         delta = 1E-16,
                         MinLambda = 0,
                         MaxIter = 1000,
                         detSort = FALSE,
                         Parallel = FALSE,
                         ProgressBar = FALSE,
                         PrintLevel = 0,
                         Digits = 3,
                         Seed = NULL){


  CALL <- match.call()

  ## Generate random seed if not supplied
  if( is.null(Seed) ){
    Seed <- as.integer((as.double(Sys.time()) * 1000 + Sys.getpid()) %% 2^31)
  } # END if (is.null(Seed))

  Matrices <- vector(length = NMatrices,
                     mode = "list")

  # initiate vector of convergence values. TRUE = Converged
  Converged <- rep(FALSE, NMatrices)


  # Are you generating a random R matrix or
  # completing a partially specified R matrix
  UpRna <- Rna[upper.tri(Rna, diag=FALSE)]
  GenRandomR <- FALSE
  if( all(is.na(UpRna) ) ) GenRandomR <- TRUE

  # Save the number of iterations at convergence for each matrix
  NIterations <- vector(length = NMatrices)

  p <- nrow(Rna)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  PrjOnSplus <- function(M){
  ## Prj on space of symmetric PSD matrices ----
    VLV <- eigen(M, symmetric = TRUE)
    V <- VLV$vectors
    L <- VLV$values

    if(RBounds == TRUE){
      L[L < 0] <- MinLambda
    } else{
      NNegL <- length(L[L < 0])
      if(NNegL > 0){
        MinPosL <- min(L[L >= 0])
        L[L < 0] <- sort(stats::runif(NNegL, 0, MinPosL),
                         decreasing = TRUE)
      }
    }#End if(RBounds == TRUE)

    M <- V %*% diag(L) %*% t(V)
    # M <- stats::cov2cor(M)
  } ## END  PrjOnSplus

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## Prj on S^n with unit diag & non NA rij ----
  PrjOnU <- function(Ri, Rna){

    Ri[!is.na(Rna)] <- Rna[!is.na(Rna)]
    Ri
  } ## END PrjOnU



 ## Start Values ----
  InsertStartValues <- function(Rna, LB, UB){
    ## Step 1 replace missing r's with uniform rndm numbers ----
    UpRna <- Rna[upper.tri(Rna, diag=FALSE)]
    num_miss <- sum(is.na(UpRna))
    UpRna[is.na(UpRna)]<- stats::runif(num_miss, LB, UB)
    Imat <- diag(nrow(Rna))
    Ri <- Imat
    Ri[upper.tri(Ri, diag=FALSE)] <- UpRna
    Ri <- Ri + t(Ri) - Imat
    Ri
  } ## END InsertStartValues

  ## ----GenerateR----
  GenerateR <- function(Rna, LB, UB){

    # Initiate vars
    k <- 0  #k = Iteration number: Start at zero
    Delta_i = 999
    LastEig = -999


    ## Insert Start Values ----
    Ri <- InsertStartValues(Rna, LB, UB)

    AllEigs <- eigen(Ri, symmetric = TRUE, only.values=TRUE)$values
    LastEig <- AllEigs[p]

    ## If Start values satisfy psd(R) save ----
    if(round(LastEig, 24) >= 0){
      Ri <- list(Ri = Ri, iter = k, Eigs = AllEigs, convergenceStatus = TRUE)
      return(Ri)
    }else {

      ## If Start values do not satisfy psd(R) then ----
      ## ___Project on PSD Splus ----
      Ri <- PrjOnSplus(Ri)

      if(PrintLevel == 2){
        cat("\nIter = ", k,"PrjOnSplus\n")
        print( round(Ri, Digits) )
      }

      ## ___Project on S^n with unit diag ----
      Ri <- PrjOnU(Ri, Rna)

      if(PrintLevel == 2){
        cat("\nIter = ", k,"PrjOnU\n")
        print( round(Ri, Digits) )
      }


      ## Constrain observed r's ----
      while( (Delta_i > delta) || (LastEig <= 1E-32) ){

        # For code checking
        #cat("\nk =", k,"Delta = ", Delta_i,"Lp = ", LastEig)

        #Update iteration number
        k <- k+1

        Ri <- PrjOnSplus(Ri)
        if(PrintLevel == 2){
           cat("\nIter = ",k,"PrjOnSplus\n")
           print(round(Ri, Digits) )
        }

        Ri <- PrjOnU(Ri, Rna)
        
        
        if(PrintLevel == 2){
          cat("\nIter = ",k,"PrjOnU\n")
          print( round(Ri, Digits) )
        }

        ## Compute matrix eigenvalues ----
        AllEigs <- eigen(Ri, symmetric = TRUE, only.values=TRUE)$values
        LastEig <- AllEigs[p]


        ## Delta_i not updated if generating random R
        if( GenRandomR ){
            Delta_i <- 0
        } else{
           Delta_i <- max( (Ri[!is.na(Rna)] - Rna[!is.na(Rna)])^2)
        }


        if(PrintLevel == 1 || PrintLevel == 2){
          cat(c( "\nTotal Iterations = ",k,
                 "\ndelta = ", round(Delta_i,8),
                 "\nMin Eig = ", LastEig))
        }



        if(k > MaxIter){
          warning("\n\n\nTry with new start values")
          # IF algorithm fails to converge return NA, MaxIter
          return(Ri =  list(NA, MaxIter, AllEigs, convergenceStatus = FALSE) )
          break
        }
      }#END while loop
      Ri <-.5 * (Ri + t(Ri))
      # return matrix and number of iterations
      Ri <- list(Ri = Ri, iter = k, Eigs = AllEigs, convergenceStatus = TRUE)
    }#END GenerateR
  }#END else

  PopulateList <- function(i) Rna
  RnaList <- lapply(1:NMatrices,PopulateList)



  if(!Parallel){
    set.seed(Seed)
    Ri<-lapply(RnaList, GenerateR, LB, UB)
  }

  if(Parallel){
    RNGkind("L'Ecuyer-CMRG")
    numcors <- parallel::detectCores(logical = FALSE)

    if(ProgressBar){
      #One cannot get reproducible results with requesting the ProgressBar
       set.seed(Seed, kind="L'Ecuyer-CMRG")
       Ri<-pbmcapply::pbmclapply(RnaList, GenerateR,
                                 LB, UB,mc.cores = numcors-1,
                                 mc.set.seed = TRUE)
    }else{
      set.seed(Seed)
      Ri<-parallel::mclapply(RnaList, GenerateR,
                                LB, UB,mc.cores = numcors-1,
                                mc.set.seed = TRUE)
    }

  }#END if parallel




  ## RETURN ----
  Ri_Matrices <- lapply(Ri, "[[", 1)
  iter <- unlist(lapply(Ri, "[[", 2) )
  Eigs <- lapply(Ri, "[[", 3)
  convergenceStatus <- unlist(lapply(Ri, "[[", 4))
  # Compute matrix determinants
  detVector <- unlist(lapply(Eigs, prod))


  if(detSort == TRUE){
     ##=================================================##
     ## SORT ALL OUTPUT BY SIZE OF MATRIX DETERMINANT ----

     # Find sort order for Matrix determinants
     detSortList <- sort.list(detVector, decreasing = TRUE)

    # Sort Ri_Matrices
    Ri_Matrices_detSorted <- Ri_Matrices[detSortList]
    #Sort Eigs
    Eigs_detSorted <- Eigs[detSortList]
    #Sort number of iterations to convergence
    iter_detSorted <- iter[detSortList]
    Ri_det_detSorted <- detVector[detSortList]
    convergenceStatus_detSorted <- convergenceStatus[detSortList]


    CROut <- list( CALL = CALL,
                   NMatrices = NMatrices,
                   Rna = Rna,
                   Ri = Ri_Matrices_detSorted,
                   iter = iter_detSorted,
                   RiEigs = Eigs_detSorted,
                   RiDet = Ri_det_detSorted,
                   converged = convergenceStatus_detSorted)

    class(CROut) <- "CompleteRmap"
    invisible(CROut)
  }#END if(detSort == TRUE)

  if(detSort == FALSE){

    CROut <- list( CALL = CALL,
                   NMatrices = NMatrices,
                   Rna = Rna,
                   Ri = Ri_Matrices,
                   iter = iter,
                   RiEigs = Eigs,
                   RiDet = detVector,
                   converged = convergenceStatus)

    class(CROut) <- "CompleteRmap"
    invisible(CROut)
  }#END if(detSort == FALSE)


}# END CompleteRmap




