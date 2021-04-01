# August 15, 2020
#' Investigate local minima in  faMain objects
#' 
#' 
#' Compute pairwise RMSD values among rotated factor patterns from 
#' an \code{faMain} object.
#' 
#' @description Compute pairwise root mean squared deviations (RMSD) 
#' among rotated factor patterns in an \code{faMain} object. 
#' Prior to computing the RMSD values, each pair of solutions is aligned to 
#' the first member of the pair.  Alignment is accomplished using the 
#' Hungarian algorithm as described in \code{faAlign}.
#' 
#' @param fout (Object from class  \code{faMain}). 
#' @param Set  (Integer) The index of the solution set (i.e., the collection of 
#'   rotated factor patterns with a common complexity value) from an 
#'   \code{faMain} object. 
#' @param HPthreshold (Scalar) A number between [0, 1] that defines the 
#'   hyperplane threshold. Factor pattern elements below \code{HPthreshold} in absolute 
#'   value are counted in the hyperplane count.
#' @param digits (Integer) Specifies the  number of significant 
#'   digits in the printed output. Default \code{digits = 5}.
#' @param PrintLevel (Integer) Determines the level of printed output.
#'  PrintLevel = 
#'  \itemize{
#'   \item \strong{0}: No output is printed. 
#'   \item \strong{1}: Print output for the six most discrepant pairs of 
#'   rotated factor patterns.
#'   \item \strong{2}: Print output for all  pairs of rotated factor patterns.
#'  }
#' @return \code{faLocalMin} function will produce the following output.
#' \itemize{
#'   \item \strong{rmsdTable}: (Matrix) A table of \code{RMSD} values for each  pair of 
#'       rotated factor patterns in  solution set \code{Set}.
#'   \item \strong{Set}: (Integer) The index of the user-specified solution set.
#'   \item \strong{complexity.val} (Numeric): The common complexity value for all members 
#'   in the user-specified solution set.
#'   \item \strong{HPcount}: (Integer) The hyperplane count for each factor pattern in the solution set.
#'   }
#' @author Niels Waller
#' @keywords fungible Statistics
#' @family Factor Analysis Routines
#' @import clue
#' @export
#' @examples
#' \dontrun{
#'   ## Generate Population Model and Monte Carlo Samples ####
#'   sout <- simFA(Model = list(NFac = 5,
#'                           NItemPerFac = 5,
#'                            Model = "orthogonal"),
#'               Loadings = list(FacLoadDist = "fixed",
#'                               FacLoadRange = .8),
#'               MonteCarlo = list(NSamples = 100, 
#'                                 SampleSize = 500),
#'               Seed = 655342)
#' 
#'   ## Population EFA loadings
#'   (True_A <- sout$loadings)
#' 
#'   ## Population Phi matrix
#'   sout$Phi
#' 
#'   ## Compute EFA on Sample 67 ####
#'   fout <- faMain (R = sout$Monte$MCData[[67]],
#'                 numFactors = 5,
#'                 targetMatrix = sout$loadings,
#'                 facMethod = "fals",
#'                 rotate= "cfT",
#'                 rotateControl = list(numberStarts = 50,
#'                                      standardize="CM",
#'                                      kappa = 1/25),
#'                 Seed=3366805)
#' 
#'   ## Summarize output from faMain
#'   summary(fout, Set = 1, DiagnosticsLevel = 2, digits=4)
#' 
#'   ## Investigate Local Solutions
#'   LMout <- faLocalMin(fout, 
#'                     Set = 1,
#'                     HPthreshold = .15,
#'                     digits= 5, 
#'                     PrintLevel = 1)
#'                     
#'   ## Print hyperplane count for each factor pattern 
#'   ## in the solution set
#'   LMout$HPcount
#'   }

faLocalMin <- function(fout, 
                       Set = 1, 
                       HPthreshold = .10,
                       digits = 5,
                       PrintLevel = 1){
 
  SolSet = Set
  

  faAlign2 <- function(F1, F2, Phi2 = NULL, MatchMethod = "LS"){
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
    ## date: Nov 4, 2017
    ## Author: Niels Waller
    ##
    ## faAlign will align the factors of F2 to F1
    ## by minimizing the least squares (LS) fit of (F1 - F2)
    ## or by maximizing the matched-factor congruence coefficients (CC).
    ##
    ## If a unique match is not found in the initial attempt then the program 
    ## uses the Hungarian algorithm to find an optimal match.
    ##
    ## see:
    ## Harold W. Kuhn, "The Hungarian Method for the assignment problem", Naval Research 
    ##   Logistics Quarterly, 2: 83-97, 1955. Kuhn's original publication.
    ## Harold W. Kuhn, "Variants of the Hungarian method for assignment problems", 
    ##   Naval Research Logistics Quarterly, 3: 253-258, 1956.
    ## C. Papadimitriou and K. Steiglitz (1982), Combinatorial Optimization:
    ##    Algorithms and Complexity. Englewood Cliffs: Prentice Hall.
    ##
    ## Requires the clue library
    ##
    ## Arguments:
    ##
    ## F1:      Order Factor Loadings Matrix
    ##
    ## F2:      Input Factor Loadings Matrix
    ##
    ## Phi2:    Optional factor correlation matrix for F2 (default = NULL)
    ##
    ##
    ## MatchMethod: "LS" (Least Squares) or "CC" (congruence coefficients)
    ##
    ##
    ## Values:
    ##
    ## F2 : The re-ordered and reflected loadings of F2
    ##
    ## Phi2: Reordered and reflected factor correlations
    ##
    ## FactorMap: a 2 x  Nfac matrix structured such that:
    ##    row 1: the original column order of F2
    ##    row 2: the sorted column order of F2  
    ##
    ## UniqueMatch: A logical indicating whether a unique match was found.
    ##
    ## MatchMethod: "LS" or "CC" (congruence coefficients)
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  
    Nfac <- ncol(F1)

    UniqueMatch <- TRUE
 
    if(MatchMethod == "LS"){
  
     # Compute modified least squares (i.e., squared distance) matrix
     A <-  matrix(colSums(F1^2), Nfac, Nfac, byrow=FALSE)
     B <- t(matrix(colSums(F2^2),Nfac, Nfac, byrow=FALSE))
     
     # When factors are optimally reflected the cross product will be positive
     LSmat <- A + B  - 2 * abs(crossprod(F1,F2))
     
     # Test for unique matches
     Qmatch <- apply(LSmat,1,which.min)
     if(length(unique(Qmatch)) != Nfac){
       UniqueMatch <- FALSE
     }
   
     # If unique match not found minimize sum of squares  
     # library(clue)
     if(UniqueMatch == FALSE){  
       Qmatch <- clue::solve_LSAP(LSmat, maximum = FALSE)
     }
     
     FactorMap <- rbind(1:Nfac, Qmatch)
     
     # Find optimal reflections
     F1noZeros <- F1
     # This allows a column of F1 to have all zeros
     F1noZeros[,colSums(F1noZeros)==0]<-1
     
     Dsgn <- diag(sign(colSums( F1noZeros*F2[,Qmatch]) ))
     F2 <- F2[ , Qmatch] %*%  Dsgn

  }#End LS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  ## Match on Congruence Coefficients   
   if(MatchMethod == "CC"){
     D.F1 <- diag(1 / sqrt(apply(F1^2, 2, sum)))
     D.F2 <- diag(1 / sqrt(apply(F2^2, 2, sum)))

   # take absolute value of CC matrix
    absCosF1F2 <- abs( D.F1 %*% t(F1) %*% F2 %*% D.F2 )
   
   
   # Test for unique matches
    Qmatch <- apply(absCosF1F2, 1, which.max)
    if(length(unique(Qmatch)) != Nfac){
      UniqueMatch <- FALSE
    }
   
   # If unique match not found maximize avg CC  
    if(UniqueMatch == FALSE){  
       Qmatch <- clue::solve_LSAP(absCosF1F2, maximum = TRUE)
    }
   
    FactorMap <- rbind(1:Nfac, Qmatch)
   
   # Find optimal reflections
    Dsgn <- diag(sign(colSums( F1*F2[,Qmatch]) ))
    F2 <- F2[ , Qmatch] %*%  Dsgn
   
  
  }#End CC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  rownames(FactorMap) <- c("Original Order", "Sorted Order")

  if(!is.null(Phi2)){
    Phi2 <- Dsgn %*% Phi2[Qmatch, Qmatch] %*% Dsgn
  }
 
  ## Compute Congrence Coefficients and 
  ## RMSE for final solution
 
  FIT <- function(F1,F2){
    D.F1 <- diag(1 / sqrt(apply(F1^2, 2, sum)))
    D.F2 <- diag(1 / sqrt(apply(F2^2, 2, sum)))
   
    # take absolute value of CC matrix
    CC <- diag( D.F1 %*% t(F1) %*% F2 %*% D.F2 )
    LS <- sqrt(apply( (F1 - F2)^2 , 2, mean) )
    RMSD <- sqrt( mean( (F1 - F2)^2) ) 
    list(CC = CC, LS = LS, RMSD = RMSD)
  }#END FIT
 
  fit <- FIT(F1,F2)

  list(F2 = F2,
       Phi2 = Phi2,
       FactorMap = FactorMap,
       UniqueMatch = UniqueMatch,
       MatchMethod = MatchMethod,
       CC = fit$CC,
       LS = fit$LS,
       RMSD = fit$RMSD)
  } # faAlign2
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

  # ---- Fnc to compute  Hyperplane Count
  Hyperplane <- function(Fload, HPthreshold) {
    sum(apply(abs(Fload) < HPthreshold, 2, sum))
  } # END Hyperplane   
  
  
  ConvergenceStatus <- sapply(fout$localSolutionSets[[SolSet]],
                                     "[[", 5)
  
  
  if(length(ConvergenceStatus) == 1){
     warning("\n\n\n **** FATAL ERROR: Only one solution in this solution set.****\n\n")
     invisible(return("error"))
    }
 
  HPcount <- sapply(lapply(fout$localSolutionSet[[SolSet]], "[[", 1),
                    Hyperplane, HPthreshold)
  
  
  ## Compute number of local minima
  num_in_Setk <-length(sapply(fout$localSolutionSet[[SolSet]], "[[", 3))
  complexity.val <- sapply(fout$localSolutionSet[[SolSet]], "[[", 3)


  #Construct results table
  len_rmsd <- num_in_Setk * (num_in_Setk - 1) / 2

  rmsdTable <- matrix(99, len_rmsd, 3)
  
  
  #Calculate rmsd values
  r <- 1
  for(i in 1: (num_in_Setk - 1) ){
    for(j in (i+1):num_in_Setk){
      rmsd_ij<- faAlign2(fout$localSolutionSet[[SolSet]][[i]]$loadings,
                         fout$localSolutionSet[[SolSet]][[j]]$loadings)$RMSD
      rmsdTable[r, 1] <- i
      rmsdTable[r, 2] <- j
      
      if(ConvergenceStatus[i]==FALSE | ConvergenceStatus[j]==FALSE){
        rmsd_ij <- 999
      }
      rmsdTable[r, 3] <- rmsd_ij
      r <- r+1
      
    }#END for j
  } #END for (i
  
  # Select converged solutions
  
  rmsdTable <- rmsdTable[sort.list(rmsdTable[,3], 
                                   decreasing = TRUE), ,drop=FALSE ]
 
  colnames(rmsdTable) <- list("Solution.i", "Solution.j", "RMSD")
  rmsdTable[, 3] <- round(rmsdTable[, 3], digits)
  
  # Print
  if(PrintLevel == 1){
    cat("\n\n Solution Set", SolSet, "\n",
        "Complexity Value", round(complexity.val[1], digits),"\n",
        "The 6 most discrepant solutions","\n\n")
    print( utils::head(rmsdTable))
  }  
  
  if(PrintLevel == 2){
    cat("\n\n Solution Set", SolSet, "\n",
        "Complexity Value", round(complexity.val[1], digits),"\n\n" )
    print( rmsdTable )
  }  
  
  
  invisible(list("rmsdTable" = rmsdTable,
       "Set" = SolSet,
       "complexity.val" = round(complexity.val[1], digits),
       "HPcount" = HPcount,
       "converged" = ConvergenceStatus)
  )     
  
} #END faLocalMin  
  


