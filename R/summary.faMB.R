#' Summary Method for an Object of Class faMB
#' 
#' This function summarizes results from a call to \strong{faMB}.
#' 
#' @param object (Object of class \code{\link{faMB}}) The returned object 
#' from a call to \strong{faMB}.
#' @param digits (Integer) Print output with user-specified number of significant digits. 
#'   Default \code{digits = 2}.   
#' @param Set The argument \code{Set} can be specified as either an integer 
#'     value (i.e., 1 through the number of unique solution sets) or a character 
#'     value (i.e., 'UnSpun'). 
#'   \itemize{
#'     \item{\strong{Integer}} Summarize the solution from the specified 
#'     solution set. If \code{Set = 1}, the "global minimum" solution is 
#'     reported. See \code{\link{faMain}} for more details about finding 
#'     the "global" and local minima. 
#'     \item{\strong{'UnSpun'}} Summarize the solution from the rotated 
#'     output that was produced by rotating from the unrotated (i.e., unspun) 
#'     factor orientation. All other solutions are rotated from a randomly 'spun' rotation 
#'     (i.e., by orientating the unrotated factor solution via a random orthonormal 
#'     matrix) . 
#'    }
#' @param HPthreshold (Numeric) User-defined threshold for declaring that the 
#'  absolute value of a factor pattern coefficient is in a hyperplane. The hyperplane count is the number of 
#'  near-zero (as defined by HPthreshold; see Cattell, 1978, p. 105) elements in the factor pattern matrix. 
#'  Default \code{HPthreshold = .05}.
#' @param  PrintLevel (Integer) Controls the level of printing. If \code{PrintLevel = 0} then no output is printed. 
#' If \code{PrintLevel = 1} then the standard output 
#' will be printed. If \code{PrintLevel = 2} more extensive output (e.g., the Factor Structure Matrix, 
#' the Residuals Matrix [i.e., Observed - fitted R]) will 
#' be printed. Default \code{PrintLevel = 1}. 
#' @param DiagnosticsLevel (Integer) Controls the amount of diagnostics information that is computed on the 
#' rotation local minima. If \code{DiagnosticsLevel = 1} then only the number 
#' of local solution sets will be reported. If \code{DiagnosticsLevel = 2} then
#' the program will determine whether all solutions within a solution set are identicial.
#' Default \code{DiagnosticsLevel = 1}.
#' @param \dots Additional arguments affecting the summary produced. 
#' 
#' @details \strong{summary.faMB} provides various criteria for judging the adequacy of 
#' the rotated factor solution(s). After reporting the number of solution sets.
#' (i.e., rotated solutions with the same complexity value) the following measures 
#'   of factor adequacy are reported for each solution set:
#' \itemize{
#'  \item \strong{Complexity Value}: The rotation complexity value (see \code{\link{faMain}} for details).
#'  \item \strong{Hyperplane Count}: The number of near-zero loadings (defined by \strong{HPthreshold}) 
#'    for all factor patterns in a solution set (if \strong{MaxWithinSetRMSD > 0} then Hyperplane Count refers to 
#'  the first factor pattern in the solution set). 
#'  \item \strong{\% Cases (x 100) in Set}: The percentage of factor patterns in each solution set.
#'  \item \strong{RMSD}: The root mean squared deviation between the first factor pattern 
#'    in each solution set with the first factor pattern  in the solution set specified by the \strong{Set} parameter. By default, \strong{Set = 1}.
#'  \item \strong{MaxWithinSetRMSD}: The maximum root mean squared deviation between all within set solutions and 
#'  the first element in the solution set. When \strong{MaxWithinSetRMSD > 0} then the solution 
#'  set contains non-identical rotated factor patterns with identical complexity values. 
#'  \item \strong{Converged}: A Logical (TRUE/FALSE) that  indicates whether all within set rotations converged.
#'  }
#' 
#' @return
#'  \itemize{
#'     \item \code{loadings} (Matrix) Factor loadings for the solution associated with the 
#'                    minimum (maximum) rotation complexity value (default) or the user-chosen solution.
#'      \item \code{Phi} (Matrix) Factor correlation matrix for the solution associated with the 
#'                    minimum (maximum) rotation complexity value (default) or the user-chosen solution.
#'      \item \code{FS} (Matrix) Factor structure matrix  for the solution associated with the 
#'                    minimum (maximum) rotation complexity value (default) or the user-chosen solution.
#'      \item \code{Set} (Integer) The returned Set number. 
#'      \item \code{facIndeterminacy} (Matrix) Factor Indeterminacy values. 
#'      \item \code{SetComplexityValues} (vector) Rotation complexity value for each solution set. 
#'      \item \code{HP_counts} (vector) Hyperplane count for each solution set.  
#'      \item \code{MaxWithinSetRMSD} (vector) If \code{DiagnosticsLevel = 2} the the program will compute
#'      within set RMSD values.  These values represent the root mean squared deviations of each 
#'      within set solution with the first solution in a set. If the \code{MaxWithinSetRMSD = 0} 
#'      for a set, then all within set solutions are identical. If  \code{MaxWithinSetRMSD > 0} 
#'      then at least one solution differs from the remaining solutions within a set (i.e., two solutions 
#'      with different factor loadings produced identical complexity values). 
#'      \item \code{ChiSq} (Numeric) Chi-square goodness of fit value. As recommended by Browne (1979), 
#'      we apply Lawley's (1959) correction when computing the chi-square value when \code{NB = 2}.
#'      \item \code{DF} (Numeric) Degrees of freedom for the estimated model. 
#'      \item \code{pvalue} (Numeric) P-value associated with the above chi-square statistic.
#'      \item \code{AIC} (Numeric) Akaike's Information Criterion where a lower value indicates better fit. 
#'      \item \code{BIC} (Numeric) Bayesian Information Criterion where a lower value indicates better fit. 
#'      \item \code{RMSEA} (Numeric) The root mean squared error of approximation (Steiger & Lind, 1980).
#'      \item \code{Resid} (Matrix) The residuals matrix (R - Rhat). 
#'      \item \code{NumberLocalSolutions} (Integer) The number of local solution sets.     
#'      \item \code{LocalSolutions} (List) A list of local solutions (factor loadings, factor correlations, etc). 
#'      \item\code{rotate} Designates which rotation method was applied.
#'    }
#'     
#' @references   Cattell, R. (1978). The scientific use of factor analysis in behavioral and life sciences. 
#'  New York, New York, Plenum. 
#'     
#' @family Factor Analysis Routines
#' 
#' @author
#' \itemize{
#'   \item Niels G. Waller (nwaller@umn.edu)
#'   \item Casey Giordano (Giord023@umn.edu)
#'   }
#'   
#' @examples
#' # These examples reproduce published multiple battery analyses. 
#' 
#' # ----EXAMPLE 1: Browne, M. W. (1979)----
#' #
#' # Data originally reported in:
#' # Thurstone, L. L. & Thurstone, T. G. (1941). Factorial studies 
#' # of intelligence. Psychometric Monograph (2), Chicago: Univ. 
#' # Chicago Press.
#' 
#' ## Load Thurstone & Thurstone's data used by Browne (1979)
#' data(Thurstone41)
#' 
#' Example1Output <-  faMB(R             = Thurstone41, 
#'                         n             = 710,
#'                         NB            = 2, 
#'                         NVB           = c(4,5), 
#'                         numFactors    = 2,
#'                         rotate        = "oblimin",
#'                         rotateControl = list(standardize = "Kaiser"))
#'                         
#' ## Call the summary function
#' summary(Example1Output)
#' 
#' # ----EXAMPLE 2: Browne, M. W. (1980)----
#' # Data originally reported in:
#' # Jackson, D. N. & Singer, J. E. (1967). Judgments, items and 
#' # personality. Journal of Experimental Research in Personality, 20, 70-79.
#' 
#' ## Load Jackson and Singer's dataset
#' data(Jackson67)
#' 
#' Example2Output <-  faMB(R             = Jackson67, 
#'                         n             = 480,
#'                         NB            = 5, 
#'                         NVB           = rep(4,5), 
#'                         numFactors    = 4,
#'                         rotate        = "varimax",
#'                         rotateControl = list(standardize = "Kaiser"),
#'                         PrintLevel    = 1)
#' 
#' ## Call the summary function
#' summary(object     = Example2Output,
#'         Set        = 1,
#'         PrintLevel = 1)
#' 
#' # ----EXAMPLE 3: Cudeck (1982)----
#' # Data originally reported by:
#' # Malmi, R. A., Underwood, B. J., & Carroll, J. B. (1979).
#' # The interrelationships among some associative learning tasks. 
#' # Bulletin of the Psychonomic Society, 13(3), 121-123. DOI: 10.3758/BF03335032 
#' 
#' ## Load Malmi et al.'s dataset
#' data(Malmi79)
#' 
#' Example3Output <- faMB(R             = Malmi79, 
#'                        n             = 97,
#'                        NB            = 3, 
#'                        NVB           = c(3, 3, 6), 
#'                        numFactors    = 2,
#'                        rotate        = "oblimin",
#'                        rotateControl = list(standardize = "Kaiser"))
#' 
#' ## Call the summary function
#' summary(object     = Example3Output,
#'         Set        = 1,
#'         PrintLevel = 2)
#'         
#' # ----Example 4: Cudeck (1982)----
#' # Data originally reported by: 
#' # Boruch, R. F., Larkin, J. D., Wolins, L. and MacKinney, A. C. (1970). 
#' #  Alternative methods of analysis: Multitrait-multimethod data. Educational 
#' #  and Psychological Measurement, 30,833-853.
#' 
#' ## Load Boruch et al.'s dataset
#' data(Boruch70)
#' 
#' Example4Output <- faMB(R             = Boruch70,
#'                        n             = 111,
#'                        NB            = 2,
#'                        NVB           = c(7,7),
#'                        numFactors    = 2,
#'                        rotate        = "oblimin",
#'                        rotateControl = list(standardize  = "Kaiser",
#'                                             numberStarts = 100))
#'
#' ## Call the summary function
#' summary(Example4Output)
#' 
#' @export

summary.faMB <- function(object, 
                           digits           = 2, 
                           Set              = 1, 
                           HPthreshold      = .05, 
                           PrintLevel       = 1, 
                           DiagnosticsLevel = 1,
                           ... ) {
  
  # is the object of the appropriate class
  stopifnot(inherits(object, "faMB"))
  
  ## Make copy of output to not disturb original object
  fout <- object
  

  maxRmsdWithFi <- NULL
  
  
  #----fnc to Print Section Headings 
  sectionHeading <- function(string){
    stringLength <- nchar(string)
    cat("\n\n")
    cat(c(paste(rep("=", stringLength ),collapse =""),"\n"))
    cat(c(string, "\n"))
    cat(c(paste(rep("=", stringLength ),collapse =""),"\n"))
  } # END sectionHeading
  
  
  ## ---Count number of local solution sets
  NumberLocalSolutions <- fout$numLocalSets
  
  
  ## -- Input Error Checking----
  if( is.numeric(Set) && Set > (NumberLocalSolutions) ){
    stop("\n\nFATAL ERROR: Invalid value for Set")
  } 
  
  # ---- Fnc to compute  Hyperplane Count
  Hyperplane <- function(Fload, salient = .05) {
    sum(apply(abs(Fload) < salient, 2, sum))
  } # END Hyperplane   
  
  
  ## ---- Gather local Solutions----
  LocalSolutionsList <- list()
  for (iSol in 1: NumberLocalSolutions) {
    #Pick first element in solution set
    LocalSolutionsList[[iSol]] <- c(iSol,
                                    fout$localSolutionSets[[iSol]][[1]])
    names(LocalSolutionsList[[iSol]]) <- c("Set", 
                                           names(fout$localSolutionSets[[iSol]][[1]]))
  } # END for(iSol in 1: NumberLocalSolutions) 
  
  # UnSpun Solution is last element in list
  LocalSolutionsList[[NumberLocalSolutions + 1]] <- c(0, fout$unSpunSolution)
  names(LocalSolutionsList[[NumberLocalSolutions + 1]])<-c("Set", names(fout$unSpunSolution))
  ## END Gather local Solutions
  
  SetArgument <- Set
  
  # ---- Which solution should be viewed?-----
  if(Set == "UnSpun") Set <- NumberLocalSolutions + 1
  
  F_view <- LocalSolutionsList[[Set]]$loadings
  Phi_view <- LocalSolutionsList[[Set]]$Phi
  
  
  
  # Align selected solution (which may not be the 
  # min rotation solution) with the fout$loadings from faMB
  
  if(ncol(fout$loadings) > 1) {
     ## Align selected solution with global minimum (easier comparison)
     F_viewAlignOut <- faAlign(F1          = fout$loadings,
                            F2          = F_view, 
                            Phi2        = Phi_view,
                            MatchMethod = "LS")
     F_view <- F_viewAlignOut$F2
     Phi_view <- F_viewAlignOut$Phi2
     Fac_indeter_view <- LocalSolutionsList[[Set]]$facIndeterminacy[F_viewAlignOut$FactorMap[2, ]]
  
     dimnames(F_view) <- dimnames(fout$loadings)
     dimnames(Phi_view) <- dimnames(fout$Phi)
     names(Fac_indeter_view) <- colnames(F_view)
  }else{
     F_view <- fout$loadings
     Phi_view <- fout$Phi
     Fac_indeter_view <-fout$facIndeterminacy
  }  #ENDif(ncol(fout$loadings) > 1) 
  
  
  RotationCoverged <- 
    LocalSolutionsList[[Set]]$RotationConverged
  # ---- END Which solution should be viewed?-----
  
  
  #Compute FI for chosen solution
  if(max(Fac_indeter_view) > 1){
    Fac_indeter_view <- rep(NA, ncol(fout$loadings))
  }
  
  # # If full targetMatrix supplied then align viewed solution 
  # if ( !is.null(fout$targetMatrix) && sum(is.na(fout$targetMatrix)) == 0) {
  #   alignOut <- faAlign(F1   = fout$targetMatrix,
  #                       F2   = F_view, 
  #                       Phi2 = Phi_view)
  #   F_view <- alignOut$F2
  #   Phi_view <- alignOut$Phi
  #   dimnames(F_view) <- dimnames(fout$loadings)
  #   dimnames(Phi_view) <- dimnames(fout$Phi)
  # } #END if( !is.null(fout$targetMatrix) && sum(is.na(fout$targetMatrix))==0) 
   
  
  # ----Compute Structure Matrix
  FS <- F_view %*% Phi_view
  
  # ---- Find set rotation complexity values
  SetComplexityValues <- 
    as.numeric(lapply(LocalSolutionsList, "[[", "RotationComplexityValue"))
  
  
  # ----Find convergence status 
  SetConvergeStatus <- lapply(LocalSolutionsList, "[[", "RotationConverged")
  
  
  
  # ---- Compute rmsd with desired set solution ----
  rmsdVec <- rep(0, NumberLocalSolutions)
  if(ncol(fout$loadings) > 1){
       for (iSet in 1:(NumberLocalSolutions + 1)) {
          # align local solution with user chosen view solution
            F_Set_Aligned <- faAlign(F1 = F_view, 
                             F2 = LocalSolutionsList[[iSet]]$loadings)$F2
    
             rmsdVec[iSet] <-  rmsd(A = F_view, 
                           B = F_Set_Aligned,
                           Symmetric = FALSE,
                           IncludeDiag = TRUE)
         } # END for (iSet in 1:(NumberLocalSolutions + 1)) 
  }#END  if(ncol(fout$loadings) > 1){   
  
  ##-------END rmsd Vec
  
  
  #----Calculate hyperplane counts----
  HP_counts <- rep(0, NumberLocalSolutions + 1)
  for (iSet in 1:(NumberLocalSolutions + 1)) {
    HP_counts[iSet] <- Hyperplane(Fload   = LocalSolutionsList[[iSet]]$loadings,
                                  salient = HPthreshold) 
  }##End count hyperplane
  
  
  #---- Num Cases in Solution Sets----
  NumCasesInSolSets <- as.numeric(lapply(fout$localSolutionSets,length))
  
  # ----Check Set Homogeneity----   
  # Are all solutions identicial within Solution Sets
  if ( DiagnosticsLevel == 2 ) {
    
    ## Check RMSD within a Set (i.e., are all loadings matrices equivalent)
    rmsdWithF1 <- function(FList, k){
      ## Ensure all solutions in a set are aligned to the first element in set
      Fk <- faAlign(F1 = FList[[1]]$loadings, 
                    F2 = FList[[k]]$loadings)$F2
      ## Compute RMSD between 1st and kth element in local solution set
      rmsd(A           = FList[[1]]$loadings, 
           B           = Fk, 
           Symmetric   = FALSE,
           IncludeDiag = TRUE)
    } # END rmsdWithF1
    
    
    ## Find the maximum within-set RMSD between 1st and kth loadings matrices
    maxRmsdWithFi <- rep(99, NumberLocalSolutions)
    for (iSet in 1:NumberLocalSolutions) {
      
      ## Declare container for rmsd values
      rmsdWithF1Out <- rep(-99, NumCasesInSolSets[iSet] )
      
      for (jCases in 1:NumCasesInSolSets[iSet]) {
        
 
        # ## For each solution within a set, order the items
        # for (kSort in 1:NumCasesInSolSets[iSet]) {
        #   fout$localSolutionSets[[iSet]][[kSort]]$loadings <- 
        #     fout$localSolutionSets[[iSet]][[kSort]]$loadings[itemOrder, ]
        # } # END for (kSort in 1:NumCasesInSolSets[iSet]) 
        
        ## compute RMSD between sorted solutions
        rmsdWithF1Out[jCases] <- rmsdWithF1(fout$localSolutionSets[[iSet]], jCases) 
        
      } # END for (jCases in 1:NumCasesInSolSets[iSet]) 
      
      ## Extract the maximum within set RMSD value 
      maxRmsdWithFi[iSet] <-  round(max(rmsdWithF1Out), 3)
    } # END for (iSet in 1:NumberLocalSolutions) 
    
  } # END Check Set Homogeneity
  
  
  
  #---Generate Rotation Fit Table
  if (DiagnosticsLevel == 1) {
    Rotation_Fit_Table <- 
      rbind( round(SetComplexityValues, digits + 1),
             HP_counts,
             c(100 * round(NumCasesInSolSets/sum(NumCasesInSolSets), 2), 0), 
             round(rmsdVec, 3),
             SetConvergeStatus)
    
    
    colnames( Rotation_Fit_Table ) <- 
      c(paste0("Set ", seq(1:NumberLocalSolutions), "  "), "UnSpun")
    rownames( Rotation_Fit_Table ) <- c("Complexity Value", 
                                        "Hyperplane Count", 
                                        "% Cases (x 100) in Set",
                                        "RMSD",
                                        "Converged")
    
  }# END if (DiagnosticsLevel == 1)
  
  if(DiagnosticsLevel == 2){
    Rotation_Fit_Table <- 
      rbind( round(SetComplexityValues, digits + 1),
             HP_counts,
             c(100 * round(NumCasesInSolSets/sum(NumCasesInSolSets), 2), 0), 
             round(rmsdVec, 3),
             c(maxRmsdWithFi, 0),
             SetConvergeStatus)
    
    colnames( Rotation_Fit_Table ) <- 
      c(paste0("Set ", seq(1:NumberLocalSolutions), "  "),"UnSpun")
    rownames( Rotation_Fit_Table ) <- c("Complexity Value", 
                                        "Hyperplane Count", 
                                        "% Cases (x 100) in Set",
                                        "RMSD",
                                        "MaxWithinSetRMSD",
                                        "Converged")
  }# END if(DiagnosticsLevel == 2)
  

    Rhat <- fout$Rhat
    ResidSave <- Resid <- fout$Resid
    rownames(Resid) <- colnames(Resid) <-  paste0("v", 1:ncol(Resid))
   
  
  
  ## ---- Print Output ----
  #----____PrintLevel == 0----
  if(PrintLevel == 0){  # No Output printed

    
    ##  RETURN VALUES
    return(invisible(
      list(loadings             = F_view,
           Phi                  = Phi_view,
           FS                   = FS,
           Set                  = Set,
           facIndeterminacy     = Fac_indeter_view,
           SetComplexityValues  = SetComplexityValues,
           HP_counts            = HP_counts,
           MaxWithinSetRMSD     = maxRmsdWithFi,
           ChiSq                = fout$fit$ChiSq,
           DF                   = fout$fit$DF,
           pvalue               = fout$fit$pvalue,
           AIC                  = fout$fit$AIC, 
           BIC                  = fout$fit$BIC,
           RMSEA                = fout$fit$RMSEA,
           Resid                = ResidSave,
           NumberLocalSolutions = NumberLocalSolutions,
           LocalSolutions       = LocalSolutionsList,
           rotate               = fout$rotate)))
  } #END if(PrintLevel == 0)
  
  # Print faMB function call
  cat("\nCall:\n")
  print(fout$Call)
  cat("\n")
  
  # Rotation information:
  sectionHeading("Rotation Diagnostics")
  
  if(ncol(fout$loadings) == 1 || fout$rotate =="none"){
       fout$rotateControl$numberStarts <- 0
  }     
  
  #----Print Rotation_Fit_Table----
  cat("\nRotation Method:", fout$rotate)
  cat("\nEpsilon value for rotation convergence and solution set clustering:", fout$rotateControl$epsilon) 
  cat("\nNumber of Rotation Random Starts:", fout$rotateControl$numberStarts)
  cat("\nNumber of Solution Sets:", NumberLocalSolutions)
  cat("\nSolution Characteristics:\n")
  cat("\n")
  print( Rotation_Fit_Table,
         right = TRUE)
  
  
  # ----____Print factor output----  
  cat(paste0("\nFactor Output for Set " , SetArgument, ":\n"))  
  cat("Factor Pattern Matrix\n\n")
  
  ## Print rounded factor loadings matrix
  print( round(F_view, digits))
  
  cat(paste0("\n\nFactor Correlation Matrix\n"))
  cat("\n")
  print( round( Phi_view, digits))
  
  #----____Print Structure Matrix----
  # PrintLevel == 2:  Sturcture Matrix
  if(PrintLevel == 2){
    cat("\n\nFactor Structure Matrix\n")
    print( round(FS, digits))
  }
  
  
  #----____Print Factor Indeterminancy----  
  
  sectionHeading("Factor Indeterminacy")
  
  Lmin <- eigen(fout$R)$values[nrow(FS)]
  if(Lmin > 1E-6){
    # 
      print( round(Fac_indeter_view, digits) )
      cat("\n")
  }# END if(Lmin > 1E-6)
  
  if(Lmin < 1E-6){
    cat("\nObserved R computationally singular. Unable to calculate FI values.\n")
    Fac_indeter_view <- rep(NA, ncol(fout$loadings))
  }
  
  #----____Print Model Fit----
  #Model Fit:
  sectionHeading("Model Fit")
  
  cat("\nSample Size:       ",   fout$fit$SampleSize)
  cat("\nChi square:        ",   round(fout$fit$ChiSq, digits = digits + 1))
  cat("\nDegrees of freedom:",   round(fout$fit$DF, digits = digits + 1))
  cat("\np-value:           ",   round(fout$fit$pvalue, digits = digits + 1))
  cat("\nAIC:               ",   round(fout$fit$AIC, digits = digits + 1))
  cat("\nBIC:               ",   round(fout$fit$BIC, digits = digits + 1))
  cat("\nRMSEA:             ",   round(fout$fit$RMSEA, digits = digits + 1))
  
  
  # END PrintLevel == 1 output 
  
  if(PrintLevel == 2){
    cat("\n\nResiduals  (R - Rhat)\n")
    if(ncol(Resid) < 16){
      #print actual labels if it will fit on screen
      print(round( ResidSave, 2))
    }else{
     print(round( Resid, 2))
    }
  }  
  
  cat("\n\n")
  
  ##  RETURN VALUES
  invisible(list(loadings = F_view,
                 Phi = Phi_view,
                 FS = FS,
                 Set = SetArgument,
                 facIndeterminacy = Fac_indeter_view,
                 SetComplexityValues = SetComplexityValues,
                 HP_counts = HP_counts,
                 MaxWithinSetRMSD = maxRmsdWithFi,
                 ChiSq = fout$fit$ChiSq,
                 DF = fout$fit$DF,
                 pvalue = fout$fit$pvalue,
                 AIC = fout$fit$AIC, 
                 BIC = fout$fit$BIC,
                 RMSEA = fout$fit$RMSEA,
                 Resid = ResidSave,
                 NumberLocalSolutions = NumberLocalSolutions,
                 LocalSolutions = LocalSolutionsList))
  
  
} ## END summary.faMB 


