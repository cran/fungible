#' Summary Method for an Object of Class faMain
#' 
#' This function summarizes results from a call to \strong{faMain}.
#' 
#' @param object (Object of class \code{\link{faMain}}) The returned object 
#' from a call to \strong{faMain}.
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
#' will be printed. If \code{PrintLevel = 2} more extensive output (e.g., the Factor Structure Matrix) will 
#' be printed. Default \code{PrintLevel = 1}. 
#' @param DiagnosticsLevel (Integer) Controls the amount of diagnostics information that is computed on the 
#' rotation local minima. If \code{DiagnosticsLevel = 1} then only the number 
#' of local solution sets will be reported. If \code{DiagnosticsLevel = 2} then
#' the program will determine whether all solutions within a solution set are identicial.
#' Default \code{DiagnosticsLevel = 1}.
#' @param itemSort (Logical) If TRUE, sort the order of the observed variables to produce
#' a "staircase"-like pattern. Note that this argument cannot handle bifactor models at this time.
#' Defaults to \code{itemSort} = FALSE.
#' @param \dots Additional arguments affecting the summary produced. 
#' 
#' @details \strong{summary.faMain} provides various criteria for judging the adequacy of 
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
#'  \item \strong{Converged}: A Logical (TRUE/FALSE) that  indicates whether the first solution 
#'  in a solution set  has a TRUE convergence status. 
#'  }
#'   Note that  the printed factor pattern is not sorted even if \strong{itemSort} 
#'  is requested in \link{faMain}.
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
#'      \item \code{h2} (Matrix) Communalities for the returned factor solution. If \code{Boostrap = TRUE} then
#'      \code{h2} also returns the bootstrap standard errors and associated confidence bounds from 
#'                the bootstrap distribution.
#'      \item{facIndeterminacy} (Vector) Factor Indeterminacy values (correlations between the factors and factor scores). If \code{Boostrap = TRUE} then
#'         \code{facIndeterminacy} also returns the bootstrap standard errors and associated confidence bounds from 
#'                the boostrap distribution.          
#'      \item \code{SetComplexityValues} (Vector) Rotation complexity value for each solution set. 
#'      \item \code{HP_counts} (Vector) Hyperplane count for each solution set.  
#'      \item \code{MaxWithinSetRMSD} (Vector) If \code{DiagnosticsLevel = 2} the the program will compute
#'      within set RMSD values.  These values represent the root mean squared deviations of each 
#'      within set solution with the first solution in a set. If the \code{MaxWithinSetRMSD = 0} 
#'      for a set, then all within set solutions are identical. If  \code{MaxWithinSetRMSD > 0} 
#'      then at least one solution differs from the remaining solutions within a set (i.e., two solutions 
#'      with different factor loadings produced identical complexity values). 
#'      \item \code{RMSD} (Numeric) The root mean squared deviation between the 
#'           observed and model-implied correlation matrix.
#'      \item \code{RMSAD} (Numeric) The root mean squared absolute deviation between the 
#'           observed and model-implied correlation matrix. 
#'      \item \code{NumberLocalSolutions} (Integer) The number of local solution sets.     
#'      \item \code{LocalSolutions} (List) A list of local solutions (factor loadings, factor correlations, etc). 
#'      \item\code{rotate} Designates which rotation method was applied.
#'      \item\code{itemOrder} The item order of the (possibly) sorted factor loadings.
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
#' ## Load Thurstone's Box data from the fungible library
#' library(fungible)
#' data(Box26)
#'
#' ## Create a matrix from Thurstone's solution
#' ## Used as a target matrix to sort columns of the estimated solution
#' ThurstoneSolution <- matrix(c(   .95,  .01,  .01,
#'                                  .02,  .92,  .01,
#'                                  .02,  .05,  .91,
#'                                  .59,  .64, -.03,
#'                                  .60,  .00,  .62,
#'                                 -.04,  .60,  .58,
#'                                  .81,  .38,  .01,
#'                                  .35,  .79,  .01,
#'                                  .79, -.01,  .41,
#'                                  .40, -.02,  .79,
#'                                 -.04,  .74,  .40,
#'                                 -.02,  .41,  .74,
#'                                  .74, -.77,  .06,
#'                                 -.74,  .77, -.06,
#'                                  .74,  .02, -.73,
#'                                 -.74, -.02,  .73,
#'                                 -.07,  .80, -.76,
#'                                  .07, -.80,  .76,
#'                                  .51,  .70, -.03,
#'                                  .56, -.04,  .69,
#'                                 -.02,  .60,  .58,
#'                                  .50,  .69, -.03,
#'                                  .52, -.01,  .68,
#'                                 -.01,  .60,  .55,
#'                                  .43,  .46,  .45,
#'                                  .31,  .51,  .46), nrow = 26, ncol = 3,
#'                                                             byrow=TRUE)
#' ## Example 1: Multiple solution sets.
#' ## Ignore warnings about non-positive definite sample correlation matrix
#' suppressWarnings(
#'   fout <- faMain(R             = Box26,
#'                  numFactors    = 3,
#'                  facMethod     = 'faregLS',
#'                  rotate        = 'infomaxQ',
#'                  targetMatrix  = ThurstoneSolution,
#'                  rotateControl = 
#'                    list(numberStarts = 25, ## increase in real problem
#'                         standardize  = 'none'),
#'                  Seed          = 123)
#' )
#' 
#' ## Summarize the factor analytic output                                     
#' summary(object           = fout, 
#'         digits           = 2,
#'         Set              = 2, 
#'         HPthreshold      = .10,
#'         PrintLevel       = 1,
#'         DiagnosticsLevel = 2)
#'           
#'      
#' ## Example 2: Bootstrap Illustration 
#' ## Step 1: In an initial analysis, confirm that all rotations converge
#'   ## to a single minimum complexity value.
#' ## Step 2: If Step 1 is satisfied then generate bootstrap samples.
#' 
#' ## Load Amazon box data             
#' data("AmzBoxes")
#' 
#' ## Convert box dimensions into Thurstone's indicators
#' BoxData <- 
#'   GenerateBoxData(AmzBoxes[, 2:4],          ## Select columns 2, 3, & 4
#'                   BoxStudy         = 26,    ## 26 indicators
#'                   Reliability      = 0.75,  ## Add unreliability
#'                   SampleSize       = 200,   ## Add sampling error
#'                   ModApproxErrVar  = 0.1,   ## Add model approx error
#'                   NMinorFac        = 50,    ## Number of minor factors
#'                   epsTKL           = 0.2,   ## Spread of minor factor influence
#'                   SeedErrorFactors = 1,     ## Reproducible starting seed
#'                   SeedMinorFactors = 2,     ## Reproducible starting seed
#'                   PRINT            = FALSE, ## Suppress some output
#'                   LB               = FALSE, ## Do not set lower-bounds
#'                   LBVal            = 1,     ## Lower bound value (ignored)
#'                   Constant         = 0)     ## Do not add constant to data
#'                            
#' ## Analyze new box data with added measurement error
#' fout <- faMain(X             = BoxData$BoxDataE,
#'                numFactors    = 3,
#'                facMethod     = 'fapa',
#'                rotate        = 'infomaxQ',
#'                targetMatrix  = ThurstoneSolution,
#'                bootstrapSE   = FALSE,
#'                rotateControl = 
#'                  list(numberStarts = 25, ## increase in real problem
#'                       standardize  = 'CM'),
#'                Seed          = 1)
#'                
#' ## Summarize factor analytic output                
#' sout <- summary(object     = fout, 
#'                 Set        = 1,
#'                 PrintLevel = 1)
#'                 
#' ## Generate bootstrap samples
#' fout <- faMain(X             = BoxData$BoxDataE,
#'                numFactors    = 3,
#'                facMethod     = 'fapa',
#'                rotate        = 'infomaxQ',
#'                targetMatrix  = ThurstoneSolution,
#'                bootstrapSE   = TRUE,
#'                numBoot       = 25,   ## increase in real problem
#'                rotateControl = 
#'                  list(numberStarts = 1,
#'                       standardize  = 'CM'),
#'                Seed          = 1)
#'
#' ## Summarize factor analytic output with bootstraps
#' sout <- summary(object     = fout, 
#'                 Set        = 1,
#'                 PrintLevel = 2)  
#'                   
#'  ## To print a specific solution without computing diagnostics and 
#'    ## summary information, use the print function.
#'  
#'    print(fout, 
#'          Set = 1)                 
#'  
#' @export

summary.faMain <- function(object, 
                           digits           = 2, 
                           Set              = 1, 
                           HPthreshold      = .05, 
                           PrintLevel       = 1, 
                           DiagnosticsLevel = 1,
                           itemSort         = FALSE, 
                           ... ) {
  
  ## TODO list
  ## add bootstrap results for h2
  ## add bootstrap results for factor indeterminancy
  
  
  
  # is the object of the appropriate class
  stopifnot(inherits(object, "faMain"))
  
  ## Make copy of output to not disturb original object
  fout <- object
  
  if(fout$rotate == "none"){
    stop("\n\n**FATAL ERROR: summary.faMain only applies to rotated, multiple factor models\n\n")
  }  
  
  ## If itemSort is false, keep original order of items
  ## See "Sort the items" for creating the itemOrder if itemSort == TRUE
  if (itemSort == FALSE) itemOrder <- seq_len(nrow(fout$loadings))
  
  maxRmsdWithFi <- NULL
  
  ## Extract confidence level for future printing
  CILevel <- round(fout$CILevel * 100, 2)
  
  # ---urLoadings = TRUE----
  #If faMain used only for rotation then set onlyLoadings to TRUE
  onlyLoadings <- TRUE
  # if "urLoadings" not in Call then all output available
  if (length(grep("urLoadings", names(fout$Call))) == 0) onlyLoadings <- FALSE
  
  ## CG note (13 sept 19): grep is > 10 times slower than is.null
  ## consider adding urLoadings as output from faMain
  # microbenchmark::microbenchmark(grep("urLoadings", names(fout$Call)),
  #                                is.null(fout$urLoadings),
  #                                times = 10000)
  
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
  
  
  ## -- Input Error Checking
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
  LocalSolutionsList[[NumberLocalSolutions+ 1]] <- c(0, fout$unSpunSolution)
  names(LocalSolutionsList[[NumberLocalSolutions+ 1]])<-c("Set", names(fout$unSpunSolution))
  ## END Gather local Solutions
  
  SetArgument <- Set
  # ---- Which solution should be viewed?-----
  if(Set == "UnSpun") Set <- NumberLocalSolutions + 1
  F_view <- LocalSolutionsList[[Set]]$loadings
  Phi_view <- LocalSolutionsList[[Set]]$Phi
  dimnames(F_view) <- dimnames(fout$loadings)
  dimnames(Phi_view) <- dimnames(fout$Phi)
  
  RotationCoverged <-
    LocalSolutionsList[[Set]]$RotationConverged
  
  # ---- END Which solution should be viewed?-----
  
  
  # If full targetMatrix supplied then align viewed solution 
  if ( !is.null(fout$targetMatrix) && sum(is.na(fout$targetMatrix)) == 0) {
    alignOut <- faAlign(F1   = fout$targetMatrix,
                        F2   = F_view, 
                        Phi2 = Phi_view)
    F_view <- alignOut$F2
    Phi_view <- alignOut$Phi
    dimnames(F_view) <- dimnames(fout$loadings)
    dimnames(Phi_view) <- dimnames(fout$Phi)
  } #END if( !is.null(fout$targetMatrix) && sum(is.na(fout$targetMatrix))==0) 
  
  ## ---- Sort the items ----
  
  ## Determine if a bifactor rotation was applied
  isBifactor <- FALSE
  if(fout$rotate == 'bifactorT' || fout$rotate == 'bifactorQ') isBifactor <- TRUE
  
  if (itemSort == TRUE) {
    ## Save all output from item sorting function
    itemSorting <- faSort(fmat     = F_view,
                          phi      = Phi_view,
                          reflect  = FALSE,
                          BiFactor = isBifactor)
    
    ## Determine the order in which the items are sorted
    itemOrder <- itemSorting$sortOrder
    
    ## Begin ordering elements
    F_view  <- F_view[itemOrder, ]
    fout$loadings <- F_view
    fout$h2 <- as.matrix(fout$h2[itemOrder,], drop=FALSE)
    fout$R  <- fout$R[itemOrder, itemOrder]
    if (!is.null(fout$Call$bootstrapSE) && fout$Call$bootstrapSE == TRUE) {
      fout$loadingsSE      <- fout$loadingsSE[itemOrder, ]
      fout$loadingsCIlower <- fout$loadingsCIlower[itemOrder, ]
      fout$loadingsCIupper <- fout$loadingsCIupper[itemOrder, ]
      
    } # END if (!is.null(fout$Call$bootstrapSE) && fout$Call$bootstrapSE == TRUE)
    
    
    ## Sort the item order of each solution
    ## Add 1 for the UnSpun solution
    for (iSort in seq_len(NumberLocalSolutions + 1)) {
      
      ## Order the items of each loading matrix
      ## Item order found from solution identified by "Set" argument
      LocalSolutionsList[[iSort]]$loadings <- 
        LocalSolutionsList[[iSort]]$loadings[itemOrder, ]
      
    } # END for (iSort in seq_len(NumberLocalSolutions + 1))
    
  } # END if (itemSort == TRUE)
  
  # ----Compute Structure Matrix
  FS <- F_view %*% Phi_view
  
  # ---- Find set rotation complexity values
  SetComplexityValues <- 
    as.numeric(lapply(LocalSolutionsList, "[[", "RotationComplexityValue"))
  
  
  # ----Find convergence status 
  SetConvergeStatus <- lapply(LocalSolutionsList, "[[", "RotationConverged")
  
  
  
  # ---- Compute rmsd with desired set solution ----
  rmsdVec <- rep(99, NumberLocalSolutions)
  for (iSet in 1:(NumberLocalSolutions + 1)) {
    # align local solution with user chosen view solution
    ## CG NOTE: all solutions are sorted (row-wise) by this point
    F_Set_Aligned <- faAlign(F1 = F_view, 
                             F2 = LocalSolutionsList[[iSet]]$loadings)$F2
    
    rmsdVec[iSet] <-  rmsd(A = F_view, 
                           B = F_Set_Aligned,
                           Symmetric = FALSE,
                           IncludeDiag = TRUE)
  } # END for (iSet in 1:(NumberLocalSolutions + 1)) 
  
  ##-------ENED rmsd Vec
  
  
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
        
 
        ## For each solution within a set, order the items
        for (kSort in 1:NumCasesInSolSets[iSet]) {
          fout$localSolutionSets[[iSet]][[kSort]]$loadings <- 
            fout$localSolutionSets[[iSet]][[kSort]]$loadings[itemOrder, ]
        } # END for (kSort in 1:NumCasesInSolSets[iSet]) 
        
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
  
  
  
  
  
  
  #----Compute Model Fit stats----  
  
  RMSD <- NULL
  RMSAD <- NULL
  
  if( onlyLoadings == FALSE ){
    
    Rhat <- F_view %*% Phi_view %*% t(F_view)
    RMSD <- rmsd(A           = fout$R, 
                 B           = Rhat,
                 Symmetric   = TRUE,
                 IncludeDiag = FALSE)
    
    RMSAD <- rmsd(A           = abs(fout$R), 
                  B           = abs(Rhat),
                  Symmetric   = TRUE,
                  IncludeDiag = FALSE)
  } ## END  if( onlyLoadings ==FALSE )
  
  
  
  ## ---- Print Output ----
  #----____PrintLevel == 0----
  if(PrintLevel == 0){  # No Output printed
    
    #Compute FI for chosen solution
    facIndeterminacy <- sqrt(diag(t(FS) %*% solve(fout$R) %*% FS))
    if(max(facIndeterminacy) > 1){
       facIndeterminacy <- rep(NA, ncol(fout$loadings))
    }
    
    ##  RETURN VALUES
    invisible(return(
      list(loadings             = F_view,
           Phi                  = Phi_view,
           FS                   = FS,
           Set                  = Set,
           h2                   = fout$h2,
           facIndeterminacy     = facIndeterminacy,
           SetComplexityValues  = SetComplexityValues,
           HP_counts            = HP_counts,
           MaxWithinSetRMSD     = maxRmsdWithFi,
           RMSD                 = RMSD,
           RMSAD                = RMSAD,
           NumberLocalSolutions = NumberLocalSolutions,
           LocalSolutions       = LocalSolutionsList,
           rotate               = fout$rotate)))
  } #END if(PrintLevel == 0)
  
  # Print faMain function call
  cat("\nCall:\n")
  print(fout$Call)
  cat("\n")
  
  # Rotation information:
  sectionHeading("Rotation Diagnostics")
  
  #----Print Rotation_Fit_Table----
  cat("\nRotation Method:", fout$rotate)
  cat("\nEpsilon value for rotation convergence and solution set clustering:", fout$rotateControl$epsilon) 
  ## CG Added (7 Sep 2019)
  cat("\nNumber of Rotation Random Starts:", fout$rotateControl$numberStarts)
  cat("\nNumber of Solution Sets:", NumberLocalSolutions)
  cat("\nSolution Characteristics:\n")
  cat("\n")
  print( Rotation_Fit_Table,
         right = TRUE)
  
  
  # ----____Print factor output----  
  cat(paste0("\nFactor Output for Set " , SetArgument, ":\n"))  
  cat("Factor Pattern Matrix\n\n")
  
  # Append communalities
  F_Prnt <- cbind(F_view, "h2" = fout$h2[,1]) 
  print( round(F_Prnt, digits))
  
  cat(paste0("\n\nFactor Correlation Matrix\n"))
  cat("\n")
  print( round( Phi_view, digits))
  
  #----____Print Structure Matrix----
  # PrintLevel == 2:  Sturcture Matrix
  if(PrintLevel == 2){
    cat("\n\nFactor Structure Matrix\n")
    print( round(FS, digits))
  }
  
  #----____Print Bootstrap Output----
  #Bootstrap Output computed for Global minimum only
  
  if( !is.null(fout$loadingsSE) && Set != 1  ){
    sectionHeading("Bootstrap Output")
    cat("\nBootstrap output only available for solutions in Set 1.\n")
  }   
  
  ##IF THERE ARE BOOTSTRAP RESULTS
  if( !is.null(fout$loadingsSE) && Set == 1  ){
    sectionHeading("Bootstrap Output")
    
    if(NumberLocalSolutions > 1){
      cat("\n*** WARNING: Bootstrap Results Are Inaccurate Due to Rotation Local Minima ***\n")
    }
    
       cat("\nFactor Pattern Bootstrap Standard Errors\n")
       print( round(fout$loadingsSE, digits))
       
    
    if(NumberLocalSolutions > 1){
      cat("\n*** WARNING: Bootstrap Results Are Inaccurate Due to Rotation Local Minima ***\n")
    }
    
    
    cat("\nFactor Pattern Bootstrap", paste0(CILevel,"% CIs: Upper Bounds\n"))
    print( round(fout$loadingsCIupper, digits))
    
    if(NumberLocalSolutions > 1){
      cat("\n*** WARNING: Bootstrap Results Are Inaccurate Due to Rotation Local Minima ***\n")
    }
    
    cat("\nFactor Pattern Bootstrap", paste0(CILevel,"% CIs: Lower Bounds\n"))
    print( round(fout$loadingsCIlower, digits))
    
    if(NumberLocalSolutions > 1){
      cat("\n*** WARNING: Bootstrap Results Are Inaccurate Due to Rotation Local Minima ***\n")
    }
    
      
  
    # If orthogonal model then do not print Phi BS results
    if( all.equal(Phi_view, diag(nrow(Phi_view)), check.attributes = FALSE) != TRUE) {
      cat("\nFactor Correlations Bootstrap SEs\n")
      print( round(fout$PhiSE, digits))

      
      if(NumberLocalSolutions > 1){
        cat("\n*** WARNING: Bootstrap Results Are Inaccurate Due to Rotation Local Minima ***\n")
      }
      
      cat("\nFactor Correlations Bootstrap", paste0(CILevel,"% CIs: Upper Bounds\n"))
      print( round(fout$PhiCIupper, digits))
      
      if(NumberLocalSolutions > 1){
        cat("\n*** WARNING: Bootstrap Results Are Inaccurate Due to Rotation Local Minima ***\n")
      }
      
      cat("\nFactor Correlations Bootstrap", paste0(CILevel,"% CIs: Lower Bounds\n"))
      print( round(fout$PhiCIlower, digits))
    }## END if(!all.equal( Phi_view,
         
    cat("\nBootstrap Communality Estimates\n")
         print(round(fout$h2, digits) )
    if(NumberLocalSolutions > 1){
      cat("\n*** WARNING: Bootstrap Results Are Inaccurate Due to Rotation Local Minima ***\n")
    }

  }## END if( !is.null(fout$loadingsSE) && Set == 1  )
  
  #----____Print Factor Indeterminancy----  
  if(!onlyLoadings){
    sectionHeading("Factor Indeterminacy")
    
    Lmin <- eigen(fout$R)$values[nrow(FS)]
    if(Lmin > 1E-6){
   
        #BOOTSTRAP RESULTS AVAILABLE
      if( !is.null(fout$loadingsSE) && Set == 1  ){ 
          cat("\nBootstrap Factor Determinancy Values\n")
          facIndeterminacy <- fout$facIndeterminacy
          print( round(fout$facIndeterminacy, digits) )
          cat("\n")
       } else{ ## IF NO BOOTSTRAP RESULTS AVAILABLE
          # These are the correlations between the factors and factor scores
           F_view_facIndeterminacy <- sqrt(diag(t(FS) %*% solve(fout$R) %*% FS))
         
             if(max(F_view_facIndeterminacy) > 1){
                cat("\nObserved R computationally singular. Unable to calculate FI values.\n")
                facIndeterminacy <- rep(NA, ncol(fout$loadings))
             }else{
               cat("\n")
                print( F_view_facIndeterminacy, digits)
                facIndeterminacy <- F_view_facIndeterminacy
             }# END else
           
       }# END if NO BOOTSTRAP RESULTS AVAILABLE
    }# END if(Lmin > 1E-6)
    
    if(Lmin < 1E-6){
      cat("\nObserved R computationally singular. Unable to calculate FI values.\n")
      facIndeterminacy <- rep(NA, ncol(fout$loadings))
    }
    
    #----____Print Model Fit
    #Model Fit:
    sectionHeading("Model Fit")
    
    
    if(max(fout$h2[,1]) > 1){
      cat("\nWARNING: Heywood case(s) encountered\n")
    }
    
    cat("\nRMSD  (R, Rhat)   = ", round(RMSD, digits = digits +1))
    cat("\nRMSAD (|R|,|Rhat|) = ", round(RMSAD, digits + 1))
    cat("\n\n")
    
  } ##END if( !onlyLoadings)  
  # END Print output 
  
  # Cannot compute Fac determinancy values without R matrix
  # Also ccanot compute Rhat fit measures
  if(onlyLoadings == TRUE){
      facIndeterminacy <- rep(NA, ncol(fout$loadings))
  }#END if(onlyLoadings == TRUE)
  
  
  ##  RETURN VALUES
  invisible(list(loadings = F_view,
                 Phi = Phi_view,
                 FS = FS,
                 Set = SetArgument,
                 h2 = fout$h2,
                 facIndeterminacy = facIndeterminacy,
                 SetComplexityValues = SetComplexityValues,
                 HP_counts = HP_counts,
                 MaxWithinSetRMSD = maxRmsdWithFi,
                 RMSD = RMSD,
                 RMSAD = RMSAD,
                 NumberLocalSolutions = NumberLocalSolutions,
                 LocalSolutions = LocalSolutionsList,
                 itemOrder = itemOrder))
  
  
} ## END summary.faMain 

