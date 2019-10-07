#' Print  Method for an Object of Class faMain
#' 
#' @param x (Object of class \strong{faMain}) The returned object 
#' from a call to \strong{faMain}.
#' @param digits (Integer) Print output with user-specified number of significant digits. Default \code{digits = 2}.   
#' @param Set 
#'   \itemize{
#'     \item{integer} (Integer) Summarize the solution from the specified solution set. 
#'     \item{'UnSpun'} (Character) Summarize the solution from the rotated output that was 
#'      produced by rotating from the unrotated (i.e., unspun) factor orientation.
#'    }
#' @param itemSort (Logical) If TRUE, sort the order of the observed variables to produce
#' a "staircase"-like pattern. In bifactor models (i.e., bifactorT and bifactorQ) item 
#' sorting is determined by  the magnitudes of the group factor loadings.
#' Defaults to \code{itemSort} = FALSE.
#' @param \dots Additional arguments affecting the summary produced.
#' 
#' @family Factor Analysis Routines
#' 
#' @export

print.faMain <- function(x, 
                         ...,
                         digits   = 2,
                         Set      = 1,
                         itemSort = FALSE){
  
  # is the object of the appropriate class
  stopifnot(inherits(x, "faMain"))
  
  
  fout <- x
  SetArgument <- Set
  
  if (itemSort == FALSE)  itemOrder <- seq_len(nrow(fout$loadings))
  
  ## ---Count number of local solution sets
  NumberLocalSolutions <- fout$numLocalSets
  if(Set == "UnSpun") Set = NumberLocalSolutions + 1
  if(Set > NumberLocalSolutions && 
     SetArgument !="UnSpun") stop("\nOut of bound value for Set\n")
  
  ## ---- Gather local Solutions----
  LocalSolutionsList <- list()
  for(iSol in 1: NumberLocalSolutions){
    #Pick first element in solution set
    LocalSolutionsList[[iSol]] <- c(iSol,
                                    fout$localSolutionSets[[iSol]][[1]])
    names(LocalSolutionsList[[iSol]])<-c("Set", 
                                         names(fout$localSolutionSets[[iSol]][[1]]))
  }
  
  # UnSpun Solution is last element in list
  LocalSolutionsList[[NumberLocalSolutions + 1]] <- c(0, fout$unSpunSolution)
  names(LocalSolutionsList[[NumberLocalSolutions + 1]])<-c("Set", names(fout$unSpunSolution))
  ## END Gather local Solutions
  
  # ---- Which solution should be viewed?-----
  F_view <- LocalSolutionsList[[Set]]$loadings
  Phi_view <- LocalSolutionsList[[Set]]$Phi
  dimnames(F_view) <- dimnames(fout$loadings)
  dimnames(Phi_view) <- dimnames(fout$Phi)
  RotationConverged <- 
    LocalSolutionsList[[Set]]$RotationConverged
  # ---- END Which solution should be viewed? 
  
  isBifactor <- FALSE
  if(fout$rotate == 'bifactorT' || fout$rotate == 'bifactorQ') isBifactor <- TRUE
  
  ## ---- Sort the items ----
  if (itemSort == TRUE) {
    ## Save all output from item sorting function
    itemSorting <- faSort(fmat     = F_view,
                          phi      = Phi_view,
                          BiFactor = isBifactor)
    
    ## Determine the order in which the items are sorted
    itemOrder <- itemSorting$sortOrder
  } # END if (itemSort == TRUE)
  
  # If full targetMatrix supplied then align viewed solution  
  if ( !is.null(fout$targetMatrix) && sum(is.na(fout$targetMatrix)) ==0 ) {
    alignOut <- faAlign(F1   = fout$targetMatrix,
                        F2   = F_view, 
                        Phi2 = Phi_view)
    
    F_view <- alignOut$F2
    Phi_view <- alignOut$Phi
    dimnames(F_view) <- dimnames(fout$loadings)
    dimnames(Phi_view) <- dimnames(fout$Phi)
  } #END  # If full targetMatrix supplied then align viewed solution  
  
  # ----Compute Structure Matrix
  FS <- F_view %*% Phi_view
  
  # ----____Print factor loadings output---- 
  cat("\nRotation Method:", fout$rotate)
  if(!is.null(fout$targetMatrix)){
    cat(paste0("\nAligned Factor Output for Set " , SetArgument, " of ", 
               NumberLocalSolutions, ":\n\n"))
  } else { 
    cat(paste0("\nFactor Output for Set " , SetArgument, " of ", 
                   NumberLocalSolutions, ":\n\n"))
  }
  
  # Append commonalities
  F_Prnt <- cbind(F_view, "h2" = fout$h2[,1])
  
  ## CG Edits (13 sept 19): Removed ifelse statement
    ## itemOrder is either a sorted vector or 1:numVariables
  ## Print a label for the matrix
  cat("Factor Pattern Matrix\n")
  print( round(F_Prnt[itemOrder, ], digits))
  
  # ----____Print factor correlation output---- 
  cat(paste0("\n\nFactor Correlation Matrix\n"))
  print( round(Phi_view, digits))
  
  #----____Print Structure Matrix----
  ## CG Edits (13 sept 19): Removed ifelse statement
    ## itemOrder is either a sorted vector or 1:numVariables
  cat("\nFactor Structure Matrix\n")
  print( round(FS[itemOrder, ], digits))
  
  cat("\nTo view a different solution (e.g., from solution Set 2), type: print(object, Set = 2)\n\n")
  
  ##  RETURN VALUES
  invisible(list(loadings = F_view[itemOrder, ], 
                 Phi      = Phi_view,
                 FS       = FS[itemOrder, ], 
                 Set      = SetArgument,
                 h2       = fout$h2[itemOrder,1])) 
  
} ## END print.faMain 

