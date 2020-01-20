#' Print  Method for an Object of Class faMB
#' 
#' @param x (Object of class \strong{faMB}) The returned object 
#' from a call to \strong{faMB}.
#' @param digits (Integer) Print output with user-specified number of significant digits. Default \code{digits = 2}.   
#' @param Set 
#'   \itemize{
#'     \item{integer} (Integer) Summarize the solution from the specified solution set. 
#'     \item{'UnSpun'} (Character) Summarize the solution from the rotated output that was 
#'      produced by rotating from the unrotated (i.e., unspun) factor orientation.
#'    }
#' @param itemSort (Logical) If TRUE, sort the order of the observed variables to produce
#' a "staircase"-like pattern.  Defaults to \code{itemSort} = FALSE.
#' @param \dots Additional arguments affecting the summary produced.
#' 
#' @family Factor Analysis Routines
#' 
#' @export

print.faMB <- function(x, 
                         ...,
                         digits   = 2,
                         Set      = 1,
                         itemSort = FALSE){
  
  # is the object of the appropriate class
  stopifnot(inherits(x, "faMB"))
  
  
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
  
  # Align selected solution (which may not be the 
  # min rotation solution) with the fout$loadings from faMB
  
  ## Only reflect and sort columns if more than 2 columns are supplied
  if ( ncol(fout$loadings) == 1 ) {
    
    ## If factor (of selected solution) is negatively oriented, flip it 
    Dsgn <- sign( sum(F_view^3) )
    
    ## Orient the lone factor to always be positively-oriented
    F_view[, 1] <- F_view * Dsgn
    
  } # END if ( ncol(fout$loadings) == 1 ) 
  
  if ( ncol(fout$loadings) > 1 ) {
    ## Align selected solution with global minimum (easier comparison)
    F_viewAlignOut <- faAlign(F1          = fout$loadings,
                              F2          = F_view, 
                              Phi2        = Phi_view,
                              MatchMethod = "LS")
    F_view <- F_viewAlignOut$F2
    Phi_view <- F_viewAlignOut$Phi2
  } # END if ( ncol(fout$loadings) > 1 ) 
  
  colnames(F_view) <- paste0("f", 1:ncol(F_view))
  rownames(Phi_view) <- colnames(Phi_view) <- paste0("f", 1:ncol(F_view))
  
  
  RotationConverged <- 
    LocalSolutionsList[[Set]]$RotationConverged
  # ---- END Which solution should be viewed? 
  
  
  ## ---- Sort the items ----
  if (itemSort == TRUE) {
    ## Save all output from item sorting function
    itemSorting <- faSort(fmat     = F_view,
                          phi      = Phi_view,
                          BiFactor = FALSE)
    
    ## Determine the order in which the items are sorted
    itemOrder <- itemSorting$sortOrder
  } # END if (itemSort == TRUE)
  
  
  # ----Compute Structure Matrix
  FS <- F_view %*% Phi_view
  
  # ----____Print factor loadings output---- 
  cat("\nMultiple Battery Factor Analysis")
  cat("\nRotation Method:", fout$rotate)
  if(!is.null(fout$targetMatrix)){
    cat(paste0("\nAligned Factor Output for Set " , SetArgument, " of ", 
               NumberLocalSolutions, ":\n\n"))
  } else { 
    cat(paste0("\nFactor Output for Set " , SetArgument, " of ", 
                   NumberLocalSolutions, ":\n\n"))
  }
  

  ## Print a label for the matrix
  cat("Factor Pattern Matrix\n")
  print( round(F_view[itemOrder, , drop = FALSE], digits))
  
  # ----____Print factor correlation output---- 
  cat(paste0("\n\nFactor Correlation Matrix\n"))
  print( round(Phi_view, digits))
  
  #----____Print Structure Matrix----
  ## itemOrder is either a sorted vector or 1:numVariables
  cat("\nFactor Structure Matrix\n")
  print( round(FS[itemOrder, , drop = FALSE], digits))
  
  if(NumberLocalSolutions > 1){
     cat("\nTo view a different solution (e.g., from solution Set 2), type: print(object, Set = 2)\n\n")
  }
  
  ##  RETURN VALUES
  invisible(list(loadings = F_view[itemOrder, ], 
                 Phi      = Phi_view,
                 FS       = FS[itemOrder, ], 
                 Set      = SetArgument)) 
  
} ## END print.faMB 

