#' Sort a factor loadings matrix
#' 
#' faSort takes an unsorted factor pattern or structure matrix and returns a
#' sorted matrix with (possibly) reflected columns. Sorting is done such that
#' variables that load on a common factor are grouped together for ease of
#' interpretation.
#' 
#' 
#' @param fmat factor loadings (pattern or structure) matrix.
#' @param phi factor correlation matrix. Default = NULL. If reflect = TRUE then
#' phi will be corrected to match the new factor orientations.
#' @param BiFactor (logical) Is the solution a bifactor model?
#' @param salient factor markers with loadings >= abs(salient) will be saved in
#' the markers list. Note that a variable can be a marker of more than one
#' factor.
#' @param reflect (logical) if reflect = TRUE then the factors will be
#' reflected such that salient loadings are mostly positive.
#' @return 
#'    \item{loadings}{sorted factor loadings matrix.} 
#'    \item{phi}{reflected factor correlation matrix when phi is given as an argument.}
#'    \item{markers}{A list of factor specific markers with loadings >=
#' abs(salient). Markers are sorted by the absolute value of the salient factor
#' loadings.} \item{sortOrder}{sorted row numbers.} 
#'    \item{SEmat}{The SEmat is a
#' so-called Start-End matrix that lists the first (start) and last (end) row
#' for each factor in the sorted pattern matrix.}
#' @author Niels Waller
#' @seealso \code{\link{fals}}
#' @keywords Statstics
#' @family Factor Analysis Routines
#' @export
#' @examples
#' 
#' set.seed(123)
#' F <- matrix( c( .5,  0, 
#'                 .6,  0,
#'                  0, .6,
#'                 .6,  0,
#'                  0, .5,
#'                 .7,  0,
#'                  0, .7,
#'                  0, .6), nrow = 8, ncol = 2, byrow=TRUE)
#' 
#' Rex1 <- F %*% t(F); diag(Rex1) <- 1
#' 
#' Items <- c("1. I am often tense.\n",
#'            "2. I feel anxious much of the time.\n",
#'            "3. I am a naturally curious individual.\n",
#'            "4. I have many fears.\n",
#'            "5. I read many books each year.\n",
#'            "6. My hands perspire easily.\n",
#'            "7. I have many interests.\n",
#'            "8. I enjoy learning new words.\n")
#' 
#' exampleOut <- fals(R = Rex1, nfactors = 2)
#' 
#' # Varimax rotation
#' Fload <- varimax(exampleOut$loadings)$loadings[]
#' 
#' # Add some row labels
#' rownames(Fload) <- paste0("V", 1:nrow(Fload))
#' 
#' cat("\nUnsorted fator loadings\n")
#' print(round( Fload, 2) )
#' 
#' # Sort items and reflect factors
#' out1 <- faSort(fmat = Fload, 
#'                salient = .25, 
#'                reflect = TRUE)
#'                
#' FloadSorted <- out1$loadings
#' 
#' cat("\nSorted fator loadings\n")
#' print(round( FloadSorted, 2) )
#' 
#' # Print sorted items
#' cat("\n Items sorted by Factor\n")
#' cat("\n",Items[out1$sortOrder])

faSort <- function(fmat, 
                   phi      = NULL, 
                   BiFactor = FALSE, 
                   salient  = .25, 
                   reflect  = TRUE) {
  
  ## Save the original loadings matrix
  fmatOriginal <- fmat
  
  ## Are we working with a bifactor solution?
  if(BiFactor) fmat <- fmat[, -1]
  
  Nfac <- ncol(fmat)
  Nrow <- nrow(fmat)
  rowNames <- rownames(fmat) 
  rowNumbers <- 1:Nrow
  
  
  ## If single factor model
  if (Nfac == 1) {
    itemOrder <- sort.list(abs(fmat), 
                           decreasing = TRUE)
    fmat <- fmat[itemOrder, 1, drop = FALSE]
    dimnames(fmat) <- list(rowNames[itemOrder], "f1")
    markers = NULL
    sortOrder = itemOrder
    se = NULL
    return(  
      list(loadings= fmat, 
           phi = phi,
           markers = markers,
           sortOrder = sortOrder,
           SEmat =se)
    )   
  } ## END if Nfac == 1
  
  # backup of unsorted loadings
  Floadings <- fmat
  
  ## append row numbers to unsorted loadings matrix
  fmat <- cbind(fmat, rowNumbers)
  
  # # locate initial factor assignments 
  # # (max abs loadings for each row)
  FacAssignments <- apply(abs(fmat[,1:Nfac]), 1, which.max)
  
  
  ## CG EDITS (6 Sept 2019)
  ## Vector to count how many salient indicators per factor, start with zero
  NumberSalient <- rep(0, Nfac)
  names(NumberSalient) <- paste(seq_len(Nfac)) ## Add names
  
  ## CG EDITS (6 Sept 2019)
  ## Count salient loadings per factor
  countSalient <- table(FacAssignments)
  
  ## CG EDITS (6 Sept 2019)
  ## Overwrite placeholders in NumberSalient with salient count
  NumberSalient[names(countSalient)] <- countSalient
  
  # Calculate number of salient loadings on each factor
  ## NOTE: This  blows up if a factor has no salient indicators
  # NumberSalient <- as.vector(table( FacAssignments ))
  
  # varsOnFacs is a vector that lists the variables that load
  # on each factor
  ## CG NOTE: pre-allocation does not work, multiple values per loop
  varsOnFacs <- NA
  for (i in 1:Nfac) {
    varsOnFacs <- c(varsOnFacs, rowNumbers[FacAssignments == i])
  } # END for (i in 1:Nfac)
  
  ## Drop the initial NA value
  varsOnFacs <- varsOnFacs[-1]
  
  # Sort loadings so that loadings
  # on a common factor are grouped together
  fmat <- fmat[varsOnFacs, ]
  
  
  # Generate an Nfac by 2 matrix, se (start/end), that
  # lists the (s)tart and (e)nd row numbers
  # for each factor
  
  ## Find variable number where each factor ends
  ## CG NOTE: start with 0, next loop adds 1 
  Frows <- c(0, cumsum(NumberSalient))
  
  # Start/End matrix
  se <- matrix(0, nrow=Nfac, ncol = 2)
  for (i in 1:Nfac) {
    
    ## CG EDITS (13 sept 19): comments added
    ## 1st element: Factor i starts where factor (i - 1) ends, plus one 
    ## 2nd element: Factor i ends in Frow[i] but we appended a zero to Frow
      ## such that we need to take Frow[i + 1] instead
    se[i, ] <- (c(Frows[i] + 1, Frows[i+1]))
  } # END for (i in 1:Nfac) 
  
  
  ## If no general factor found then sort items 
    ## within factors
  if (se[1, 2] != Nrow ) { ## If all items load on the first factor...
    
    # Sort loadings by abs value within factors
    fmatAbs <- abs(fmat)
    for (i in 1:Nfac) {
      
      ## CG EDITS: If factor has no salient items, do nothing 
      if (NumberSalient[i] == 0) {
        ## If a factor has no salient markers, do not sort it
        ## Do nothing, ignore this factor
        fmat <- fmat
      } else {
        
        ## CG EDITS (13 sept 19): Added comments
        
        ## For factor i, only sort items with salient markers on that factor
        ## Vector is sorted order of items loading on factor i
        r_order <- Frows[i] + sort.list(fmatAbs[(se[i, 1]:se[i, 2]), i], 
                                        decreasing = TRUE)
        ## Updated loadings matrix with new item order
        fmat[(se[i,1]:se[i,2]), ] <- fmat[r_order, ]
      } # END if (NumberSalient[i] == 0) 
      
    } # END for (i in 1:Nfac) 
  } # END  if (se[1,2] != Nrow )
  
  
  
  
  # sortOrder gives the sort order of the 
  # original variables
  ## CG EDITS (13 sept 19): Nfac + 1 is "rowNumbers"
  sortOrder <- fmat[, Nfac + 1]
  
  ## Drop rowNumbers column from fmat, no longer needed
  fmat <- fmat[, -(Nfac+1)]
  
  # find factor markers (variables that load 
  # above a user-defined salient threhold) 
  # Note that due to cross loadings, a variable
  # can be a marker of more than one factor
  markers <- as.list(rep(NA, Nfac))
  for (iMark in 1:Nfac) {
    markers[[iMark]] <- sortOrder[ abs(fmat[, iMark]) >= salient ] 
  } # END for (iMark in 1:Nfac) 
  
  # If reflect = TRUE then reflect factors s.t.
  # salient loadings are positive
  if (reflect == TRUE) {
    ## Determine whether factors are negatively oriented
    Dsgn <- diag(sign(colSums(fmat^3))) ## CG (13 sept 19): Is the cube needed??
    ## If factors negatively oriented, multiply by -1, else multiply by 1
    fmat <- fmat %*% Dsgn
    if (!is.null(phi)) {
      ## If factor is negative, reverse corresponding factor correlations
      phi <- Dsgn %*% phi %*% Dsgn
    } # END if (!is.null(phi)) 
  } # END if (reflect == TRUE) 
  
  # Sort markers by |factor loading|
  for (i in 1:Nfac) {
    markers[[i]] <- markers[[i]][sort.list(abs(Floadings[markers[[i]],i]), 
                                           decreasing = TRUE)]
  } # END for (i in 1:Nfac) 
  
  ## CG EDITS (13 sept 19): Add names to list of factor markers
  names(markers) <- paste("Factor", 1:Nfac, "markers")
  
  ## Add dimension names to the Start/End matrix
  rownames(se) <- paste0("f", 1:Nfac)  
  colnames(se) <-c("Start", "End")
  
  if(!BiFactor){
    fmatReturn <- fmat
  }
  if(BiFactor){
    fmatReturn <- cbind(fmatOriginal[sortOrder, 1], fmat)
  }
  
  ## Return list of output
  list(loadings  = fmatReturn, 
       phi       = phi,
       markers   = markers,
       sortOrder = sortOrder,
       SEmat     = se)
  
} ## END faSort function


