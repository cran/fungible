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
#' @param salient factor markers with loadings >= abs(salient) will be saved in
#' the markers list. Note that a variable can be a marker of more than one
#' factor.
#' @param reflect (logical) if reflect = TRUE then the factors will be
#' reflected such that salient loadings are mostly positive.
#' @return \item{loadings}{sorted factor loadings matrix.} \item{phi}{reflected
#' factor correlation matrix when phi is given as an argument.}
#' \item{markers}{A list of factor specific markers with loadings >=
#' abs(salient). Markers are sorted by the absolute value of the salient factor
#' loadings.} \item{sortOrder}{sorted row numbers.} \item{SEmat}{The SEmat is a
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
#' 
faSort<-function(fmat, phi = NULL, salient = .25, reflect = TRUE){
#------------------------------------------------------#
# faSort: Niels G Waller 
#
# Vers. February 26, 2018
#
# faSort takes an unsorted factor pattern or structure matrix
# and returns a sorted matrix with (possibly)
# reflected columns. Sorting is done so that variables that 
# load on a common factor are grouped together for ease of 
# interpretation.  
#

#
# Arguments:
#
#  fmat:   factor loadings (pattern or structure) matrix
#
#  phi:    factor correlation matrix. If reflect = TRUE
#          then phi will be corrected to match the new factor 
#          orientations.
#
#  salient:  |loadings| >= salient will be saved in markers.
#
#  reflect: (logical) if reflect = TRUE then the factors
#          will be reflected such that salient loadings
#          are mostly positive.
#  
#  Value:
#
#  loadings: sorted factor loadings matrix. 
#
#  phi:  possible reflected factor correlation matrix.  
#
#  markers: a list of salient markers (sorted) for each factor.
#
#  sortOrder: sorted row numbers.
#
#  SEmat: The so-called Start-End matrix lists the first 
#         (start) and last (end) row for each factor in 
#         the sorted pattern matrix.
#------------------------------------------------------#

  Nfac <- ncol(fmat)
  rowNumbers <- 1:nrow(fmat)
  
  # backup of unsorted loadings
  Floadings <- fmat
  
  # append row numbers to unsorted loadings matrix
  fmat <- cbind(fmat, rowNumbers)
  
  # locate initial factor assignments 
  # (max abs loadings for each row)
  FacAssignments <- apply(abs(fmat[,1:Nfac]), 1, which.max)
  
  # Calculate number of salient loadings on each factor
  NumberSalient <- as.vector(table( FacAssignments ))
  
  # varsOnFacs is a vector that lists the variables that load
  # on each factor
  varsOnFacs <- 999
  for(i in 1:Nfac){
    varsOnFacs <- c(varsOnFacs, rowNumbers[FacAssignments == i])
  }
  varsOnFacs <- varsOnFacs[-1]

  # Sort loadings so that loadings
  # on a common factor are grouped together
  fmat <- fmat[varsOnFacs, ]
  
  # Generate an Nfac by 2 matrix, se, that
  # lists the (s)tart and (e)nd row numbers
  # for each factor
  Frows <- c(0, cumsum(NumberSalient))
  # start End matrix
  se <- matrix(0, nrow=Nfac, ncol = 2)
  for(i in 1:Nfac){
    se[i, ] <- (c(Frows[i] + 1, Frows[i+1]))
  }
  
  # Sort loadings by abs value within factors
  fmatAbs <- abs(fmat)
  for(i in 1:Nfac){
    r_order <- Frows[i]+ sort.list(fmatAbs[(se[i,1]:se[i,2]), i], 
                                   decreasing=TRUE)
    fmat[(se[i,1]:se[i,2]), ] <- fmat[r_order, ]
  }
  
  # sortOrder gives the sort order of the 
  # original variables
  sortOrder <- fmat[, Nfac + 1]
  
  fmat <- fmat[, -(Nfac+1)]
  
  # find factor markers (variables that load 
  # above a user-defined salient threhold) 
  # Note that due to cross loadings, a variable
  # can be a marker of more than one factor
  markers <- as.list(rep(9, Nfac))
  for(i in 1:Nfac){
    markers[[i]] <- sortOrder[ abs(fmat[,i]) >= salient ] 
  }
  
  # If reflect = TRUE then reflect factors s.t.
  # salient loadings are positive
  if(reflect == TRUE){
   Dsgn <- diag(sign(colSums(fmat^3)))
   fmat <- fmat %*% Dsgn
    if(!is.null(phi)){
      phi <- Dsgn %*% phi %*% Dsgn
    }
  }
  
  # # Sort markers by |factor loading|
   for(i in 1:Nfac){
     markers[[i]] <- markers[[i]][sort.list(abs(Floadings[markers[[i]],i]), 
                                            decreasing = TRUE)]
   }

 rownames(se) <- paste0("f", 1:Nfac)  
 colnames(se) <-c("Start", "End")
 
 list(loadings= fmat, 
      phi = phi,
      markers = markers,
      sortOrder = sortOrder,
      SEmat =se)
} ## END faSort function
      
###~~~~~~~~~~~~EXAMPLE

# 
# fmat <- matrix(c(-.5,    -.5,   0,    0,
#                  .5,     .5,    0,    0,
#                  .4,      0,    0,    0,
#                  .3,      0,    0,    0,
#                  0,     .5,    0,    0,
#                  0,     .5,    0,    0,
#                  0,      0,   .5,    0,
#                  0,      0,   .35,    0,
#                  0,      0,   .5,    0,
#                  0,      0,    0,   .4,
#                  0,      0,    0,   -.5,
#                  0,      0,    0,    .5), 12, 4, byrow=TRUE)
# 
# 
# #scramble rows
# fmat <- fmat[sample(1:nrow(fmat)), ]
# 
# out <- faSort(fmat, phi = NULL, salient = .25, reflect = TRUE)
# 
# print(out)
