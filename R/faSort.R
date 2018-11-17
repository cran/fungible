

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
