#' Multi-Trait Multi-Method correlation matrix reported by Malmi, Underwood, and Carroll (1979).
#' 
#' The original study assessed six variables across three separate 
#' assessment methods. Note that only the last method included six variables whereas 
#' the other two methods included three variables. 
#' 
#' The sample size is \emph{n} = 97.
#'  
#' The following variables were assessed (abbreviations in parentheses):
#' \strong{Variables}:
#' \enumerate{
#'   \item Words (Words)
#'   \item Triads (Triads)
#'   \item Sentences (Sentences)
#'   \item 12 stimuli with 2 responses each (12s.2r)
#'   \item 4 stimuli with 6 responses each (4s.6r)
#'   \item 2 stimuli with 12 responses each (2s.12r)
#' }
#' 
#' The above variables were assessed from the following three assessment methods 
#' (abbreviations in parentheses): 
#' \strong{Test Structure}:
#' \itemize{
#' 
#'   \item \strong{Free Recall} (FR)
#'   \itemize{
#'     \item Words
#'     \item Triads
#'     \item Sentences
#'   }
#'   
#'   \item \strong{Serial List} (SL)
#'   \itemize{
#'     \item Words
#'     \item Triads
#'     \item Sentences
#'   }
#'   
#'   \item \strong{Paired Association} (PA)
#'   \itemize{
#'     \item Words
#'     \item Triads
#'     \item Sentences
#'     \item 12 stimuli with 4 responses
#'     \item 4 stimuli with 6 responses
#'     \item 2 stimuli with 12 responses
#'   }
#'   
#' }
#'
#' @source Malmi, R. A., Underwood, 3. J. & Carroll, J. B. The interrelationships among
#' some associative learning tasks. \emph{Bulletin of the Psychrmomic Society, 13}(3), 121-123.
#' https://doi.org/10.3758/BF03335032
#'
#' @usage data(Malmi79)
#'
#' @docType data
#' @keywords Multiple battery
#' @name Malmi79
#' @format A 12 by 12 correlation matrix with dimension names
#' 
#' @examples
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
#' summary(Example3Output)                        
"Malmi79"