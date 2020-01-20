#' Multi-Trait Multi-Method correlation matrix reported by Thurstone and Thurstone (1941).
#' 
#' The original study assessed a total of 63 variables.
#'  However, we report the 9 variables, across 2 tests,
#'  used to reproduce the multiple battery factor analyses of Browne (1979). 
#' 
#' The sample size is \emph{n} = 710.
#'  
#' The following variables were assessed (abbreviations in parentheses):
#' \strong{Variables}:
#' \itemize{
#'   \item \strong{Test #1} (X)
#'   \itemize{
#'     \item Prefixes (Prefix)
#'     \item Suffixes (Suffix)
#'     \item Sentences (Sentences)
#'     \item Chicago Reading Test: Vocabulary (Vocab)
#'     \item Chicago Reading Test: Sentences (Sentence)
#'   }
#'   \item \strong{Test #2} (Y)
#'   \itemize{
#'     \item First and Last Letters (FLLetters)
#'     \item First Letters (Letters)
#'     \item Four-Letter Words (Words)
#'     \item Completion (Completion)
#'     \item Same and Opposite (SameOpposite)
#'   }
#' }
#' 
#' @source Thurstone, L. L. and Thurstone, T. G. (1941). Factorial studies of 
#' intelligence. \emph{Psychometric Monographs, 2}. Chicago: University Chicago Press. 
#'
#' @usage data(Thurstone41)
#'
#' @docType data
#' @keywords Multiple battery
#' @name Thurstone41
#' @format A 9 by 9 correlation matrix with dimension names
#' 
#' @examples
#' ## Load Thurstone & Thurstone's data used by Browne (1979)
#' data(Thurstone41)
#' Example1Output <-  faMB(R             = Thurstone41, 
#'                         n             = 710,
#'                         NB            = 2, 
#'                         NVB           = c(4,5), 
#'                         numFactors    = 2,
#'                         rotate        = "oblimin",
#'                         rotateControl = list(standardize = "Kaiser"))
#' 
#' summary(Example1Output, PrintLevel = 2)                         
#' 

"Thurstone41"