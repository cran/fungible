#' Multi-Trait Multi-Method correlation matrix reported by Boruch, Larkin, Wolins, and MacKinney (1970)
#' 
#' The original study assessed supervisors on seven dimensions (i.e., 7 
#' variables) from two sources (i.e., their least effective and most effective subordinate). 
#' 
#' 
#' The sample size is \emph{n} = 111.
#' 
#' The following variables were assessed:
#' \strong{Variables}:
#' \enumerate{
#'   \item Consideration
#'   \item Structure
#'   \item Satisfaction with the supervisor
#'   \item Job satisfaction
#'   \item General effectiveness
#'   \item Human relations skill
#'   \item Leadership
#' }
#' 
#' The test structure is as follows:
#' \strong{Test Structure}:
#' \itemize{
#'   \item Test One: variables 1 through 7
#'   \item Test Two: variables 8 through 14
#' }
#'
#' @source Boruch, R. F., Larkin, J. D., Wolins, L., and MacKinney, A. C. (1970). 
#' Alternative methods of analysis: Multitrait multimethod data. \emph{Educational 
#' and Psychological Measurement, 30}, 833-853.
#'
#' @usage data(Boruch70)
#'
#' @docType data
#' @keywords Multiple battery
#' @name Boruch70
#' @format A 14 by 14 correlation matrix with dimension names
#' 
#' @examples 
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
#' summary(Example4Output, digits = 3)                                             

"Boruch70"