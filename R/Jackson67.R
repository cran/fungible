#' Multi-Trait Multi-Method correlation matrix reported by Jackson and Singer (1967)
#' 
#' The original study assessed four personality traits (i.e., femininity, 
#' anxiety, somatic complaints, and socially-deviant attitudes) from five 
#' judgemental perspectives (i.e., ratings about (a) desirability in self, 
#' (b) desirability in others, (c) what others find desirable, (d) frequency,
#'  and (e) harmfulness). The harmfulness variable was reverse coded. 
#'  
#'  The sample size is \emph{n} = 480.
#'  
#' The following four variables were assessed (abbreviations in parentheses):
#' \strong{Variables}:
#' \enumerate{
#'   \item Femininity (Fem)
#'   \item Anxiety (Anx)
#'   \item Somatic Complaints (SomatComplaint)
#'   \item Socially-Deviant Attitudes (SDAttitude)
#' }
#' 
#' The above variables were assessed from the following methodological judgement
#' perspectives (abbreviations in parentheses): 
#' \strong{Test Structure}:
#' \itemize{
#'   \item Desirability in the Self (DiS)
#'   \item Desirability in Others (DiO)
#'   \item What Others Find Desirable (WOFD)
#'   \item Frequency (Freq)
#'   \item Harmfulness (Harm)
#' }
#'
#' @source Jackson, D. N., & Singer, J. E. (1967). Judgments, items, and 
#' personality. \emph{Journal of Experimental Research in Personality, 2}(1), 70-79.
#'
#' @usage data(Jackson67)
#'
#' @docType data
#' @keywords Multiple battery
#' @name Jackson67
#' @format A 20 by 20 correlation matrix with dimension names
#' 
#' @examples 
#' ## Load Jackson and Singer's dataset
#' data(Jackson67)
#' 
#' 
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
#' summary(Example2Output)                         

"Jackson67"