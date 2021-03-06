#' Factor Pattern and Factor Correlations for Thurstone's 20 hypothetical box attributes. 
#'
#' Factor Pattern and Factor Correlations for Thurstone's 20 hypothetical box attributes. 
#'
#' @docType data
#' 
#' @usage data(ThurstoneBox20)
#' 
#' @format This is a list containing the \code{Loadings} (original factor pattern) and \code{Phi}
#' matrix (factor correlation matrix) from Thurstone's 20 Box problem (Thurstone, 1940, p. 227). 
#' The original 20-variable Box problem  contains measurements on the following score
#'  functions of box length (x), width (y), and height (z). 
#' \strong{Box20} variables:
#' \enumerate{
#'      \item   x^2
#'      \item   y^2
#'      \item   z^2
#'      \item   xy
#'      \item   xz
#'      \item   yz
#'      \item   sqrt(x^2 + y^2)
#'      \item   sqrt(x^2 + z^2)
#'      \item   sqrt(y^2 + z^2)
#'      \item   2x + 2y
#'      \item   2x + 2z
#'      \item   2y + 2z 
#'      \item   log(x)
#'      \item   log(y)
#'      \item   log(z)
#'      \item   xyz
#'      \item   sqrt(x^2 + y^2 + z^2)
#'      \item   exp(x)
#'      \item   exp(y)
#'      \item   exp(z)
#'    }
#' 
#' @keywords datasets
#' 
#' @details 
#' Two data sets have been described in the literature as Thurstone's Box Data 
#' (or Thurstone's Box Problem). The first consists of 20 measurements on a set of 20 
#' hypothetical boxes (i.e., Thurstone made up the data).  Those data are available 
#' in \strong{Box20}. 
#' 
#' @references 
#' 
#' Thurstone, L. L. (1940). Current issues in factor analysis. Psychological Bulletin, 37(4), 189. 
#' Thurstone, L. L.  (1947).  Multiple factor analysis.  Chicago: University of Chicago Press. 
#' 
#' @seealso \code{\link{AmzBoxes}}, \code{\link{Box20}}, \code{\link{Box26}},
#' \code{\link{GenerateBoxData}} 
#' 
#' @examples 
#' data(ThurstoneBox20)  
#' ThurstoneBox20
"ThurstoneBox20"


