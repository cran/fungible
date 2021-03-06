#' Adjective Checklist Data. 
#'
#' Adjective checklist data from the California Twin Registry.
#'
#'
#' @docType data
#' 
#' @usage data(ACL)
#' 
#' @format Adjective Checklist data from the California 
#' Twin Registry (see Waller, Bouchard, Lykken, Tellegen, A., & Blacker,  
#'  1993). 
#'   \strong{ACL} variables:
#'   \enumerate{
#'   \item id
#'   \item sex
#'   \item age
#'   \item items 1 ... 300 
#'   }
  
#' 
#' @details 
#' This is a de-identified subset of the ACL data from the 
#' California Twin Registry (data collected by Waller in the 1990s).  This 
#' data set of 257 cases includes complete (i.e., no missing data) ACL 
#' item responses from a random member of each twin pair.  The item 
#' response vectors are independent.  
#' 
#' @references 
#' Gough, H. G. & Heilbrun, A. B. (1980). The Adjective Checklist 
#' Manual: 1980 Edition. Consulting Psychologists Press.
#' 
#' Waller, N. G., Bouchard, T. J., Lykken, D. T., Tellegen, A., and Blacker, D. 
#'  (1993).  Creativity, heritability, familiarity: 
#'  Which word does not belong?.  Psychological Inquiry, 
#'  4(3), 235--237.  

#' @keywords datasets
#' 
#' 
#' @examples 
#'\dontrun{
#'
#'  data(ACL)
#'  
#' # Factor analyze a random subset of ACL items
#' # for illustrative purposes
#' 
#' set.seed(1)
#' RandomItems <- sample(1:300, 
#'                       50, 
#'                       replace = FALSE)
#' 
#' ACL50 <- ACL[, RandomItems + 3]
#' 
#' tetR_ACL50 <- tetcor(x = ACL50)$r
#' 
#' fout <- faMain(R     = tetR_ACL50,
#'                numFactors    = 5,
#'                facMethod     = "fals",
#'                rotate        = "oblimin",
#'                bootstrapSE   = FALSE,
#'         rotateControl = list(
#'                numberStarts = 100,  
#'                standardize  = "none"),
#'                Seed = 123)
#'
#' summary(fout, itemSort = TRUE)  
#' }   
"ACL"

