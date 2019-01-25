#' Factor Extraction (faX) Routines
#'
#' This function can be used to extract an unrotated factor structure matrix using the following algorithms: (a) unweighted least squares ("fals"); (b) maximum likelihood ("faml"); and (c) iterated principal axis factoring ("fapa").
#'
#' @param R (Matrix) A correlation matrix used for factor extraction.
#' @param numFactors (Numeric) The number of factors to extract for subsequent rotation.
#' @param facMethod (Character) The method used for factor extraction. The supported options are "fals" for unweighted least squares, "faml" for maximum likelihood, and "fapa" for iterated principal axis factoring. The default method is "fals".
#' \itemize{
#'   \item \strong{"fals"}: Factors are extracted using the unweighted least squares estimation procedure using the \code{\link[fungible]{fals}} function.
#'   \item \strong{"faml"}: Factors are extracted using the maximum likelihood estimation procedure using the \code{\link[stats]{factanal}} function.
#'   \item \strong{"fapa"}: Factors are extracted using the iterated principal axis factoring estimation procedure using the \code{\link{fapa}} function.
#' }
#' @param faControl (List) A list of optional parameters passed to the factor extraction routines.
#' \itemize{
#'   \item \strong{treatHeywood}: (Logical) In fals, if treatHeywood is true, a penalized least squares function is used to bound the communality estimates below 1.0. The default is TRUE.
#'   \item \strong{nStart}: (Numeric) In faml, determine the number of starting values to try. The default is 10 start values.
#'   \item \strong{maxCommunality}: (Numeric) In faml, set the maximum communality value for the estimated solution. The default maximum is .995.
#'   \item \strong{epsilon}: (Numeric) In fapa, the numeric threshold designating when the algorithm has converged. The default value is 1e-4.
#'   \item \strong{communality}: (Character) In fapa, the routine requires an initial communality estimate. Select how communalities are initially estimated. The default is squared multiple correlation ("SMC").
#'   \itemize{
#'     \item \strong{"SMC"}: Initial communalities are estimated by taking the squared multiple correlations of each indicator after regressing the indicator on the remaining variables. The following equation is employed to find the squared multiple correlation: \eqn{1 - 1 / diag(R^-1)}.
#'     \item \strong{"maxRsqr"}: Initial communalities equal the largest squared correlation in each column of the correlation matrix.
#'     \item \strong{"unity"}: Initial communalities equal 1.0 for all variables.
#'   }
#'   \item \strong{maxITR}: (Numeric) In fapa, the maximum number of iterations to reach convergence. The default is 15,000
#' }
#'
#' @param digits (Numeric) The number of digits to round all output to.
#'
#' @details
#' \itemize{
#'   \item \strong{Initial communality estimate}: According to Widaman and Herringer (1985), the initial communality estimate does not have much bearing on the resulting solution \emph{when the a stringent convergence criterion is used}. In their analyses, a convergence criterion of .001 (i.e., slightly less stringent than the default of 1e-4) is sufficiently stringent to produce virtually identical communality estimates irrespective of the initial estimate used. It should be noted that all four methods for estimating the initial communality in Widaman and Herringer (1985) are the exact same used in this function. Based on their findings, it is not recommended to use a convergence criterion lower than 1e-3.
#' }
#'
#' @return This function returns a list of output relating to the extracted factor loadings.
#' \itemize{
#'   \item \strong{loadings}: (Matrix) An unrotated factor structure matrix.
#'   \item \strong{h2}: (Vector) A vector containing the item communality estimates.
#'   \item \strong{uniqueness}: (Vector) A vector of the item uniqueness estimates (1 - h2).
#'   \item \strong{Heywood}: (Logical) Whether a Heywood case (i.e., a communality value > 1.0) was detected.
#' }
#'
#' @author
#' \itemize{
#'   \item Casey Giordano (Giord023@umn.edu)
#'   \item Niels G. Waller (nwaller@umn.edu)
#' }
#'
#' @references Widaman, K. F., & Herringer, L. G. (1985). Iterative least squares estimates of communality: Initial estimate need not affect stabilized value. \emph{Psychometrika, 50}(4), 469-477.
#'
#' @examples
#' ## Generate an example factor structure matrix
#' lambda <- matrix(c(.62, .00, .00,
#'                    .54, .00, .00,
#'                    .41, .00, .00,
#'                    .00, .31, .00,
#'                    .00, .58, .00,
#'                    .00, .62, .00,
#'                    .00, .00, .38,
#'                    .00, .00, .43,
#'                    .00, .00, .37),
#'                  nrow = 9, ncol = 3, byrow = TRUE)
#'
#' ## Find the model implied correlation matrix
#' R <- lambda %*% t(lambda)
#' diag(R) <- 1
#'
#' ## Extract (principal axis) factors using the factExtract function
#' Out1 <- faX(R          = R,
#'             numFactors = 3,
#'             facMethod  = "fapa",
#'             faControl  = list(communality = "maxRsqr",
#'                               epsilon     = 1e-4))
#'
#' ## Extract (least squares) factors using the factExtract function
#' Out2 <- faX(R          = R,
#'             numFactors = 3,
#'             facMethod  = "fals",
#'             faControl  = list(treatHeywood = TRUE))
#'
#' @export

faX <- function(R,
                numFactors = NULL,
                facMethod  = "fals",
                faControl  = NULL,
                digits     = NULL) {

  ## ~~~~~~~~~~~~~~~~~~ ##
  #### Error Checking ####
  ## ~~~~~~~~~~~~~~~~~~ ##

  ## Check R

  ## Symmetrical?
  if ( !isSymmetric(R) ) {
    stop("'R' is not a symmetric correlation matrix.")
  } # END if ( !isSymmetric(R) ) {

  ## Unit diagonal?
  if ( all.equal(diag(R), rep(1, ncol(R)), check.attributes = FALSE) != TRUE ) {
    stop("'R' must have a unit diagonal")
  } # END if ( any(diag(R) != 1) )

  ## Positive definite
  eigs <- eigen(R)$values
  if ( min(eigs) <= 0 ) {
    stop("'R' is not a positive definite correlation matrix. ")
  } # END if ( min(eigs) <= 0 )

  ## Check facMethod

  ## Specified correctly?
  facMethodOptions <- c("fals", "faml", "fapa")
  if ( (facMethod %in% facMethodOptions) == FALSE ) {
    stop("The method of factor extraction is incorrectly specified. Select either 'fals', 'faml', or 'fapa'.")
  } # END if ( (facMethod %in% facMethodOptions) == FALSE )

  ## Check numFactors

  ## Specified?
  if ( is.null(numFactors) ) {
    stop("The 'numFactors' argument must be specified.")
  } # END if ( is.null(numFactors) )

  ## Check faControl

  ## Assert the default options
  cnFA <- list(treatHeywood   = TRUE,
               nStart         = 10,
               maxCommunality = .995,
               epsilon        = 1e-4,
               communality    = "SMC",
               maxITR         = 15000)

  ## Correctly specified?
  if ( !is.null(faControl) ) {

    ## Total number of correct names
    cnFALength <- length( names(cnFA) )

    ## Total number of all names supplied
    allcnFALength <- length( unique( c( names(faControl), names(cnFA) ) ) )

    ## If lengths differ, something is mis-specified
    if (cnFALength != allcnFALength) {

      ## Find which are incorrect
      incorrectFAArgs <- which( (names(faControl) %in% names(cnFA)) == FALSE)

      stop(paste("The following arguments are not valid inputs for the list of factor analysis control arguments:", paste(c(names(faControl)[incorrectFAArgs]), collapse = ", "), collapse = ": " ) )

    } # END if (cnFALength != allcnFALength)

  } # END if ( !is.null(faControl) )

  ## Check digits

  ## If not specified, give it a value
  if ( is.null(digits) ) {
    digits <- 100
  } # END if ( is.null(digits) )

  ## ~~~~~~~~~~~~~~~~~~ ##
  #### Begin Function ####
  ## ~~~~~~~~~~~~~~~~~~ ##

  Out <- switch(facMethod,
                "fals" = {
                  fals(R            = R,
                       nfactors     = numFactors,
                       TreatHeywood = cnFA$treatHeywood)
                },
                "faml" = {
                  factanal(covmat   = R,
                           factors  = numFactors,
                           rotation = "none",
                           control  = list(nstart = cnFA$nStart,
                                           lower  = 1 - cnFA$maxCommunality))
                },
                "fapa" = {
                  fapa(R            = R,
                       numFactors   = numFactors,
                       epsilon      = cnFA$epsilon,
                       communality  = cnFA$communality,
                       maxITR       = cnFA$maxITR)
                })

  ## Compute the estimated communality values
  hsq <- apply(Out$loadings^2, 1, sum)

  ## Compute the uniqueness values
  uniqueness <- (1 - hsq)

  ## Check for Heywood cases
  Heywood <- FALSE
  if ( any(hsq > 1) ) {
    Heywood <- TRUE
    warning("A Heywood case was detected in the extracted factor structure matrix.")

  } # END if ( any(hsq > 1) )


  list(loadings   = round(Out$loadings, digits),
       h2         = round(hsq,          digits),
       uniqueness = round(uniqueness,   digits),
       Heywood    = Heywood)

} # END factExtract















