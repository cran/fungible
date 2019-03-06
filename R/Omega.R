#' Compute Omega hierarchical
#'
#' This function computes McDonald's Omega hierarchical to determine the proportions of variance (for a given test) associated with the latent factors and with the general factor.
#'
#' @param lambda (Matrix) A factor pattern matrix to be analyzed.
#' @param genFac (Scalar, Vector) Which column(s) contains the general factor(s). The default value is the first column.
#' @param digits (Scalar) The number of digits to round all output to.
#'
#' @details
#' \itemize{
#'   \item \strong{Omega Hierarchical}: For a reader-friendly description (with some examples), see the Rodriguez et al., (2016) \emph{Psychological Methods} article. Most of the relevant equations and descriptions are found on page 141.
#' }
#'
#' @return
#' \itemize{
#'   \item \strong{omegaTotal}: (Scalar) The total reliability of the latent, common factors for the given test.
#'   \item \strong{omegaGeneral}: (Scalar) The proportion of total variance that is accounted for by the general factor(s).
#' }
#'
#' @author
#' \itemize{
#'   \item Casey Giordano (Giord023@umn.edu)
#'   \item Niels G. Waller (nwaller@umn.edu)
#'}
#'
#' @references McDonald, R. P. (1999). \emph{Test theory: A unified approach}. Mahwah, NJ:Erlbaum.
#' @references Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor models: Calculating and interpreting statistical indices. \emph{Psychological Methods, 21}(2), 137.
#' @references Zinbarg, R.E., Revelle, W., Yovel, I., & Li. W. (2005). Cronbach's Alpha, Revelle's Beta, McDonald's Omega: Their relations with each and two alternative conceptualizations of reliability. \emph{Psychometrika. }70, 123-133. https://personality-project.org/revelle/publications/zinbarg.revelle.pmet.05.pdf
#'
#' @examples
#' ## Create a bifactor structure
#' bifactor <- matrix(c(.21, .49, .00, .00,
#'                      .12, .28, .00, .00,
#'                      .17, .38, .00, .00,
#'                      .23, .00, .34, .00,
#'                      .34, .00, .52, .00,
#'                      .22, .00, .34, .00,
#'                      .41, .00, .00, .42,
#'                      .46, .00, .00, .47,
#'                      .48, .00, .00, .49),
#'                    nrow = 9, ncol = 4, byrow = TRUE)
#'
#' ## Compute Omega
#' Out1 <- Omega(lambda = bifactor)
#'
#' @export

Omega <- function(lambda,
                  genFac = 1,
                  digits = NULL) {

  ## ~~~~~~~~~~~~~~~~~~ ##
  #### Error Checking ####
  ## ~~~~~~~~~~~~~~~~~~ ##

  ## If digits is not specified, give it arbitrarily large value
  if ( is.null(digits) ) {
    digits <- 100
  } # END if ( is.null(digits) )

  ## Is lambda a factor structure?
  hsq <- apply(lambda^2, 1, sum)

  if ( any(hsq > 1.0) ) {
    stop("The supplied 'lambda' factor structure has a Heywood case (i.e., at least one item communality exceeds 1.0.")
  } # END if ( any(hsq > 1.0) )

  ## ~~~~~~~~~~~~~~ ##
  #### Begin Code ####
  ## ~~~~~~~~~~~~~~ ##

  ## Compute the indicator uniqueness values (1 - communality)
  Uniqueness <- 1 - diag( tcrossprod( lambda ) )

  ## Compute the squared column sums
  SquaredSum <- colSums(lambda)^2

  ## Compute overall omega (genFac + GrpFac in numerator)
  ## Rodriguez, Reise, & Haviland, 2016, Psych. Methods., equation 4, p. 141
  Total <- sum( SquaredSum ) / sum( SquaredSum, Uniqueness)

  ## Compute omega hierarchical: general factor (genFac in numerator)
  ## Rodriguez, Reise, & Haviland, 2016, Psych. Methods., equation 6, p. 141
  General <- sum( SquaredSum[genFac] ) / sum( SquaredSum, Uniqueness)

  ## Return with overall and general factor omega
  list(OmegaTotal   = round(Total, digits),
       OmegaGeneral = round(General, digits))

}# END Omega


