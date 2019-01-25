#' Compute the cosine(s) between either 2 matrices or 2 vectors.
#'
#' This function will compute the cosines (i.e., the angle) between two vectors or matrices. When applied to matrices, it will compare the two matrices one vector (i.e., column) at a time. For instance, the cosine (angle) between factor 1 in matrix A and factor 1 in matrix B.
#'
#' @param A (Matrix, Vector) Either a matrix or vector.
#' @param B (Matrix, Vector) Either a matrix or vector (must be of the same dimensions as A).
#' @param align (Logical) Whether to run a factor alignment before computing the cosine.
#' @param digits (Numeric) The number of digits to round the output to.
#'
#' @details
#' \itemize{
#'   \item \strong{Chance Congruence}: Factor cosines were originally described by Burt (1948) and later popularized by Tucker (1951). Several authors have noted the tendency for two factors to have spuriously large factor cosines. Paunonen (1997) provides a good overview and describes how factor cosines between two vectors of random numbers can appear to be congruent.
#'   \item \strong{Effect Size Benchmarks}: When computing congruence coefficients (cosines) in factor analytic studies, it can be useful to know what constitutes large versus small congruence. Lorenzo-Seva and ten Berge (2006) currently provide the most popular (i.e., most frequently cited) recommended benchmarks for congruence. ``A value in the range .85-.94 means that the two factors compared display \emph{fair} similarity. This result should prevent congruence below .85 from being interpreted as indicative of any factor similarity at all. A value higher than .95 means that the two factors or components compared can be considered equal. That is what we have called a \emph{good} similarity in our study'' (Lorenzo-Seva & ten Berge, 2006, p. 61, emphasis theirs).
#' }
#'
#' @return A vector of cosines will be returned. When comparing two vectors, only one cosine can be computed. When comparing matrices, one cosine is computed per column.
#' \itemize{
#'   \item \strong{cosine}: (Matrix) A matrix of cosines between the two inputs. 
#'   \item \strong{A}: (Matrix) The A input matrix.
#'   \item \strong{B}: (Matrix) The B input matrix. 
#'   \item \strong{align}: (Logical) Whether Matrix B was aligned to A.
#' }
#'
#' @author
#' \itemize{
#'   \item Casey Giordano (Giord023@umn.edu)
#'   \item Niels G. Waller (nwaller@umn.edu)
#'}
#'
#' @references Burt, C. (1948). The factorial study of temperament traits. \emph{British Journal of Psychology, Statistical Section, 1}, 178-203.
#' @references Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tuckers Congruence Coefficient as a meaningful index of factor similarity. \emph{Methodology, 2}(2), 57-64.
#' @references Paunonen, S. V. (1997). On chance and factor congruence following orthogonal Procrustes rotation. \emph{Educational and Psychological Measurement}, 57, 33-59.
#' @references Tucker, L. R. (1951). \emph{A method for synthesis of factor analysis studies} (Personnel Research Section Report No. 984). Washington, DC: Department of the Army.
#'
#' @examples
#' ## Cosine between two vectors
#' A <- rnorm(5)
#' B <- rnorm(5)
#'
#' cosMat(A, B)
#'
#' ## Cosine between the columns of two matrices
#' A <- matrix(rnorm(5 * 5), 5, 5)
#' B <- matrix(rnorm(5 * 5), 5, 5)
#'
#' cosMat(A, B)
#'
#' @export

cosMat <- function(A,
                   B,
                   align  = FALSE,
                   digits = NULL) {



  ## If digits argument is not supplied, give arbitrarily big value
  if ( is.null(digits) ) digits <- options()$digits

  ## Ensure proper class to detect whether vectors or matricies are given
  A <- as.matrix(A)
  B <- as.matrix(B)

  ## IF the input is a vector, use vector formula
  if (nrow(A) == 1 | ncol(A) == 1) {

    ## Check to ensure vectors have the same length
    if ( all.equal(dim(A), dim(B)) != TRUE ) {

      stop("The vectors A and B do not have the same dimensions.")
    } #END if (length(A) != length(B)) {

    ## Cannot find cosine of non-numeric values
    if ( !all(is.numeric(A)) | !all(is.numeric(B)) ) {

      ## Give an error message, stop the function
      stop("The matrices A and/or B must only contain numeric values")

    } # END if (!all(is.numeric(A)) | !all(is.numeric(B))) {

    ## Compute the vector cosine
    cosine <- crossprod(A, B) / (sqrt(crossprod(A, A)) * sqrt(crossprod(B, B)))

  } # END if (is.vector(A)) {

  ## Use matrix multiplication if matrices are specified
  if (nrow(A) > 1 & ncol(A) > 1) {

    ## Check to make sure matrices have the same dimensions
    if ( all.equal(dim(A), dim(B)) != TRUE ) {
      stop("The matrices A and B do not have the same dimensions.")
    }# END if (all.equal(dim(A), dim(B)) == FALSE) {

    ## Ensure only numeric values are given
    if (!all(is.numeric(A)) | !all(is.numeric(B))) {
      stop(" The matrices A and/or B must only contain numeric values")
    } # END if (!all(is.numeric(A)) | !all(is.numeric(B))) {

    ## If specified, align the matrices using a least-squares criterion
    if (align == TRUE) {
      B <- fungible::faAlign(F1 = A,
                             F2 = B)$F2
    } #END  if (align == TRUE)

    ## Diagonal matrices to find the norms of matrices A and B
    A1 <- diag(1 / sqrt(apply(A^2, 2, sum)))
    B1 <- diag(1 / sqrt(apply(B^2, 2, sum)))

    ## Normed A matrix multiplied by the normed B matrix
    cosine <- A1 %*% t(A) %*% B %*% B1

  }#ENDif (nrow(A) > 1 & ncol(A) > 1)

  list(cosine    = round(as.matrix(cosine), digits),
       A         = round(A, digits),
       B         = round(B, digits),
       align      = align)

} # END cosMat <- function(A, B, align = FALSE) {


