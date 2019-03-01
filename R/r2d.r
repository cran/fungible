# convert radians to degrees


#' Convert Radians to Degrees
#' 
#' Convert radian measure to degrees.
#' 
#' 
#' @param radian Radian measure of an angle
#' @return \item{}{Degree measure of an angle}
#' @keywords Statistics
#' @export
#' @examples
#' 
#'  r2d(.5*pi)
#' 
  r2d <- function(radian) radian*180/pi
