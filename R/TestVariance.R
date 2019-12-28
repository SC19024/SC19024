#' @title TestVariance
#' @description compute the variance of an estimation
#' @param TESTFUN the number of samples
#' @param M the number of simulations
#' @param ... other parameters
#' @return real value number.
#' @examples
#' \dontrun{
#' TestVariance(RMonte,10,exp,1,5,100)
#' }
#' @export
TestVariance <- function(TESTFUN,M,...){
  result = numeric(M)
  for(i in 1:M){
    result[i] = TESTFUN(...)
  }
  return(var(result))
}