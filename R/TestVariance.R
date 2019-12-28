#' @title mediation analysis using R
#' @description mediation analysis using R
#' @param TESTFUN the number of samples
#' @param M the number of between-sample random numbers
#' @param ... hhh
#' @return a random sample of size \code{n}
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