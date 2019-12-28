#' @title mediation analysis using R
#' @description mediation analysis using R
#' @param N the number of samples
#' @param thin the number of between-sample random numbers
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' TestVariance(RMonte,1000,exp,1,5,1000)
#' }
#' @export
TestVariance <- function(TESTFUN,M,...){
  result = numeric(M)
  for(i in 1:M){
    result[i] = FUN(...)
  }
  return(var(result))
}