#' @title mediation analysis using R
#' @description mediation analysis using R
#' @param N the number of samples
#' @param thin the number of between-sample random numbers
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' RMonte(exp,0,1,10)
#' }
#' @export
RMonte <- function(FUN,low,up,N){
  sep = C_unif(N,low,up)
  return(mean(FUN(sep)))
}