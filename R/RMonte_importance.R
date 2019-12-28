#' @title mediation analysis using R
#' @description mediation analysis using R
#' @param N the number of samples
#' @param thin the number of between-sample random numbers
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' RMonte_importance(exp,0,1,10)
#' }
#' @export
RMonte_inf = function(FUN,low=-Inf,up=Inf,N){
  X = rnorm(N,0,100)
  X = X[(X>low) & (X<up)]
  return(mean(FUN(X)/dnorm(X,0,100))/N*length(X))
}