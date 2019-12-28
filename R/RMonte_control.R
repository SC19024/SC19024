#' @title mediation analysis using R
#' @description mediation analysis using R
#' @param N the number of samples
#' @param thin the number of between-sample random numbers
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' RMonte_control(exp,0,1,10)
#' }
#' @export
RMonte_control <- function(FUN,low,up,N){
  X = C_unif(N,low,up)
  m_X = (low+up)/2
  v_X = (up-low)^2/12
  g_X = FUN(X)
  C_star = cov(X,g_X)/v_X
  g_X_control = g_X - C_star*(X-m_X)
  return (mean(g_X_control))
}
