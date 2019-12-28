#' @title mediation analysis using R
#' @description mediation analysis using R
#' @param N the number of samples
#' @param thin the number of between-sample random numbers
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' RMonte_antithetic(exp,0,1,10)
#' }
#' @export
RMonte_antithetic <- function(FUN,low,up,N){
  N_2 = ceiling(N/2)
  X_1 = C_unif(N_2,low,up)
  X_2 = up+low - X_1
  X = c(X_1,X_2)
  g_X_antithetic = FUN(X)
  return (mean(g_X_antithetic))
}