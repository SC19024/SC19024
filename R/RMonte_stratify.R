#' @title mediation analysis using R
#' @description mediation analysis using R
#' @param N the number of samples
#' @param thin the number of between-sample random numbers
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' RMonte_stratify(exp,0,1,10,1)
#' }
#' @export
RMonte_stratify <- function(FUN,low,up,N,layer){
  N_0 = ceiling(N/layer)
  R = (up-low)/layer
  temp = numeric(layer)
  for(i in 1:layer){
    X = C_unif(N,low+(i-1)*R,low+i*R)
    temp[i] = mean(FUN(X))
  }
  return(mean(temp))
}