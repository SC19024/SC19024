#' @title mediation analysis using R
#' @description mediation analysis using R
#' @param FUN the number of samples
#' @param low the number of between-sample random numbers
#' @param up hhh
#' @param N hhh
#' @param layer hhh
#' @param Method hhh
#' @return a random sample of size \code{n}
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm runif var cov
#' @examples
#' \dontrun{
#' RMonte_stratify(exp,0,1,10)
#' }
#' @export
RMonte <- function(FUN,low,up,N,layer=1,Method='None'){
  N_0 = ceiling(N/layer)
  R = (up-low)/layer
  temp = numeric(layer)
  for(i in 1:layer){
    if(Method=="None"){
      X = runif(N_0,low+(i-1)*R,low+i*R)
      temp[i] = mean(FUN(X)*R)     
    }
    if(Method=='control'){
      X = runif(N_0,low+(i-1)*R,low+i*R)
      m_X = (low*2+(2*i-1)*R)/2
      v_X = R^2/12
      g_X = FUN(X)
      C_star = cov(X,g_X)/v_X
      g_X_control = g_X - C_star*(X-m_X)
      temp[i] = mean(g_X_control*R)
    }
    if(Method=="Antithetic"){
      N_2 = ceiling(N_0/2)
      X_1 = runif(N_2,low,up)
      X_2 = low*2+(2*i-1)*R - X_1
      X = c(X_1,X_2)
      g_X_antithetic = FUN(X)
      temp[i] = mean(g_X_antithetic*R)
    }
    else{
      break
    }
  }
  return(sum(temp))
}