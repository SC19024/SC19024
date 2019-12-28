#' @title RMonte
#' @description Monte Carlo integration 
#' @param FUN the function want to integrate
#' @param low the lower bound of the integration interval
#' @param up the up bound of the integration interval
#' @param N time of simulations
#' @param layer the layer of stratification
#' @param Method the variance-reduction method,None,control,Antithetic
#' @return the value of integration
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm runif var cov sd qnorm
#' @useDynLib SC19024
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
  }
  return(sum(temp))
}