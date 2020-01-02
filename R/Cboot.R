#' @title Cboot
#' @description bootstrap estimation of mean of a vector
#' @param DATA a vector which contain the data
#' @param sims the number of simulations
#' @param Method algorithm for computing confidence interval. sohuld be Norm, perc, basic, stud
#' @param alpha 1 - confidence level
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' Cboot(1:100,10,'Norm',0.05)
#' }
#' @export
Cboot = function(DATA, sims = 10,Method = 'Norm',alpha=0.05){
  N = length(DATA)
  MEAN = numeric(sims)
  VAR = numeric(sims)
  for(i in 1:sims){
    SP = sample(DATA,N,replace = TRUE)
    MV = C_mean(SP)
    MEAN[i] = MV[1]
    VAR[i] = MV[2]
  }
  MEAN_mean = mean(MEAN)
  MEAN_sd = sd(MEAN)
  result = list()
  result[['MEAN']] = MEAN_mean
  result[['SD']] = MEAN_sd
  result[['BIAS']] = MEAN_mean - mean(DATA)
  if(Method=='Norm'){
    result[['INTERVAL']] = c(MEAN_mean - qnorm(alpha/2)*MEAN_sd, MEAN_mean+qnorm(alpha/2)*MEAN_sd)
  }
  if(Method=='Perc'){
    MEAN = sort(MEAN)
    result[['INTERVAL']] = c(MEAN[ceiling(alpha/2*sims)],MEAN[ceiling((1-alpha/2)*sims)])
  }
  if(Method=='Basic'){
    MEAN = sort(MEAN)
    result[['INTERVAL']] = c(2*MEAN_mean - MEAN[ceiling((1-alpha/2)*sims)],2*MEAN_mean - MEAN[ceiling((alpha/2)*sims)])
  }
  if(Method=='stud'){
    t = sort((MEAN -mean(DATA))/sqrt(VAR))
    result[['INTERVAL']] = c(MEAN_mean - MEAN_sd*t[ceiling(1-alpha/2)*sims],MEAN_mean - MEAN_sd*t[ceiling(alpha/2*sims)])
  }
  return(result)
}