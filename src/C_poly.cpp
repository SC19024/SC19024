#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;
// [[Rcpp::export]]
double compute_poly(double x, NumericVector COEF){
  int L_COEF = COEF.size();
  double result = 0;
  double temp = 1;
  for(int i = 0;i<L_COEF;i++){
    result += COEF[i] * temp;
    temp*=x;
  }
  return result;
}
//' @title generate a random variable from an exp-family distribution
//' @description generate a random variable from an exp-family distribution by MCMC
//' @param COEF_h polynomial of hx
//' @param COEF polynomial of the fucntion in exp item
//' @param sigma variance of the random walk
//' @param L length of the markov chain 
//' @param low lower bound of the support set
//' @param up up bound of the support set
//' @return a random sample of size \code{L}
//' @export
//' @examples
//'\dontrun{
//' C_poly(c(1),c(0,0,-0.5),1,10000,-1000,1000)
//' }
// [[Rcpp::export]]
NumericVector C_poly(NumericVector COEF_h,NumericVector COEF, double sigma, int L,
                     double low, double up){
  NumericVector Series (L);
  NumericVector Uniform = runif(L);
  int k=0;
  for (int i =1; i<L; i++){
    NumericVector y = rnorm(1,Series[i-1],sigma);
    if(y[0]<up && y[0]>low){
      if(Uniform[i]<=compute_poly(y[0],COEF_h)/compute_poly(Series[i-1],COEF_h)*exp(compute_poly(y[0],COEF)-compute_poly(Series[i-1],COEF))){
        Series[i] = y[0];
        k=k+1;
      }
      else{
        Series[i] = Series[i-1];
      }
    }
    else{
      Series[i] = Series[i-1];
    }
  }
  return Series;
}

