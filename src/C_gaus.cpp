#include <Rcpp.h>
using namespace Rcpp;

double gaussrand()
{
  static double V1, V2, S;
  static int phase = 0;
  double X;
  
  if ( phase == 0 ) {
    do {
      double U1 = (double)rand() / RAND_MAX;
      double U2 = (double)rand() / RAND_MAX;
      
      V1 = 2 * U1 - 1;
      V2 = 2 * U2 - 1;
      S = V1 * V1 + V2 * V2;
    } while(S >= 1 || S == 0);
    
    X = V1 * sqrt(-2 * log(S) / S);
  } else
    X = V2 * sqrt(-2 * log(S) / S);
  
  phase = 1 - phase;
  
  return X;
}

//' @title A random variable sampler using Rcpp
//' @description A random variable sampler using Rcpp
//' @param N the number of samples
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' rnC = C_gaus(100)
//' }
//' @export
// [[Rcpp::export]]
NumericVector C_gaus(int n) {
  NumericVector Series(n);
  for (int i = 0; i < n; i++) {
    Series[i] = gaussrand();
  }
  return(Series);
}
