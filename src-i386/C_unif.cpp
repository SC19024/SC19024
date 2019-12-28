#include <Rcpp.h>
#include<vector>
using namespace Rcpp;
//' @title A random variable sampler using Rcpp
//' @description A random variable sampler using Rcpp
//' @param N the number of samples
//' @return a random sample of size \code{n}
//' @importFrom Rcpp evalCpp
//' @useDynLib SC19024
//' @examples
//' \dontrun{
//' rnC = C_unif(100)
//' }
//' @export
// [[Rcpp::export]]
NumericVector C_unif(int n, double low, double up) {
  NumericVector Series(n);
  double R = up-low;
  for (int i = 0; i < n; i++) {
    Series[i] = (rand() / double(RAND_MAX) )*R + low;
  }
  return(Series);
}
