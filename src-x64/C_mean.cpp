#include <Rcpp.h>
#include<iostream>
using namespace Rcpp;
//' @title compute the mean of a vector
//' @description A random variable sampler using Rcpp
//' @param N the number of samples
//' @return a random sample of size \code{n}
//' \dontrun{
//' rnC = C_mean(c(1,2,3,4,5))
//' }
//' @export
// [[Rcpp::export]]
NumericVector C_mean(NumericVector X){
  int SIZE = X.size();
  double SUM = 0;
  double SQUARE_SUM = 0;
  for(int i=0;i<SIZE;i++){
    SUM+=X[i];
    SQUARE_SUM+=X[i]*X[i];
  }
  NumericVector result(2);
  result[0] = SUM/SIZE;
  result[1] = (SQUARE_SUM/SIZE - result[0]*result[0])*SIZE/(SIZE-1);
  return(result);
}