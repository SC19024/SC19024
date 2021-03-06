# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title compute the mean of a vector
#' @description A random variable sampler using Rcpp
#' @param X the number of samples
#' @return a random sample of size \code{n}
#' @export
C_mean <- function(X) {
    .Call('_SC19024_C_mean', PACKAGE = 'SC19024', X)
}

compute_poly <- function(x, COEF) {
    .Call('_SC19024_compute_poly', PACKAGE = 'SC19024', x, COEF)
}

#' @title generate a random variable from an exp-family distribution
#' @description generate a random variable from an exp-family distribution by MCMC
#' @param COEF_h polynomial of hx
#' @param COEF polynomial of the fucntion in exp item
#' @param sigma variance of the random walk
#' @param L length of the markov chain 
#' @param low lower bound of the support set
#' @param up up bound of the support set
#' @return a random sample of size \code{L}
#' @export
#' @examples
#'\dontrun{
#' C_poly(c(1),c(0,0,-0.5),1,10000,-1000,1000)
#' }
C_poly <- function(COEF_h, COEF, sigma, L, low, up) {
    .Call('_SC19024_C_poly', PACKAGE = 'SC19024', COEF_h, COEF, sigma, L, low, up)
}

rcpp_hello_world <- function() {
    .Call('_SC19024_rcpp_hello_world', PACKAGE = 'SC19024')
}

