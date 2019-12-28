// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// C_gaus
NumericVector C_gaus(int n);
RcppExport SEXP _SC19024_C_gaus(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(C_gaus(n));
    return rcpp_result_gen;
END_RCPP
}
// C_mean
NumericVector C_mean(NumericVector X);
RcppExport SEXP _SC19024_C_mean(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(C_mean(X));
    return rcpp_result_gen;
END_RCPP
}
// C_unif
NumericVector C_unif(int n, double low, double up);
RcppExport SEXP _SC19024_C_unif(SEXP nSEXP, SEXP lowSEXP, SEXP upSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type low(lowSEXP);
    Rcpp::traits::input_parameter< double >::type up(upSEXP);
    rcpp_result_gen = Rcpp::wrap(C_unif(n, low, up));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _SC19024_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SC19024_C_gaus", (DL_FUNC) &_SC19024_C_gaus, 1},
    {"_SC19024_C_mean", (DL_FUNC) &_SC19024_C_mean, 1},
    {"_SC19024_C_unif", (DL_FUNC) &_SC19024_C_unif, 3},
    {"_SC19024_rcpp_hello_world", (DL_FUNC) &_SC19024_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_SC19024(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
