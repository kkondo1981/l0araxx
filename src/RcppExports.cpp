// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// l0araxxC
List l0araxxC(arma::mat X, arma::vec y, arma::vec weights, arma::vec offset, String family, double lambda, int maxit, double eps);
RcppExport SEXP _libaglm_l0araxxC(SEXP XSEXP, SEXP ySEXP, SEXP weightsSEXP, SEXP offsetSEXP, SEXP familySEXP, SEXP lambdaSEXP, SEXP maxitSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< String >::type family(familySEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(l0araxxC(X, y, weights, offset, family, lambda, maxit, eps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_libaglm_l0araxxC", (DL_FUNC) &_libaglm_l0araxxC, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_libaglm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
