// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// get_nlines
long get_nlines(std::string fn);
RcppExport SEXP _Siccuracy_get_nlines(SEXP fnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fn(fnSEXP);
    rcpp_result_gen = Rcpp::wrap(get_nlines(fn));
    return rcpp_result_gen;
END_RCPP
}
// convert_phases
long convert_phases(std::string fnin, std::string fnout, long nlines, long na, Rcpp::IntegerVector range);
RcppExport SEXP _Siccuracy_convert_phases(SEXP fninSEXP, SEXP fnoutSEXP, SEXP nlinesSEXP, SEXP naSEXP, SEXP rangeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fnin(fninSEXP);
    Rcpp::traits::input_parameter< std::string >::type fnout(fnoutSEXP);
    Rcpp::traits::input_parameter< long >::type nlines(nlinesSEXP);
    Rcpp::traits::input_parameter< long >::type na(naSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type range(rangeSEXP);
    rcpp_result_gen = Rcpp::wrap(convert_phases(fnin, fnout, nlines, na, range));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Siccuracy_get_nlines", (DL_FUNC) &_Siccuracy_get_nlines, 1},
    {"_Siccuracy_convert_phases", (DL_FUNC) &_Siccuracy_convert_phases, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_Siccuracy(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
