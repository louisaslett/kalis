// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// CacheAllSequences2
int CacheAllSequences2(CharacterVector seqs, int bufsize);
RcppExport SEXP StatGen_CacheAllSequences2(SEXP seqsSEXP, SEXP bufsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type seqs(seqsSEXP);
    Rcpp::traits::input_parameter< int >::type bufsize(bufsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(CacheAllSequences2(seqs, bufsize));
    return rcpp_result_gen;
END_RCPP
}
// QueryCache2
IntegerVector QueryCache2(int idx);
RcppExport SEXP StatGen_QueryCache2(SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(QueryCache2(idx));
    return rcpp_result_gen;
END_RCPP
}
// ClearSequenceCache2
void ClearSequenceCache2();
RcppExport SEXP StatGen_ClearSequenceCache2() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    ClearSequenceCache2();
    return R_NilValue;
END_RCPP
}
// Exact_ComputeTable_naive_cpp
NumericMatrix Exact_ComputeTable_naive_cpp(int l, NumericMatrix Pi, NumericVector mu, NumericVector rho);
RcppExport SEXP StatGen_Exact_ComputeTable_naive_cpp(SEXP lSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(Exact_ComputeTable_naive_cpp(l, Pi, mu, rho));
    return rcpp_result_gen;
END_RCPP
}
// ExactBackwardNaiveC_cpp
NumericMatrix ExactBackwardNaiveC_cpp(int t, int L, int N, NumericMatrix Pi, NumericVector mu, NumericVector rho);
RcppExport SEXP StatGen_ExactBackwardNaiveC_cpp(SEXP tSEXP, SEXP LSEXP, SEXP NSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(ExactBackwardNaiveC_cpp(t, L, N, Pi, mu, rho));
    return rcpp_result_gen;
END_RCPP
}
// ExactForwardNaiveC_cpp
NumericMatrix ExactForwardNaiveC_cpp(int t, int L, int N, NumericMatrix Pi, NumericVector mu, NumericVector rho);
RcppExport SEXP StatGen_ExactForwardNaiveC_cpp(SEXP tSEXP, SEXP LSEXP, SEXP NSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(ExactForwardNaiveC_cpp(t, L, N, Pi, mu, rho));
    return rcpp_result_gen;
END_RCPP
}
// ExactForwardYepppExpC_cpp
NumericMatrix ExactForwardYepppExpC_cpp(int t, int L, int N, NumericMatrix Pi, NumericVector mu, NumericVector rho);
RcppExport SEXP StatGen_ExactForwardYepppExpC_cpp(SEXP tSEXP, SEXP LSEXP, SEXP NSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(ExactForwardYepppExpC_cpp(t, L, N, Pi, mu, rho));
    return rcpp_result_gen;
END_RCPP
}
// ExactForwardYepppExpAVX_cpp
NumericMatrix ExactForwardYepppExpAVX_cpp(int t, int L, int N, NumericMatrix Pi, NumericVector mu, NumericVector rho);
RcppExport SEXP StatGen_ExactForwardYepppExpAVX_cpp(SEXP tSEXP, SEXP LSEXP, SEXP NSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(ExactForwardYepppExpAVX_cpp(t, L, N, Pi, mu, rho));
    return rcpp_result_gen;
END_RCPP
}
// ExactForwardYepppExpAVX2_cpp
NumericMatrix ExactForwardYepppExpAVX2_cpp(int t, int L, int N, NumericMatrix Pi, NumericVector mu, NumericVector rho);
RcppExport SEXP StatGen_ExactForwardYepppExpAVX2_cpp(SEXP tSEXP, SEXP LSEXP, SEXP NSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(ExactForwardYepppExpAVX2_cpp(t, L, N, Pi, mu, rho));
    return rcpp_result_gen;
END_RCPP
