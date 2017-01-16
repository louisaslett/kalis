// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// Backward
void Backward(List bck, const int t, NumericMatrix Pi, NumericVector mu, NumericVector rho, const int nthreads);
RcppExport SEXP StatGen_Backward(SEXP bckSEXP, SEXP tSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type bck(bckSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Backward(bck, t, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
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
// ComputeStatus
void ComputeStatus();
RcppExport SEXP StatGen_ComputeStatus() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    ComputeStatus();
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
void ExactBackwardNaiveC_cpp(NumericMatrix beta, NumericVector beta_g, NumericVector beta_g2, const int beta_from_rec, const int beta_t, const int t, const int from_rec, const int to_rec, const int L, const int N, NumericMatrix Pi, NumericVector mu, NumericVector rho);
RcppExport SEXP StatGen_ExactBackwardNaiveC_cpp(SEXP betaSEXP, SEXP beta_gSEXP, SEXP beta_g2SEXP, SEXP beta_from_recSEXP, SEXP beta_tSEXP, SEXP tSEXP, SEXP from_recSEXP, SEXP to_recSEXP, SEXP LSEXP, SEXP NSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_g(beta_gSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_g2(beta_g2SEXP);
    Rcpp::traits::input_parameter< const int >::type beta_from_rec(beta_from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int >::type from_rec(from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type to_rec(to_recSEXP);
    Rcpp::traits::input_parameter< const int >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    ExactBackwardNaiveC_cpp(beta, beta_g, beta_g2, beta_from_rec, beta_t, t, from_rec, to_rec, L, N, Pi, mu, rho);
    return R_NilValue;
END_RCPP
}
// ParExactBackwardNaiveC_cpp
void ParExactBackwardNaiveC_cpp(NumericMatrix beta, NumericVector beta_g, NumericVector beta_g2, const int beta_from_rec, const int beta_t, const int t, const int from_rec, const int to_rec, const int L, const int N, NumericMatrix Pi, NumericVector mu, NumericVector rho, const int nthreads);
RcppExport SEXP StatGen_ParExactBackwardNaiveC_cpp(SEXP betaSEXP, SEXP beta_gSEXP, SEXP beta_g2SEXP, SEXP beta_from_recSEXP, SEXP beta_tSEXP, SEXP tSEXP, SEXP from_recSEXP, SEXP to_recSEXP, SEXP LSEXP, SEXP NSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_g(beta_gSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_g2(beta_g2SEXP);
    Rcpp::traits::input_parameter< const int >::type beta_from_rec(beta_from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int >::type from_rec(from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type to_rec(to_recSEXP);
    Rcpp::traits::input_parameter< const int >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    ParExactBackwardNaiveC_cpp(beta, beta_g, beta_g2, beta_from_rec, beta_t, t, from_rec, to_rec, L, N, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
// ExactBackwardNoExpAVX3_cpp
void ExactBackwardNoExpAVX3_cpp(NumericMatrix beta, NumericVector beta_g, NumericVector beta_g2, const int beta_from_rec, const int beta_t, const int t, const int from_rec, const int to_rec, const int L, const int N, NumericMatrix Pi, NumericVector mu, NumericVector rho);
RcppExport SEXP StatGen_ExactBackwardNoExpAVX3_cpp(SEXP betaSEXP, SEXP beta_gSEXP, SEXP beta_g2SEXP, SEXP beta_from_recSEXP, SEXP beta_tSEXP, SEXP tSEXP, SEXP from_recSEXP, SEXP to_recSEXP, SEXP LSEXP, SEXP NSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_g(beta_gSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_g2(beta_g2SEXP);
    Rcpp::traits::input_parameter< const int >::type beta_from_rec(beta_from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int >::type from_rec(from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type to_rec(to_recSEXP);
    Rcpp::traits::input_parameter< const int >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    ExactBackwardNoExpAVX3_cpp(beta, beta_g, beta_g2, beta_from_rec, beta_t, t, from_rec, to_rec, L, N, Pi, mu, rho);
    return R_NilValue;
END_RCPP
}
// ParExactBackwardNoExpAVX3_cpp
void ParExactBackwardNoExpAVX3_cpp(NumericMatrix beta, NumericVector beta_g, NumericVector beta_g2, const int beta_from_rec, const int beta_t, const int t, const int from_rec, const int to_rec, const int L, const int N, NumericMatrix Pi, NumericVector mu, NumericVector rho, const int nthreads);
RcppExport SEXP StatGen_ParExactBackwardNoExpAVX3_cpp(SEXP betaSEXP, SEXP beta_gSEXP, SEXP beta_g2SEXP, SEXP beta_from_recSEXP, SEXP beta_tSEXP, SEXP tSEXP, SEXP from_recSEXP, SEXP to_recSEXP, SEXP LSEXP, SEXP NSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_g(beta_gSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_g2(beta_g2SEXP);
    Rcpp::traits::input_parameter< const int >::type beta_from_rec(beta_from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int >::type from_rec(from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type to_rec(to_recSEXP);
    Rcpp::traits::input_parameter< const int >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    ParExactBackwardNoExpAVX3_cpp(beta, beta_g, beta_g2, beta_from_rec, beta_t, t, from_rec, to_rec, L, N, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
// ExactBackwardNoExpAVX3single_cpp
void ExactBackwardNoExpAVX3single_cpp(NumericMatrix beta, NumericVector beta_g, NumericVector beta_g2, const int beta_from_rec, const int beta_t, const int t, const int from_rec, const int to_rec, const int L, const int N, NumericMatrix Pi, NumericVector mu, NumericVector rho);
RcppExport SEXP StatGen_ExactBackwardNoExpAVX3single_cpp(SEXP betaSEXP, SEXP beta_gSEXP, SEXP beta_g2SEXP, SEXP beta_from_recSEXP, SEXP beta_tSEXP, SEXP tSEXP, SEXP from_recSEXP, SEXP to_recSEXP, SEXP LSEXP, SEXP NSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_g(beta_gSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_g2(beta_g2SEXP);
    Rcpp::traits::input_parameter< const int >::type beta_from_rec(beta_from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int >::type from_rec(from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type to_rec(to_recSEXP);
    Rcpp::traits::input_parameter< const int >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    ExactBackwardNoExpAVX3single_cpp(beta, beta_g, beta_g2, beta_from_rec, beta_t, t, from_rec, to_rec, L, N, Pi, mu, rho);
    return R_NilValue;
END_RCPP
}
// ExactForwardNaiveC_cpp
void ExactForwardNaiveC_cpp(NumericMatrix alpha, NumericVector alpha_f, NumericVector alpha_f2, const int alpha_from_rec, const int alpha_t, const int t, const int from_rec, const int to_rec, const int L, const int N, NumericMatrix Pi, NumericVector mu, NumericVector rho);
RcppExport SEXP StatGen_ExactForwardNaiveC_cpp(SEXP alphaSEXP, SEXP alpha_fSEXP, SEXP alpha_f2SEXP, SEXP alpha_from_recSEXP, SEXP alpha_tSEXP, SEXP tSEXP, SEXP from_recSEXP, SEXP to_recSEXP, SEXP LSEXP, SEXP NSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha_f(alpha_fSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha_f2(alpha_f2SEXP);
    Rcpp::traits::input_parameter< const int >::type alpha_from_rec(alpha_from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type alpha_t(alpha_tSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int >::type from_rec(from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type to_rec(to_recSEXP);
    Rcpp::traits::input_parameter< const int >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    ExactForwardNaiveC_cpp(alpha, alpha_f, alpha_f2, alpha_from_rec, alpha_t, t, from_rec, to_rec, L, N, Pi, mu, rho);
    return R_NilValue;
END_RCPP
}
// ParExactForwardNaiveC_cpp
void ParExactForwardNaiveC_cpp(NumericMatrix alpha, NumericVector alpha_f, NumericVector alpha_f2, const int alpha_from_rec, const int alpha_t, const int t, const int from_rec, const int to_rec, const int L, const int N, NumericMatrix Pi, NumericVector mu, NumericVector rho, const int nthreads);
RcppExport SEXP StatGen_ParExactForwardNaiveC_cpp(SEXP alphaSEXP, SEXP alpha_fSEXP, SEXP alpha_f2SEXP, SEXP alpha_from_recSEXP, SEXP alpha_tSEXP, SEXP tSEXP, SEXP from_recSEXP, SEXP to_recSEXP, SEXP LSEXP, SEXP NSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha_f(alpha_fSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha_f2(alpha_f2SEXP);
    Rcpp::traits::input_parameter< const int >::type alpha_from_rec(alpha_from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type alpha_t(alpha_tSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int >::type from_rec(from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type to_rec(to_recSEXP);
    Rcpp::traits::input_parameter< const int >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    ParExactForwardNaiveC_cpp(alpha, alpha_f, alpha_f2, alpha_from_rec, alpha_t, t, from_rec, to_rec, L, N, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
// ExactForwardNoExpAVX3_cpp
void ExactForwardNoExpAVX3_cpp(NumericMatrix alpha, NumericVector alpha_f, NumericVector alpha_f2, const int alpha_from_rec, const int alpha_t, const int t, const int from_rec, const int to_rec, const int L, const int N, NumericMatrix Pi, NumericVector mu, NumericVector rho);
RcppExport SEXP StatGen_ExactForwardNoExpAVX3_cpp(SEXP alphaSEXP, SEXP alpha_fSEXP, SEXP alpha_f2SEXP, SEXP alpha_from_recSEXP, SEXP alpha_tSEXP, SEXP tSEXP, SEXP from_recSEXP, SEXP to_recSEXP, SEXP LSEXP, SEXP NSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha_f(alpha_fSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha_f2(alpha_f2SEXP);
    Rcpp::traits::input_parameter< const int >::type alpha_from_rec(alpha_from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type alpha_t(alpha_tSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int >::type from_rec(from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type to_rec(to_recSEXP);
    Rcpp::traits::input_parameter< const int >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    ExactForwardNoExpAVX3_cpp(alpha, alpha_f, alpha_f2, alpha_from_rec, alpha_t, t, from_rec, to_rec, L, N, Pi, mu, rho);
    return R_NilValue;
END_RCPP
}
// ParExactForwardNoExpAVX3_cpp
void ParExactForwardNoExpAVX3_cpp(NumericMatrix alpha, NumericVector alpha_f, NumericVector alpha_f2, const int alpha_from_rec, const int alpha_t, const int t, const int from_rec, const int to_rec, const int L, const int N, NumericMatrix Pi, NumericVector mu, NumericVector rho, const int nthreads);
RcppExport SEXP StatGen_ParExactForwardNoExpAVX3_cpp(SEXP alphaSEXP, SEXP alpha_fSEXP, SEXP alpha_f2SEXP, SEXP alpha_from_recSEXP, SEXP alpha_tSEXP, SEXP tSEXP, SEXP from_recSEXP, SEXP to_recSEXP, SEXP LSEXP, SEXP NSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha_f(alpha_fSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha_f2(alpha_f2SEXP);
    Rcpp::traits::input_parameter< const int >::type alpha_from_rec(alpha_from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type alpha_t(alpha_tSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int >::type from_rec(from_recSEXP);
    Rcpp::traits::input_parameter< const int >::type to_rec(to_recSEXP);
    Rcpp::traits::input_parameter< const int >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    ParExactForwardNoExpAVX3_cpp(alpha, alpha_f, alpha_f2, alpha_from_rec, alpha_t, t, from_rec, to_rec, L, N, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
// ResetForwardTable
void ResetForwardTable(List fwd);
RcppExport SEXP StatGen_ResetForwardTable(SEXP fwdSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type fwd(fwdSEXP);
    ResetForwardTable(fwd);
    return R_NilValue;
END_RCPP
}
// Forward
void Forward(List fwd, const int t, NumericMatrix Pi, NumericVector mu, NumericVector rho, const int nthreads);
RcppExport SEXP StatGen_Forward(SEXP fwdSEXP, SEXP tSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type fwd(fwdSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Forward(fwd, t, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
// FillTableCache
void FillTableCache(List cache, NumericMatrix Pi, NumericVector mu, NumericVector rho, const int nthreads, int from, int to);
RcppExport SEXP StatGen_FillTableCache(SEXP cacheSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type cache(cacheSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< int >::type from(fromSEXP);
    Rcpp::traits::input_parameter< int >::type to(toSEXP);
    FillTableCache(cache, Pi, mu, rho, nthreads, from, to);
    return R_NilValue;
END_RCPP
}
// CopyForwardTable
void CopyForwardTable(List to, List from);
RcppExport SEXP StatGen_CopyForwardTable(SEXP toSEXP, SEXP fromSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type to(toSEXP);
    Rcpp::traits::input_parameter< List >::type from(fromSEXP);
    CopyForwardTable(to, from);
    return R_NilValue;
END_RCPP
}
