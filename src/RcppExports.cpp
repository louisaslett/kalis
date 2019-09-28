// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ResetBackwardTable
void ResetBackwardTable(List bck);
RcppExport SEXP _kalis_ResetBackwardTable(SEXP bckSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type bck(bckSEXP);
    ResetBackwardTable(bck);
    return R_NilValue;
END_RCPP
}
// Backward_densePi_densemu_cpp
void Backward_densePi_densemu_cpp(List bck, const int t, NumericMatrix Pi, NumericVector mu, NumericVector rho, const int nthreads);
RcppExport SEXP _kalis_Backward_densePi_densemu_cpp(SEXP bckSEXP, SEXP tSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type bck(bckSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Backward_densePi_densemu_cpp(bck, t, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
// Backward_scalarPi_densemu_cpp
void Backward_scalarPi_densemu_cpp(List bck, const int t, const double Pi, NumericVector mu, NumericVector rho, const int nthreads);
RcppExport SEXP _kalis_Backward_scalarPi_densemu_cpp(SEXP bckSEXP, SEXP tSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type bck(bckSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Backward_scalarPi_densemu_cpp(bck, t, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
// Backward_densePi_scalarmu_cpp
void Backward_densePi_scalarmu_cpp(List bck, const int t, NumericMatrix Pi, const double mu, NumericVector rho, const int nthreads);
RcppExport SEXP _kalis_Backward_densePi_scalarmu_cpp(SEXP bckSEXP, SEXP tSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type bck(bckSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Backward_densePi_scalarmu_cpp(bck, t, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
// Backward_scalarPi_scalarmu_cpp
void Backward_scalarPi_scalarmu_cpp(List bck, const int t, const double Pi, const double mu, NumericVector rho, const int nthreads);
RcppExport SEXP _kalis_Backward_scalarPi_scalarmu_cpp(SEXP bckSEXP, SEXP tSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type bck(bckSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Backward_scalarPi_scalarmu_cpp(bck, t, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
// CacheHaplotypes_matrix_2
int CacheHaplotypes_matrix_2(IntegerMatrix x, int N, int L, int transpose);
RcppExport SEXP _kalis_CacheHaplotypes_matrix_2(SEXP xSEXP, SEXP NSEXP, SEXP LSEXP, SEXP transposeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type transpose(transposeSEXP);
    rcpp_result_gen = Rcpp::wrap(CacheHaplotypes_matrix_2(x, N, L, transpose));
    return rcpp_result_gen;
END_RCPP
}
// CacheHaplotypes_hdf5_2
int CacheHaplotypes_hdf5_2(Function nexthaps, int N, int L);
RcppExport SEXP _kalis_CacheHaplotypes_hdf5_2(SEXP nexthapsSEXP, SEXP NSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type nexthaps(nexthapsSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(CacheHaplotypes_hdf5_2(nexthaps, N, L));
    return rcpp_result_gen;
END_RCPP
}
// QueryCache2_ind
IntegerVector QueryCache2_ind(int idx);
RcppExport SEXP _kalis_QueryCache2_ind(SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(QueryCache2_ind(idx));
    return rcpp_result_gen;
END_RCPP
}
// QueryCache2_loc
IntegerVector QueryCache2_loc(int idx);
RcppExport SEXP _kalis_QueryCache2_loc(SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(QueryCache2_loc(idx));
    return rcpp_result_gen;
END_RCPP
}
// ClearHaplotypeCache2
void ClearHaplotypeCache2();
RcppExport SEXP _kalis_ClearHaplotypeCache2() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    ClearHaplotypeCache2();
    return R_NilValue;
END_RCPP
}
// ComputeStatus
void ComputeStatus();
RcppExport SEXP _kalis_ComputeStatus() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    ComputeStatus();
    return R_NilValue;
END_RCPP
}
// Dedip_min
NumericVector Dedip_min(NumericMatrix fwd, NumericMatrix bck, NumericVector s);
RcppExport SEXP _kalis_Dedip_min(SEXP fwdSEXP, SEXP bckSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type fwd(fwdSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bck(bckSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(Dedip_min(fwd, bck, s));
    return rcpp_result_gen;
END_RCPP
}
// Dedip_2nd_min
NumericVector Dedip_2nd_min(NumericMatrix fwd, NumericMatrix bck, NumericVector s);
RcppExport SEXP _kalis_Dedip_2nd_min(SEXP fwdSEXP, SEXP bckSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type fwd(fwdSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bck(bckSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(Dedip_2nd_min(fwd, bck, s));
    return rcpp_result_gen;
END_RCPP
}
// Dedip_max
NumericVector Dedip_max(NumericMatrix fwd, NumericMatrix bck, NumericVector s);
RcppExport SEXP _kalis_Dedip_max(SEXP fwdSEXP, SEXP bckSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type fwd(fwdSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bck(bckSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(Dedip_max(fwd, bck, s));
    return rcpp_result_gen;
END_RCPP
}
// Dedip_dom
NumericVector Dedip_dom(NumericMatrix fwd, NumericMatrix bck, NumericVector s);
RcppExport SEXP _kalis_Dedip_dom(SEXP fwdSEXP, SEXP bckSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type fwd(fwdSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bck(bckSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(Dedip_dom(fwd, bck, s));
    return rcpp_result_gen;
END_RCPP
}
// Dedip_add
NumericVector Dedip_add(NumericMatrix fwd, NumericMatrix bck, NumericVector s);
RcppExport SEXP _kalis_Dedip_add(SEXP fwdSEXP, SEXP bckSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type fwd(fwdSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bck(bckSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(Dedip_add(fwd, bck, s));
    return rcpp_result_gen;
END_RCPP
}
// Dedip_mean
NumericVector Dedip_mean(NumericMatrix fwd, NumericMatrix bck, NumericVector s);
RcppExport SEXP _kalis_Dedip_mean(SEXP fwdSEXP, SEXP bckSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type fwd(fwdSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bck(bckSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(Dedip_mean(fwd, bck, s));
    return rcpp_result_gen;
END_RCPP
}
// Dedip_all
List Dedip_all(NumericMatrix fwd, NumericMatrix bck, NumericVector s);
RcppExport SEXP _kalis_Dedip_all(SEXP fwdSEXP, SEXP bckSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type fwd(fwdSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bck(bckSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(Dedip_all(fwd, bck, s));
    return rcpp_result_gen;
END_RCPP
}
// Dedip2_min
NumericVector Dedip2_min(NumericMatrix M);
RcppExport SEXP _kalis_Dedip2_min(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(Dedip2_min(M));
    return rcpp_result_gen;
END_RCPP
}
// Dedip2_2nd_min
NumericVector Dedip2_2nd_min(NumericMatrix M);
RcppExport SEXP _kalis_Dedip2_2nd_min(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(Dedip2_2nd_min(M));
    return rcpp_result_gen;
END_RCPP
}
// Dedip2_dom
NumericVector Dedip2_dom(NumericMatrix M);
RcppExport SEXP _kalis_Dedip2_dom(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(Dedip2_dom(M));
    return rcpp_result_gen;
END_RCPP
}
// Dedip2_all
List Dedip2_all(NumericMatrix M);
RcppExport SEXP _kalis_Dedip2_all(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(Dedip2_all(M));
    return rcpp_result_gen;
END_RCPP
}
// DedipAndMul
NumericVector DedipAndMul(NumericMatrix M, NumericMatrix alpha, NumericMatrix beta, NumericVector x, int from_recipient, int nthreads);
RcppExport SEXP _kalis_DedipAndMul(SEXP MSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP xSEXP, SEXP from_recipientSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type from_recipient(from_recipientSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(DedipAndMul(M, alpha, beta, x, from_recipient, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// DedipOnlyMul
NumericVector DedipOnlyMul(NumericMatrix M, NumericVector x, size_t from_recipient, size_t nthreads);
RcppExport SEXP _kalis_DedipOnlyMul(SEXP MSEXP, SEXP xSEXP, SEXP from_recipientSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< size_t >::type from_recipient(from_recipientSEXP);
    Rcpp::traits::input_parameter< size_t >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(DedipOnlyMul(M, x, from_recipient, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// ResetForwardTable
void ResetForwardTable(List fwd);
RcppExport SEXP _kalis_ResetForwardTable(SEXP fwdSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type fwd(fwdSEXP);
    ResetForwardTable(fwd);
    return R_NilValue;
END_RCPP
}
// Forward_densePi_densemu_cpp
void Forward_densePi_densemu_cpp(List fwd, const int t, NumericMatrix Pi, NumericVector mu, NumericVector rho, const int nthreads);
RcppExport SEXP _kalis_Forward_densePi_densemu_cpp(SEXP fwdSEXP, SEXP tSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type fwd(fwdSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Forward_densePi_densemu_cpp(fwd, t, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
// Forward_scalarPi_densemu_cpp
void Forward_scalarPi_densemu_cpp(List fwd, const int t, const double Pi, NumericVector mu, NumericVector rho, const int nthreads);
RcppExport SEXP _kalis_Forward_scalarPi_densemu_cpp(SEXP fwdSEXP, SEXP tSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type fwd(fwdSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Forward_scalarPi_densemu_cpp(fwd, t, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
// Forward_densePi_scalarmu_cpp
void Forward_densePi_scalarmu_cpp(List fwd, const int t, NumericMatrix Pi, const double mu, NumericVector rho, const int nthreads);
RcppExport SEXP _kalis_Forward_densePi_scalarmu_cpp(SEXP fwdSEXP, SEXP tSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type fwd(fwdSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Forward_densePi_scalarmu_cpp(fwd, t, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
// Forward_scalarPi_scalarmu_cpp
void Forward_scalarPi_scalarmu_cpp(List fwd, const int t, const double Pi, const double mu, NumericVector rho, const int nthreads);
RcppExport SEXP _kalis_Forward_scalarPi_scalarmu_cpp(SEXP fwdSEXP, SEXP tSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type fwd(fwdSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Forward_scalarPi_scalarmu_cpp(fwd, t, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
// Forward1step_scalarPi_scalarmu_cpp
void Forward1step_scalarPi_scalarmu_cpp(List fwd, const int t, const double Pi, const double mu, NumericVector rho, const int nthreads);
RcppExport SEXP _kalis_Forward1step_scalarPi_scalarmu_cpp(SEXP fwdSEXP, SEXP tSEXP, SEXP PiSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type fwd(fwdSEXP);
    Rcpp::traits::input_parameter< const int >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Forward1step_scalarPi_scalarmu_cpp(fwd, t, Pi, mu, rho, nthreads);
    return R_NilValue;
END_RCPP
}
// ResetTable
void ResetTable(List tbl);
RcppExport SEXP _kalis_ResetTable(SEXP tblSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type tbl(tblSEXP);
    ResetTable(tbl);
    return R_NilValue;
END_RCPP
}
// CopyForwardTable_cpp
void CopyForwardTable_cpp(List to, List from);
RcppExport SEXP _kalis_CopyForwardTable_cpp(SEXP toSEXP, SEXP fromSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type to(toSEXP);
    Rcpp::traits::input_parameter< List >::type from(fromSEXP);
    CopyForwardTable_cpp(to, from);
    return R_NilValue;
END_RCPP
}
// CopyBackwardTable_cpp
void CopyBackwardTable_cpp(List to, List from);
RcppExport SEXP _kalis_CopyBackwardTable_cpp(SEXP toSEXP, SEXP fromSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type to(toSEXP);
    Rcpp::traits::input_parameter< List >::type from(fromSEXP);
    CopyBackwardTable_cpp(to, from);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_kalis_ResetBackwardTable", (DL_FUNC) &_kalis_ResetBackwardTable, 1},
    {"_kalis_Backward_densePi_densemu_cpp", (DL_FUNC) &_kalis_Backward_densePi_densemu_cpp, 6},
    {"_kalis_Backward_scalarPi_densemu_cpp", (DL_FUNC) &_kalis_Backward_scalarPi_densemu_cpp, 6},
    {"_kalis_Backward_densePi_scalarmu_cpp", (DL_FUNC) &_kalis_Backward_densePi_scalarmu_cpp, 6},
    {"_kalis_Backward_scalarPi_scalarmu_cpp", (DL_FUNC) &_kalis_Backward_scalarPi_scalarmu_cpp, 6},
    {"_kalis_CacheHaplotypes_matrix_2", (DL_FUNC) &_kalis_CacheHaplotypes_matrix_2, 4},
    {"_kalis_CacheHaplotypes_hdf5_2", (DL_FUNC) &_kalis_CacheHaplotypes_hdf5_2, 3},
    {"_kalis_QueryCache2_ind", (DL_FUNC) &_kalis_QueryCache2_ind, 1},
    {"_kalis_QueryCache2_loc", (DL_FUNC) &_kalis_QueryCache2_loc, 1},
    {"_kalis_ClearHaplotypeCache2", (DL_FUNC) &_kalis_ClearHaplotypeCache2, 0},
    {"_kalis_ComputeStatus", (DL_FUNC) &_kalis_ComputeStatus, 0},
    {"_kalis_Dedip_min", (DL_FUNC) &_kalis_Dedip_min, 3},
    {"_kalis_Dedip_2nd_min", (DL_FUNC) &_kalis_Dedip_2nd_min, 3},
    {"_kalis_Dedip_max", (DL_FUNC) &_kalis_Dedip_max, 3},
    {"_kalis_Dedip_dom", (DL_FUNC) &_kalis_Dedip_dom, 3},
    {"_kalis_Dedip_add", (DL_FUNC) &_kalis_Dedip_add, 3},
    {"_kalis_Dedip_mean", (DL_FUNC) &_kalis_Dedip_mean, 3},
    {"_kalis_Dedip_all", (DL_FUNC) &_kalis_Dedip_all, 3},
    {"_kalis_Dedip2_min", (DL_FUNC) &_kalis_Dedip2_min, 1},
    {"_kalis_Dedip2_2nd_min", (DL_FUNC) &_kalis_Dedip2_2nd_min, 1},
    {"_kalis_Dedip2_dom", (DL_FUNC) &_kalis_Dedip2_dom, 1},
    {"_kalis_Dedip2_all", (DL_FUNC) &_kalis_Dedip2_all, 1},
    {"_kalis_DedipAndMul", (DL_FUNC) &_kalis_DedipAndMul, 6},
    {"_kalis_DedipOnlyMul", (DL_FUNC) &_kalis_DedipOnlyMul, 4},
    {"_kalis_ResetForwardTable", (DL_FUNC) &_kalis_ResetForwardTable, 1},
    {"_kalis_Forward_densePi_densemu_cpp", (DL_FUNC) &_kalis_Forward_densePi_densemu_cpp, 6},
    {"_kalis_Forward_scalarPi_densemu_cpp", (DL_FUNC) &_kalis_Forward_scalarPi_densemu_cpp, 6},
    {"_kalis_Forward_densePi_scalarmu_cpp", (DL_FUNC) &_kalis_Forward_densePi_scalarmu_cpp, 6},
    {"_kalis_Forward_scalarPi_scalarmu_cpp", (DL_FUNC) &_kalis_Forward_scalarPi_scalarmu_cpp, 6},
    {"_kalis_Forward1step_scalarPi_scalarmu_cpp", (DL_FUNC) &_kalis_Forward1step_scalarPi_scalarmu_cpp, 6},
    {"_kalis_ResetTable", (DL_FUNC) &_kalis_ResetTable, 1},
    {"_kalis_CopyForwardTable_cpp", (DL_FUNC) &_kalis_CopyForwardTable_cpp, 2},
    {"_kalis_CopyBackwardTable_cpp", (DL_FUNC) &_kalis_CopyBackwardTable_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_kalis(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
