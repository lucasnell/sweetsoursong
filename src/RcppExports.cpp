// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// one_flower_ode
NumericMatrix one_flower_ode(const double& dt, const double& max_t, const double& Y_delay, const double& B_delay, const double& Y0, const double& B0, const double& A0, const double& H0, const double& D, double A_0, const double& r_Y, const double& r_B, const double& m_Y, const double& m_B, const double& e_B, const double& q_Y, const double& q_B, const double& c_Y, const double& c_B, const double& h_B, const double& h_Y);
RcppExport SEXP _sweetsoursong_one_flower_ode(SEXP dtSEXP, SEXP max_tSEXP, SEXP Y_delaySEXP, SEXP B_delaySEXP, SEXP Y0SEXP, SEXP B0SEXP, SEXP A0SEXP, SEXP H0SEXP, SEXP DSEXP, SEXP A_0SEXP, SEXP r_YSEXP, SEXP r_BSEXP, SEXP m_YSEXP, SEXP m_BSEXP, SEXP e_BSEXP, SEXP q_YSEXP, SEXP q_BSEXP, SEXP c_YSEXP, SEXP c_BSEXP, SEXP h_BSEXP, SEXP h_YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< const double& >::type max_t(max_tSEXP);
    Rcpp::traits::input_parameter< const double& >::type Y_delay(Y_delaySEXP);
    Rcpp::traits::input_parameter< const double& >::type B_delay(B_delaySEXP);
    Rcpp::traits::input_parameter< const double& >::type Y0(Y0SEXP);
    Rcpp::traits::input_parameter< const double& >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< const double& >::type A0(A0SEXP);
    Rcpp::traits::input_parameter< const double& >::type H0(H0SEXP);
    Rcpp::traits::input_parameter< const double& >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type A_0(A_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type r_Y(r_YSEXP);
    Rcpp::traits::input_parameter< const double& >::type r_B(r_BSEXP);
    Rcpp::traits::input_parameter< const double& >::type m_Y(m_YSEXP);
    Rcpp::traits::input_parameter< const double& >::type m_B(m_BSEXP);
    Rcpp::traits::input_parameter< const double& >::type e_B(e_BSEXP);
    Rcpp::traits::input_parameter< const double& >::type q_Y(q_YSEXP);
    Rcpp::traits::input_parameter< const double& >::type q_B(q_BSEXP);
    Rcpp::traits::input_parameter< const double& >::type c_Y(c_YSEXP);
    Rcpp::traits::input_parameter< const double& >::type c_B(c_BSEXP);
    Rcpp::traits::input_parameter< const double& >::type h_B(h_BSEXP);
    Rcpp::traits::input_parameter< const double& >::type h_Y(h_YSEXP);
    rcpp_result_gen = Rcpp::wrap(one_flower_ode(dt, max_t, Y_delay, B_delay, Y0, B0, A0, H0, D, A_0, r_Y, r_B, m_Y, m_B, e_B, q_Y, q_B, c_Y, c_B, h_B, h_Y));
    return rcpp_result_gen;
END_RCPP
}
// plant_metacomm_stoch_cpp
NumericMatrix plant_metacomm_stoch_cpp(const uint32_t& n_reps, const std::vector<double>& m, const std::vector<double>& d_yp, const std::vector<double>& d_b0, const std::vector<double>& d_bp, const std::vector<double>& g_yp, const std::vector<double>& g_b0, const std::vector<double>& g_bp, const std::vector<double>& L_0, const double& u, const double& X, const std::vector<double>& Y0, const std::vector<double>& B0, const double& n_sigma, const double& season_len, const double& season_surv, const int& rand_season, const bool& open_sys, const double& dt, const double& max_t, const double& burnin);
RcppExport SEXP _sweetsoursong_plant_metacomm_stoch_cpp(SEXP n_repsSEXP, SEXP mSEXP, SEXP d_ypSEXP, SEXP d_b0SEXP, SEXP d_bpSEXP, SEXP g_ypSEXP, SEXP g_b0SEXP, SEXP g_bpSEXP, SEXP L_0SEXP, SEXP uSEXP, SEXP XSEXP, SEXP Y0SEXP, SEXP B0SEXP, SEXP n_sigmaSEXP, SEXP season_lenSEXP, SEXP season_survSEXP, SEXP rand_seasonSEXP, SEXP open_sysSEXP, SEXP dtSEXP, SEXP max_tSEXP, SEXP burninSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const uint32_t& >::type n_reps(n_repsSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type d_yp(d_ypSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type d_b0(d_b0SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type d_bp(d_bpSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type g_yp(g_ypSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type g_b0(g_b0SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type g_bp(g_bpSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type L_0(L_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type Y0(Y0SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< const double& >::type n_sigma(n_sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type season_len(season_lenSEXP);
    Rcpp::traits::input_parameter< const double& >::type season_surv(season_survSEXP);
    Rcpp::traits::input_parameter< const int& >::type rand_season(rand_seasonSEXP);
    Rcpp::traits::input_parameter< const bool& >::type open_sys(open_sysSEXP);
    Rcpp::traits::input_parameter< const double& >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< const double& >::type max_t(max_tSEXP);
    Rcpp::traits::input_parameter< const double& >::type burnin(burninSEXP);
    rcpp_result_gen = Rcpp::wrap(plant_metacomm_stoch_cpp(n_reps, m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, u, X, Y0, B0, n_sigma, season_len, season_surv, rand_season, open_sys, dt, max_t, burnin));
    return rcpp_result_gen;
END_RCPP
}
// plant_metacomm_cpp
arma::mat plant_metacomm_cpp(const std::vector<double>& m, const std::vector<double>& d_yp, const std::vector<double>& d_b0, const std::vector<double>& d_bp, const std::vector<double>& g_yp, const std::vector<double>& g_b0, const std::vector<double>& g_bp, const std::vector<double>& L_0, const double& u, const double& X, const std::vector<double>& Y0, const std::vector<double>& B0, const bool& open_sys, const double& dt, const double& max_t);
RcppExport SEXP _sweetsoursong_plant_metacomm_cpp(SEXP mSEXP, SEXP d_ypSEXP, SEXP d_b0SEXP, SEXP d_bpSEXP, SEXP g_ypSEXP, SEXP g_b0SEXP, SEXP g_bpSEXP, SEXP L_0SEXP, SEXP uSEXP, SEXP XSEXP, SEXP Y0SEXP, SEXP B0SEXP, SEXP open_sysSEXP, SEXP dtSEXP, SEXP max_tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type d_yp(d_ypSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type d_b0(d_b0SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type d_bp(d_bpSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type g_yp(g_ypSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type g_b0(g_b0SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type g_bp(g_bpSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type L_0(L_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const double& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type Y0(Y0SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< const bool& >::type open_sys(open_sysSEXP);
    Rcpp::traits::input_parameter< const double& >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< const double& >::type max_t(max_tSEXP);
    rcpp_result_gen = Rcpp::wrap(plant_metacomm_cpp(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, u, X, Y0, B0, open_sys, dt, max_t));
    return rcpp_result_gen;
END_RCPP
}
// dissimilarity
double dissimilarity(NumericVector yeast, NumericVector bact);
RcppExport SEXP _sweetsoursong_dissimilarity(SEXP yeastSEXP, SEXP bactSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type yeast(yeastSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bact(bactSEXP);
    rcpp_result_gen = Rcpp::wrap(dissimilarity(yeast, bact));
    return rcpp_result_gen;
END_RCPP
}
// dissimilarity_vector
NumericVector dissimilarity_vector(NumericVector yeast, NumericVector bact, const size_t& group_size, const bool& overall_mean);
RcppExport SEXP _sweetsoursong_dissimilarity_vector(SEXP yeastSEXP, SEXP bactSEXP, SEXP group_sizeSEXP, SEXP overall_meanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type yeast(yeastSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bact(bactSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type group_size(group_sizeSEXP);
    Rcpp::traits::input_parameter< const bool& >::type overall_mean(overall_meanSEXP);
    rcpp_result_gen = Rcpp::wrap(dissimilarity_vector(yeast, bact, group_size, overall_mean));
    return rcpp_result_gen;
END_RCPP
}
// diversity
double diversity(NumericVector yeast, NumericVector bact, double zero_threshold);
RcppExport SEXP _sweetsoursong_diversity(SEXP yeastSEXP, SEXP bactSEXP, SEXP zero_thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type yeast(yeastSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bact(bactSEXP);
    Rcpp::traits::input_parameter< double >::type zero_threshold(zero_thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(diversity(yeast, bact, zero_threshold));
    return rcpp_result_gen;
END_RCPP
}
// diversity_vector
double diversity_vector(NumericVector yeast, NumericVector bact, const size_t& group_size, double zero_threshold);
RcppExport SEXP _sweetsoursong_diversity_vector(SEXP yeastSEXP, SEXP bactSEXP, SEXP group_sizeSEXP, SEXP zero_thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type yeast(yeastSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bact(bactSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type group_size(group_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type zero_threshold(zero_thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(diversity_vector(yeast, bact, group_size, zero_threshold));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sweetsoursong_one_flower_ode", (DL_FUNC) &_sweetsoursong_one_flower_ode, 21},
    {"_sweetsoursong_plant_metacomm_stoch_cpp", (DL_FUNC) &_sweetsoursong_plant_metacomm_stoch_cpp, 21},
    {"_sweetsoursong_plant_metacomm_cpp", (DL_FUNC) &_sweetsoursong_plant_metacomm_cpp, 15},
    {"_sweetsoursong_dissimilarity", (DL_FUNC) &_sweetsoursong_dissimilarity, 2},
    {"_sweetsoursong_dissimilarity_vector", (DL_FUNC) &_sweetsoursong_dissimilarity_vector, 4},
    {"_sweetsoursong_diversity", (DL_FUNC) &_sweetsoursong_diversity, 3},
    {"_sweetsoursong_diversity_vector", (DL_FUNC) &_sweetsoursong_diversity_vector, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_sweetsoursong(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
