
/*
 Lanscape of plants where the total number of flowers (F) at each plant is
 constant and where pollinators are not affected by flower densities.
 We don't even include F and instead model Y and B as proportions
 of total flowers with non-colonized N = 1 - Y - B.
 */

#include <RcppArmadillo.h>
#include <vector>

#include "ode.h"
#include "plant-metacomm.h"

using namespace Rcpp;




// [[Rcpp::export]]
arma::mat plant_metacomm_cpp(const std::vector<double>& m,
                             const std::vector<double>& d_yp,
                             const std::vector<double>& d_b0,
                             const std::vector<double>& d_bp,
                             const std::vector<double>& g_yp,
                             const std::vector<double>& g_b0,
                             const std::vector<double>& g_bp,
                             const std::vector<double>& L_0,
                             const double& u,
                             const double& X,
                             const std::vector<double>& Y0,
                             const std::vector<double>& B0,
                             const bool& open_sys,
                             const double& dt,
                             const double& max_t) {

    size_t np = m.size();
    /*
     I can't just use 'stop()' because it causes a segfault (or similar).
     The system below is my workaround.
     */
    bool err = lanscape_constF_arg_checks(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                                          L_0, u, X, Y0, B0, dt, max_t);
    if (err) return arma::mat(0,0);

    size_t n_states = 2U;
    MatType x(np, n_states);
    for (size_t i = 0; i < np; i++) {
        x(i,0) = Y0[i];
        x(i,1) = B0[i];
    }


    MetaObs obs;
    LandscapeConstF system(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, u, X,
                           open_sys, max_t, max_t + 1.0);

    boost::numeric::odeint::integrate_const(
        MatStepperType(), std::ref(system),
        x, 0.0, max_t, dt, std::ref(obs));

    MatType output;
    obs.fill_output(output, system);

    return output;
}
