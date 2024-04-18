#define _USE_MATH_DEFINES

#include <Rcpp.h>
#include <vector>
#include <cmath>

#include "ode.h"

using namespace Rcpp;


//[[Rcpp::export]]
NumericVector test_R(NumericVector time, const double& mu, const double& sigma) {
    NumericVector R(time.size());
    for (size_t i = 0; i < time.size(); i++) {
        R[i] = (1 / (sigma * std::sqrt(2 * M_PI))) *
            std::exp(-0.5 * std::pow((time[i] - mu) / sigma, 2U));
    }
    return R;
}

//[[Rcpp::export]]
NumericVector landscape_weights(const NumericMatrix& x,
                                const double& S_0,
                                const double& q,
                                const std::vector<double>& X,
                                const double& w,
                                const NumericMatrix& z) {

    NumericVector out;

    size_t n_plants = z.nrow();

    bool err = false;
    if (z.ncol() != n_plants) {
        Rcout << "z needs to be square!" << std::endl;
        err = true;
    }
    len_check<double>(err, X, "X", n_plants);
    if (x.nrow() != n_plants) {
        Rcout << "x needs to have " << std::to_string(n_plants) << " rows!" << std::endl;
        err = true;
    }
    if (x.ncol() != 3U) {
        Rcout << "x needs to have 3 cols!" << std::endl;
        err = true;
    }
    if (err) return out;


    std::vector<double> wts_vec(n_plants);

    MatType x_boost(x.nrow(), x.ncol());
    MatType exp_wz(z.nrow(), z.ncol());

    for (size_t i = 0; i < n_plants; i++) {
        // fill x_boost:
        for (size_t k = 0; k < 3U; k++) x_boost(i,k) = x(i,k);
        // also fill exp_wz:
        for (size_t j = 0; j < n_plants; j++) {
            if (i == j) {
                exp_wz(i,j) = 1;
            } else {
                exp_wz(i,j) = std::exp(-w * z(i, j));
            }
        }
    }

    landscape_weights__(wts_vec, x_boost, n_plants, S_0, q, X, exp_wz);

    out = wrap(wts_vec);

    return out;

}
