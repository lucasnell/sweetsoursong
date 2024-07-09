#define _USE_MATH_DEFINES

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>

#include "ode.h"

using namespace Rcpp;



//[[Rcpp::export]]
NumericMatrix make_dist_mat_rcpp(const NumericVector& x,
                                 const NumericVector& y) {

    size_t n = x.size();
    NumericMatrix dm(n, n);

    double xdiff, ydiff;

    for (size_t i = 0; i < (n-1); i++) {
        for (size_t j = i+1; j < n; j++) {
            xdiff = x(i) - x(j);
            ydiff = y(i) - y(j);
            dm(j,i) = std::sqrt(xdiff * xdiff + ydiff * ydiff);
            dm(i,j) = dm(j,i);
        }
    }

    return dm;

}


//[[Rcpp::export]]
NumericMatrix make_spat_wts_rcpp(const NumericMatrix& dm, const double& m) {

    size_t n = dm.nrow();

    NumericMatrix sw(n, n);

    for (size_t i = 0; i < (n-1); i++) {
        for (size_t j = i+1; j < n; j++) {
            sw(j,i) = 1 / std::pow(dm(i,j), m);
            sw(i,j) = sw(j,i);
        }
    }

    return sw;

}



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

    Rcout << "this function is outdated. Returning nothing." << std::endl;

    // size_t n_plants = z.nrow();
    //
    // bool err = false;
    // if (z.ncol() != n_plants) {
    //     Rcout << "z needs to be square!" << std::endl;
    //     err = true;
    // }
    // len_check(err, X, "X", n_plants);
    // if (x.nrow() != n_plants) {
    //     Rcout << "x needs to have " << std::to_string(n_plants) << " rows!" << std::endl;
    //     err = true;
    // }
    // if (x.ncol() != 3U) {
    //     Rcout << "x needs to have 3 cols!" << std::endl;
    //     err = true;
    // }
    // if (err) return out;
    //
    //
    // std::vector<double> wts_vec(n_plants);
    //
    // MatType x_boost(x.nrow(), x.ncol());
    // MatType exp_wz(z.nrow(), z.ncol());
    //
    // for (size_t i = 0; i < n_plants; i++) {
    //     // fill x_boost:
    //     for (size_t k = 0; k < 3U; k++) x_boost(i,k) = x(i,k);
    //     // also fill exp_wz:
    //     for (size_t j = 0; j < n_plants; j++) {
    //         if (i == j) {
    //             exp_wz(i,j) = 1;
    //         } else {
    //             exp_wz(i,j) = std::exp(-w * z(i, j));
    //         }
    //     }
    // }
    //
    // landscape_weights__(wts_vec, x_boost, n_plants, S_0, q, X, exp_wz);
    //
    // out = wrap(wts_vec);

    return out;

}




/*
 =====================================================================================
 =====================================================================================
 Community metrics
 =====================================================================================
 =====================================================================================
 */

//' Bray–Curtis dissimilarity.
//'
//' @param yeast Vector of yeast abundances.
//' @param bact Vector of bacteria abundances.
//'
//' @name dissimilarity
//'
//' @return A single number indicating the mean dissimilarity across the
//'     two vectors.
//'
//' @export
//'
//[[Rcpp::export]]
double dissimilarity(NumericVector yeast, NumericVector bact) {
    size_t n = yeast.size();
    if (n != bact.size()) stop("lengths do not match");
    size_t n_combos = static_cast<size_t>(Rf_choose(static_cast<double>(n), 2.0));
    NumericVector bc_vec(n_combos);
    size_t k = 0;
    double min_y, min_b, denom;
    for (size_t i = 0; i < (n-1U); i++) {
        for (size_t j = i+1U; j < n; j++) {
            min_y = (yeast(i) < yeast(j)) ? yeast(i) : yeast(j);
            min_b = (bact(i) < bact(j)) ? bact(i) : bact(j);
            denom = yeast(i) + yeast(j) + bact(i) + bact(j);
            bc_vec(k) = 1 - (2 * (min_y + min_b)) / denom;
            k++;
        }
    }
    double bc_mean = mean(bc_vec);
    return bc_mean;
}


//' Bray–Curtis dissimilarity on grouped vectors.
//'
//' @inheritParams dissimilarity
//' @param group_size Size of groups within which dissimilarities should be
//'   calculated. It's assumed that the vectors are
//'   all sorted by this grouping such that every `group_size` elements
//'   in the vectors belong to the same group.
//'   The `yeast` and `bact` vectors' lengths should be divisible by
//'   `group_size`.
//' @param overall_mean Single logical for whether to return the overall mean
//'   dissimilarity after grouping. Defaults to `FALSE`.
//'
//'
//' @return If `overall_mean` is `FALSE`, a numeric vector of length
//'   `length(yeast) %/% group_size` with the mean dissimilarities across
//'   the vectors for each group.
//'   Otherwise, a single number indicating the mean dissimilarity across the
//'   two vectors after grouping.
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector dissimilarity_vector(NumericVector yeast,
                                   NumericVector bact,
                                   const size_t& group_size,
                                   const bool& overall_mean = false) {

    if (yeast.size() != bact.size()) stop("lengths do not match");
    if (yeast.size() % group_size != 0) stop("vector lengths aren't divisible by group_size");

    size_t total_groups = yeast.size() / group_size;
    NumericVector bc_mean_vec;
    if (overall_mean) {
        bc_mean_vec = NumericVector(1);
    } else bc_mean_vec = static_cast<NumericVector>(no_init(total_groups));

    const size_t& n(group_size);
    double n_combos = Rf_choose(static_cast<double>(n), 2.0);
    double bc_sum, min_y, min_b, denom;
    size_t i0;

    for (size_t k = 0; k < total_groups; k++) {
        bc_sum = 0;
        i0 = k * group_size;
        for (size_t i = i0; i < (i0+n-1U); i++) {
            for (size_t j = i+1U; j < (i0+n); j++) {
                if (j >= yeast.size()) stop("vector lengths aren't divisible by group_size!!");
                min_y = (yeast(i) < yeast(j)) ? yeast(i) : yeast(j);
                min_b = (bact(i) < bact(j)) ? bact(i) : bact(j);
                denom = yeast(i) + yeast(j) + bact(i) + bact(j);
                bc_sum += 1 - (2 * (min_y + min_b)) / denom;
            }
        }
        if (overall_mean) {
            bc_mean_vec(0) += bc_sum;
        } else bc_mean_vec(k) = bc_sum / n_combos;
    }

    if (overall_mean) bc_mean_vec(0) /= (n_combos * total_groups);

    return bc_mean_vec;
}






//' Shannon diversity index.
//'
//' @param yeast Vector of yeast abundances.
//' @param bact Vector of bacteria abundances.
//'
//' @name diversity
//'
//' @return A single number indicating the mean diversity across the
//'     two vectors.
//'
//' @export
//'
//[[Rcpp::export]]
double diversity(NumericVector yeast, NumericVector bact,
                 double zero_threshold = 2.220446e-16) {
    size_t n = yeast.size();
    if (n != bact.size()) stop("lengths do not match");
    // Do calculation while accounting for zeros:
    double p_yeast, p_bact, total, H_i, H_mean = 0;
    for (size_t i = 0; i < n; i++) {
        if (bact(i) > zero_threshold && yeast(i) > zero_threshold) {
            total = yeast(i) + bact(i);
            p_yeast = yeast(i) / total;
            p_bact = bact(i) / total;
            H_i = - p_yeast * log(p_yeast) - p_bact * log(p_bact);
            H_mean += H_i;
        }
    }
    H_mean /= static_cast<double>(n);
    return H_mean;
}
