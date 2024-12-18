#define _USE_MATH_DEFINES

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>

#include "ode.h"

using namespace Rcpp;




// logit and inverse logit functions

//' Logit function
//'
//' @param p Numeric vector of values in range `[0, 1]`.
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector logit(NumericVector p) {
    NumericVector x(p.size());
    for (size_t i = 0; i < p.size(); i++) {
        x(i) = std::log(p(i) / (1 - p(i)));
    }
    return x;
}


//' Inverse logit function
//'
//' @param x Numeric vector.
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector inv_logit(NumericVector x) {
    NumericVector p(x.size());
    for (size_t i = 0; i < p.size(); i++) {
        p(i) = 1 / (1 + std::exp(- x(i)));
    }
    return p;
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
    double n_combos = 0; // only count combos where both patches aren't empty
    double bc_sum = 0;
    double min_y, min_b, denom;
    for (size_t i = 0; i < (n-1U); i++) {
        for (size_t j = i+1U; j < n; j++) {
            denom = yeast(i) + yeast(j) + bact(i) + bact(j);
            if (denom > 0) {
                min_y = (yeast(i) < yeast(j)) ? yeast(i) : yeast(j);
                min_b = (bact(i) < bact(j)) ? bact(i) : bact(j);
                bc_sum += (1 - (2 * (min_y + min_b)) / denom);
                n_combos += 1.0;
            }
        }
    }
    double bc_mean = bc_sum / n_combos;
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
    double n_combos, bc_sum, min_y, min_b, denom;
    size_t i0;

    for (size_t k = 0; k < total_groups; k++) {
        n_combos = 0; // only count combos where both patches aren't empty
        bc_sum = 0;
        i0 = k * group_size;
        for (size_t i = i0; i < (i0+n-1U); i++) {
            for (size_t j = i+1U; j < (i0+n); j++) {
                if (j >= yeast.size()) stop("vector lengths aren't divisible by group_size!!");
                denom = yeast(i) + yeast(j) + bact(i) + bact(j);
                if (denom > 0) {
                    min_y = (yeast(i) < yeast(j)) ? yeast(i) : yeast(j);
                    min_b = (bact(i) < bact(j)) ? bact(i) : bact(j);
                    bc_sum += (1 - (2 * (min_y + min_b)) / denom);
                    n_combos += 1.0;
                }
            }
        }
        if (overall_mean) {
            bc_mean_vec(0) += (bc_sum / n_combos);
        } else bc_mean_vec(k) = bc_sum / n_combos;
    }

    if (overall_mean) bc_mean_vec(0) /= total_groups;

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


//' Shannon diversity index after summing by groups.
//'
//'
//' @inheritParams diversity
//' @param group_size Size of groups within which abundances should be
//'   summed before diversities are calculated.
//'   It's assumed that the vectors are all sorted by this grouping such
//'   that every `group_size` elements in the vectors belong to the same group.
//'   The `yeast` and `bact` vectors' lengths should be divisible by
//'   `group_size`.
//'
//'
//' @return A single number indicating the mean diversities across the
//'   two vectors after summing species by groups.
//'
//' @export
//'
//[[Rcpp::export]]
double diversity_vector(NumericVector yeast,
                        NumericVector bact,
                        const size_t& group_size,
                        double zero_threshold = 2.220446e-16) {

    if (yeast.size() != bact.size()) stop("lengths do not match");
    if (yeast.size() % group_size != 0) stop("vector lengths aren't divisible by group_size");

    size_t total_groups = yeast.size() / group_size;

    double p_yeast, p_bact, total, H_i, y, b;
    double H_mean = 0;

    const size_t& n(group_size);

    for (size_t i0 = 0; i0 < yeast.size(); i0+=group_size) {
        // First sum within group:
        y = 0;
        b = 0;
        for (size_t i = i0; i < (i0+n); i++) {
            y += yeast(i);
            b += bact(i);
        }
        // Now calculate H:
        if (b > zero_threshold && y > zero_threshold) {
            total = y + b;
            p_yeast = y / total;
            p_bact = b / total;
            H_i = - p_yeast * log(p_yeast) - p_bact * log(p_bact);
            H_mean += H_i;
        }
    }
    H_mean /= static_cast<double>(total_groups);

    return H_mean;
}
