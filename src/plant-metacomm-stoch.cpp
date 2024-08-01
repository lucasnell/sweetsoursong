
/*
 Lanscape of plants where the total number of flowers (F) at each plant is
 constant and where pollinators are not affected by flower densities.
 We don't even include F and instead model Y and B as proportions
 of total flowers with non-colonized N = 1 - Y - B.
 We also add stochasticity.
 */

#define _USE_MATH_DEFINES


#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <random>

#include "ode.h"
#include "plant-metacomm.h"

#include <RcppParallel.h>
#include <pcg_random.hpp>


using namespace Rcpp;


// this deals with std::remainder's rounding issues
inline bool zero_remainder(const double& numer, const double& denom) {
    return std::abs(std::remainder(numer, denom)) < 1e-10;
}

// logit and inverse logit functions
inline void logit(const double& p, double& x) {
    x = std::log(p / (1 - p));
    return;
}
inline void inv_logit(const double& x, double& p) {
    p = 1 / (1 + std::exp(- x));
    return;
}
// overloaded for changing doubles in place
inline void logit(double& x) {
    x = std::log(x / (1 - x));
    return;
}
inline void inv_logit(double& p) {
    p = 1 / (1 + std::exp(- p));
    return;
}



class StochLandscapeStepper
{
public:

    typedef boost::numeric::odeint::stepper_tag stepper_category;

    static unsigned short order( void ) { return 1; }

    StochLandscapeStepper(const size_t& n_plants,
                          const size_t& n_states,
                          const double& season_len_,
                          const double& season_surv_,
                          const double& season_sigma_)
        : det(n_plants, n_states),
          stoch(n_plants, n_states),
          season_len(season_len_),
          season_surv(season_surv_),
          season_sigma(season_sigma_) {}

    template< class System >
    void do_step(System system, MatType& x, double t, double dt) {
        // New season:
        if (t > 0 && zero_remainder(t, season_len)) {
            // If no seasonal variation included:
            if (season_sigma <= 0) {

                x *= season_surv;

            } else {

                // If seasonal variation included:

                x.elem( arma::find(x < 1) ) += 1e-6; // to remove zeros

                pcg32& rng(system.second.m_rng);
                std::normal_distribution<double>& norm(system.second.m_dist);
                double YB, xij;
                for (size_t i = 0 ; i < x.n_rows ; i++) {
                    YB = arma::accu(x.row(i)); // Y+B for this plant
                    for (size_t j = 0 ; j < x.n_cols ; j++) {
                        xij = x(i,j) / YB;
                        logit(xij);
                        xij += (season_sigma * norm(rng));
                        inv_logit(xij);
                        x(i,j) = xij;
                    }
                    // This makes sure it always sums to (YB * season_surv):
                    x.row(i) /= arma::accu(x.row(i));
                    x.row(i) *= (YB * season_surv);
                }
            }
            return;
        }
        // Standard iteration:
        system.first(x, det);
        system.second(x, stoch);
        double sqrt_dt = std::sqrt(dt);
        for (size_t i = 0 ; i < x.n_rows ; i++) {
            for (size_t j = 0 ; j < x.n_cols ; j++) {
                x(i,j) += dt * det(i,j) + sqrt_dt * stoch(i,j);
                if (x(i,j) > 1) x(i,j) = 1;
                if (x(i,j) < 0) x(i,j) = 0;
            }
        }
        return;
    }

private:
    MatType det;
    MatType stoch;
    double season_len;
    double season_surv;
    double season_sigma;
};






// Stochastic process of the stochastic landscape
struct StochLandscapeStochProcess
{
    pcg32& m_rng;
    std::normal_distribution<double> m_dist;
    double n_sigma;

    StochLandscapeStochProcess(pcg32& rng,
                               double n_sigma_)
        : m_rng(rng),
          m_dist(0.0, 1.0),
          n_sigma(n_sigma_) {}

    void operator()(const MatType &x, MatType &dxdt) {

        double stdev;

        for (size_t i = 0 ; i < x.n_rows ; i++) {
            for (size_t j = 0 ; j < x.n_cols ; j++) {
                stdev = std::sqrt(x(i,j) * (1 - x(i,j)) / n_sigma);
                dxdt(i,j) = stdev * m_dist(m_rng);
            }
        }

        return;
    }
};



// RcppParallel Worker to do runs for a single thread:
struct StochLandCFWorker : public RcppParallel::Worker {

    std::vector<MatType> output;
    std::vector<std::vector<uint64_t>> seeds;
    MatType x0;
    LandscapeConstF determ_sys0;
    double n_sigma;
    double dt;
    double max_t;
    double season_len;
    double season_surv;
    double season_sigma;

    StochLandCFWorker(const uint32_t& n_reps,
                      const std::vector<double>& m,
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
                      const double& n_sigma_,
                      const double& season_len_,
                      const double& season_surv_,
                      const double& season_sigma_,
                      const double& dt_,
                      const double& max_t_)
        : output(n_reps, MatType(0,0)),
          seeds(n_reps, std::vector<uint64_t>(2)),
          x0(m.size(), 2U),
          determ_sys0(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, u, X),
          n_sigma(n_sigma_),
          dt(dt_),
          max_t(max_t_),
          season_len(season_len_),
          season_surv(season_surv_),
          season_sigma(season_sigma_) {

        std::vector<uint64_t> tmp_seeds(4);
        for (uint32_t i = 0; i < n_reps; i++) {
            // These are 32-bit integers cast as 64-bit for downstream compatibility
            tmp_seeds = as<std::vector<uint64_t>>(Rcpp::runif(4,0,4294967296));
            // Converting to two 64-bit seeds for pcg32
            seeds[i][0] = (tmp_seeds[0]<<32) + tmp_seeds[1];
            seeds[i][1] = (tmp_seeds[2]<<32) + tmp_seeds[3];
        }

        for (size_t i = 0; i < x0.n_rows; i++) {
            x0(i,0) = Y0[i];
            x0(i,1) = B0[i];
        }

    };

    void operator()(size_t begin, size_t end) {

        pcg32 rng;
        const size_t& np(determ_sys0.n_plants);
        MatType x;
        Observer<MatType> obs;
        std::vector<double> wts(np);

        for (size_t rep = begin; rep < end; rep++) {

            rng.seed(seeds[rep][0], seeds[rep][1]);

            x = x0;
            obs.data.clear();
            obs.time.clear();

            boost::numeric::odeint::integrate_const(
                StochLandscapeStepper(np, 2U, season_len, season_surv, season_sigma),
                std::make_pair(determ_sys0,
                               StochLandscapeStochProcess(rng, n_sigma)),
                               x, 0.0, max_t, dt, std::ref(obs));

            size_t n_steps = obs.data.size();
            output[rep].set_size(n_steps * np, 6U);
            // colnames(output[rep]) = CharacterVector::create("rep", "t", "p", "Y", "B", "P");
            size_t i = 0;
            double dbl_rep = static_cast<double>(rep) + 1;
            for (size_t t = 0; t < n_steps; t++) {
                determ_sys0.make_weights(wts, obs.data[t]);
                for (size_t k = 0; k < np; k++) {
                    output[rep](i,0) = dbl_rep;
                    output[rep](i,1) = obs.time[t];
                    output[rep](i,2) = k;
                    output[rep](i,3) = obs.data[t](k,0);
                    output[rep](i,4) = obs.data[t](k,1);
                    output[rep](i,5) = wts[k];
                    i++;
                }

            }

        }
        return;
    }
};




//' @export
// [[Rcpp::export]]
NumericMatrix plant_metacomm_stoch(const uint32_t& n_reps,
                                   const std::vector<double>& m,
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
                                   const double& n_sigma,
                                   SEXP season_len = R_NilValue,
                                   const double& season_surv = 0.01,
                                   const double& season_sigma = 0,
                                   const double& dt = 0.1,
                                   const double& max_t = 100.0) {

    /*
     I can't just use 'stop()' because it causes a segfault (or similar).
     The system below is my workaround.
     */
    bool err = lanscape_constF_arg_checks(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                                          L_0, u, X, Y0, B0, dt, max_t);
    min_val_check(err, n_sigma, "n_sigma", 0, false);
    double season_len_ = (season_len == R_NilValue) ? max_t + 1.0 : as<double>(season_len);
    min_val_check(err, season_len_, "season_len", 0, false);
    min_val_check(err, season_surv, "season_surv", 0);
    max_val_check(err, season_surv, "season_surv", 1);
    if (season_len_ < max_t && ! zero_remainder(season_len_, dt)) {
        Rcout << "season_len is " << std::to_string(season_len_);
        Rcout << " but should be divisible by dt (";
        Rcout << std::to_string(dt) << ")!" << std::endl;
        err = true;
    }
    min_val_check(err, season_sigma, "season_sigma", 0);
    if (err) return NumericMatrix(0,0);


    StochLandCFWorker worker(n_reps, m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                             L_0, u, X, Y0, B0, n_sigma,
                             season_len_, season_surv, season_sigma,
                             dt, max_t);

    RcppParallel::parallelFor(0, n_reps, worker);

    size_t n_rows = 0;
    for (const MatType& m : worker.output) n_rows += m.n_rows;

    NumericMatrix output(n_rows, 6U);
    colnames(output) = CharacterVector::create("rep", "t", "p", "Y", "B", "P");
    size_t i = 0;
    for (size_t rep = 0; rep < worker.output.size(); rep++) {
        const MatType& m(worker.output[rep]);
        for (size_t k = 0; k < m.n_rows; k++) {
            for (size_t j = 0; j < m.n_cols; j++) {
                output(i,j) = m(k,j);
            }
            i++;
        }

    }
    return output;
}
