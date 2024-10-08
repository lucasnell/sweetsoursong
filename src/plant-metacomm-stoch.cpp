
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

// Generate from ~U(0,1) using a pcg32 object:
namespace pcg {
    const double max = static_cast<double>(pcg32::max());
}
inline double runif_01(pcg32& eng) {
    return (static_cast<double>(eng()) + 1) / (pcg::max + 2);
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
                          const int& rand_season_)
        : det(n_plants, n_states),
          stoch(n_plants, n_states),
          season_rands(n_plants, n_states),
          season_len(season_len_),
          season_surv(season_surv_),
          rand_season(rand_season_) {}

    template< class System >
    void do_step(System system, MatType& x, double t, double dt) {
        // New season:
        if (t > 0 && zero_remainder(t, season_len)) {

            if (rand_season == 0) {

                // If starting abundances for each flower are based on
                // that flower's abundances from the previous season:

                x *= season_surv;

            } else if (rand_season == 1) {

                // If each flower's starting abundances are NOT based on
                // previous season, but still result in a CLOSED system
                // (i.e., the landscape-level abundances do NOT change):

                double n_plants = static_cast<double>(x.n_rows);

                pcg32& rng(system.second.m_rng);
                // Fill `season_rands` with ~ U(0,1), then make each column
                // sum to `n`:
                for (size_t j = 0 ; j < season_rands.n_cols; j++) {
                    double rand_sum = 0;
                    for (double& rnd : season_rands.col(j)) {
                        rnd = runif_01(rng);
                        rand_sum += rnd;
                    }
                    season_rands.col(j) *= (n_plants / rand_sum);
                }
                // Average abundance of each microbe type:
                double Ybar = arma::mean(x.col(0));
                double Bbar = arma::mean(x.col(1));
                // Now fill in new season's starting abundances:
                for (size_t i = 0 ; i < x.n_rows ; i++) {
                    x(i,0) = season_surv * Ybar * season_rands(i,0);
                    x(i,1) = season_surv * Bbar * season_rands(i,1);
                }

            } else {

                // If starting abundances are NOT based on previous season,
                // but result in an OPEN system
                // (i.e., the landscape-level abundances do change):

                pcg32& rng(system.second.m_rng);
                double YB, rnd_u;
                for (size_t i = 0 ; i < x.n_rows ; i++) {
                    YB = arma::accu(x.row(i)); // Y+B for this plant
                    rnd_u = runif_01(rng);
                    x(i,0) = season_surv * rnd_u * YB;
                    x(i,1) = season_surv * (1 - rnd_u) * YB;
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
    MatType season_rands;
    double season_len;
    double season_surv;
    int rand_season;
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
    double burnin;
    double season_len;
    double season_surv;
    int rand_season;

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
                      const int& rand_season_,
                      const bool& open_sys,
                      const double& dt_,
                      const double& max_t_,
                      const double& burnin_)
        : output(n_reps, MatType(0,0)),
          seeds(n_reps, std::vector<uint64_t>(2)),
          x0(m.size(), 2U),
          determ_sys0(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, u, X, open_sys),
          n_sigma(n_sigma_),
          dt(dt_),
          max_t(max_t_),
          burnin(burnin_),
          season_len(season_len_),
          season_surv(season_surv_),
          rand_season(rand_season_) {

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
        ObserverBI<MatType> obs(burnin);
        std::vector<double> wts(np);

        for (size_t rep = begin; rep < end; rep++) {

            rng.seed(seeds[rep][0], seeds[rep][1]);

            x = x0;
            obs.data.clear();
            obs.time.clear();

            boost::numeric::odeint::integrate_const(
                StochLandscapeStepper(np, 2U, season_len, season_surv, rand_season),
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




// [[Rcpp::export]]
NumericMatrix plant_metacomm_stoch_cpp(const uint32_t& n_reps,
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
                                       const double& season_len,
                                       const double& season_surv,
                                       const int& rand_season,
                                       const bool& open_sys,
                                       const double& dt,
                                       const double& max_t,
                                       const double& burnin) {

    /*
     I can't just use 'stop()' because it causes a segfault (or similar).
     The system below is my workaround.
     */
    bool err = lanscape_constF_arg_checks(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                                          L_0, u, X, Y0, B0, dt, max_t);
    min_val_check(err, n_sigma, "n_sigma", 0, false);
    min_val_check(err, season_len, "season_len", 0, false);
    min_val_check(err, season_surv, "season_surv", 0);
    max_val_check(err, season_surv, "season_surv", 1);
    if (season_len < max_t && ! zero_remainder(season_len, dt)) {
        Rcout << "season_len is " << std::to_string(season_len);
        Rcout << " but should be divisible by dt (";
        Rcout << std::to_string(dt) << ")!" << std::endl;
        err = true;
    }
    if (err) return NumericMatrix(0,0);


    StochLandCFWorker worker(n_reps, m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                             L_0, u, X, Y0, B0, n_sigma,
                             season_len, season_surv, rand_season, open_sys,
                             dt, max_t, burnin);

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
