
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
inline double logit(const double& p) {
    double x = std::log(p / (1 - p));
    return x;
}
inline double inv_logit(const double& x) {
    double p = 1 / (1 + std::exp(- x));
    return p;
}
// overloaded for creating new vectors:
inline arma::vec logit(const arma::vec& p) {
    arma::vec x(arma::size(p));
    for (size_t i = 0; i < p.n_elem; i++) x(i) = std::log(p(i) / (1 - p(i)));
    return x;
}
inline arma::vec inv_logit(const arma::vec& x) {
    arma::vec p(arma::size(x));
    for (size_t i = 0; i < x.n_elem; i++) p(i) = 1 / (1 + std::exp(- x(i)));
    return p;
}
// overloaded for filling vectors:
inline void logit(const arma::vec& p, arma::vec& x) {
    if (arma::size(x) != arma::size(p)) x.set_size(p.n_elem);
    for (size_t i = 0; i < p.n_elem; i++) x(i) = std::log(p(i) / (1 - p(i)));
    return;
}
inline void inv_logit(const arma::vec& x, arma::vec& p) {
    if (arma::size(p) != arma::size(x)) p.set_size(x.n_elem);
    for (size_t i = 0; i < x.n_elem; i++) p(i) = 1 / (1 + std::exp(- x(i)));
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
                          const double& q_)
        : det(n_plants, n_states),
          stoch(n_plants, n_states),
          zeta(n_plants),
          season_len(season_len_),
          season_surv(season_surv_),
          q(q_) {}

    template< class System >
    void do_step(System system, MatType& x, double t, double dt) {

        // No use continuing if x has NaNs or infinity
        if (x.has_nan() || x.has_inf()) return;

        // New season:
        if (t > 0 && zero_remainder(t, season_len)) {

            if (q >= 1) {

                // If there is no between-season noise, don't bother generating
                // random numbers:

                x *= season_surv;

            } else {

                pcg32& rng(system.second.m_rng);

                // This is to avoid doing this multiple times
                arma::rowvec x_sum = arma::sum(x);

                for (size_t j = 0; j < x.n_cols; j++) {

                    // to avoid processing columns with all zeros bc they
                    // don't need it:
                    if (x_sum(j) <= 0) continue;

                    // Vector of ~U(0,1) normalized to have the same sum
                    // as x.col(j)
                    double z_sum = 0;
                    for (double& z : zeta) {
                        z = runif_01(rng);
                        z_sum += z;
                    }
                    double z_mult = x_sum(j) / z_sum;
                    for (double& z : zeta) z *= z_mult;

                    double qxz;
                    for (size_t i = 0; i < zeta.n_rows; i++) {
                        qxz = q * x(i,j) + (1 - q) * zeta(i);
                        x(i,j) = season_surv * qxz;
                    }

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
    arma::vec zeta;
    double season_len;
    double season_surv;
    double q;
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
    int summarize;
    double season_len;
    double season_surv;
    double q;
    int status = 0;

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
                      const double& q_,
                      const bool& open_sys,
                      const double& dt_,
                      const double& max_t_,
                      const double& burnin_,
                      const int& summarize_)
        : output(n_reps, MatType(0,0)),
          seeds(n_reps, std::vector<uint64_t>(2)),
          x0(m.size(), 2U),
          determ_sys0(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, u, X,
                      open_sys, max_t_, season_len_),
          n_sigma(n_sigma_),
          dt(dt_),
          max_t(max_t_),
          burnin(burnin_),
          summarize(summarize_),
          season_len(season_len_),
          season_surv(season_surv_),
          q(q_) {

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

        if (summarize == 0) {
            status = do_work<MetaObsStoch>(begin, end);
        } else if (summarize == 1) {
            status = do_work<MetaObsStochSumm>(begin, end);
        } else {
            status = do_work<MetaObsStochSummRep>(begin, end);
        }

        return;
    }

    // make final output matrix and move `output` items there (removing them
    // from this object) as you go
    MatType make_output() {

        size_t n_rows = 0;
        size_t n_cols = this->output.front().n_cols;
        for (const MatType& m : this->output) {
            n_rows += m.n_rows;
            if (m.n_cols > n_cols) n_cols = m.n_cols;
        }

        MatType final_out(n_rows, n_cols, arma::fill::none);
        size_t i = 0;
        for (size_t rep = 0; rep < this->output.size(); rep++) {
            MatType& m(this->output[rep]);
            for (size_t k = 0; k < m.n_rows; k++) {
                for (size_t j = 0; j < m.n_cols; j++) {
                    final_out(i,j) = m(k,j);
                }
                i++;
            }
            m.reset();

        }
        return final_out;
    }

private:

    /*
     Main function that does most of the work, can be used for the
     MetaObsStoch, MetaObsStochSumm, or MetaObsStochSummRep classed
     defined in `plant-metacomm.h`
     */
    template <class C>
    int do_work(const size_t& begin, const size_t& end) {
        pcg32 rng;
        const size_t& np(determ_sys0.n_plants);
        MatType x;
        C obs(burnin);

        int status = 0;


        for (size_t rep = begin; rep < end; rep++) {

            rng.seed(seeds[rep][0], seeds[rep][1]);

            x = x0;
            obs.clear();

            boost::numeric::odeint::integrate_const(
                StochLandscapeStepper(np, 2U, season_len, season_surv, q),
                std::make_pair(determ_sys0,
                               StochLandscapeStochProcess(rng, n_sigma)),
                               x, 0.0, max_t, dt, std::ref(obs));

            obs.fill_output(output[rep], determ_sys0,
                            static_cast<double>(rep) + 1.0);

            if (x.has_inf()) {
                status = 1;
                break;
            }
        }

        return status;
    }
};




// [[Rcpp::export]]
arma::mat plant_metacomm_stoch_cpp(const uint32_t& n_reps,
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
                                   const double& q,
                                   const bool& open_sys,
                                   const double& dt,
                                   const double& max_t,
                                   const double& burnin,
                                   const int& summarize) {

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
    min_val_check(err, q, "q", 0);
    if (season_len < max_t && ! zero_remainder(season_len, dt)) {
        Rcout << "season_len is " << std::to_string(season_len);
        Rcout << " but should be divisible by dt (";
        Rcout << std::to_string(dt) << ")!" << std::endl;
        err = true;
    }
    if (err) return arma::mat(0,0);


    StochLandCFWorker worker(n_reps, m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                             L_0, u, X, Y0, B0, n_sigma,
                             season_len, season_surv, q, open_sys,
                             dt, max_t, burnin, summarize);

    RcppParallel::parallelFor(0, n_reps, worker);

    if (worker.status == 1) {
        stop("Value of c too high in StochLandscapeStepper::do_step function");
    }

    arma::mat output = worker.make_output();

    return output;
}
