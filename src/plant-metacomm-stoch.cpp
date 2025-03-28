
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






// Generate from ~U(0,1) using a pcg32 object:
namespace pcg {
    const double max = static_cast<double>(pcg32::max());
}
inline double runif_01(pcg32& eng) {
    return (static_cast<double>(eng()) + 1) / (pcg::max + 2);
}



class StochLandscapeStepper {

    MatType det;
    MatType stoch;
    arma::vec zeta;
    size_t season_len;
    double season_surv;
    double q;
    size_t iters = 0;

public:

    typedef boost::numeric::odeint::stepper_tag stepper_category;

    static unsigned short order( void ) { return 1; }

    StochLandscapeStepper(const size_t& n_plants,
                          const size_t& n_states,
                          const size_t& season_len_,
                          const double& season_surv_,
                          const double& q_)
        : det(n_plants, n_states),
          stoch(n_plants, n_states),
          zeta(n_plants),
          season_len(season_len_),
          season_surv(season_surv_),
          q(q_) {}

    template< class System >
    void do_step(System& system, MatType& x, double t, double dt) {

        // No use continuing if x has NaNs or infinity
        if (x.has_nan() || x.has_inf()) return;

        // New season:
        if (iters >= season_len) {

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

            iters = 1; // this counts as the first day of this season

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

        iters++;

        return;
    }


};






// Stochastic process of the stochastic landscape
struct StochLandscapeStochProcess {

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
    size_t save_every;
    bool begin_end;
    int summarize;
    size_t season_len;
    double season_surv;
    double q;

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
                      const size_t& season_len_,
                      const double& season_surv_,
                      const double& q_,
                      const bool& open_sys,
                      const double& dt_,
                      const double& max_t_,
                      const double& burnin_,
                      const size_t& save_every_,
                      const bool& begin_end_,
                      const int& summarize_)
        : output(n_reps, MatType(0,0)),
          seeds(n_reps, std::vector<uint64_t>(2)),
          x0(m.size(), 2U),
          determ_sys0(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, u, X,
                      open_sys, max_t_, dt_, season_len_),
          n_sigma(n_sigma_),
          dt(dt_),
          max_t(max_t_),
          burnin(burnin_),
          save_every(save_every_),
          begin_end(begin_end_),
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
            do_work<MetaObsStoch>(begin, end);
        } else if (summarize == 1) {
            do_work<MetaObsStochSumm>(begin, end);
        } else {
            do_work<MetaObsStochSummRep>(begin, end);
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
     MetaObsStoch, MetaObsStochSumm, or MetaObsStochSummRep classes
     defined in `plant-metacomm.h`
     */
    template <class C>
    void do_work(const size_t& begin, const size_t& end) {
        pcg32 rng;
        const size_t& np(determ_sys0.n_plants);
        MatType x;
        C obs(burnin, save_every, season_len, begin_end);

        using System = std::pair<LandscapeConstF, StochLandscapeStochProcess>;

        for (size_t rep = begin; rep < end; rep++) {

            rng.seed(seeds[rep][0], seeds[rep][1]);
            StochLandscapeStochProcess stoch_sys(rng, n_sigma);

            x = x0;
            obs.clear();

            System system = std::make_pair(determ_sys0, stoch_sys);

            boost::numeric::odeint::integrate_const(
                StochLandscapeStepper(np, 2U, season_len, season_surv, q),
                system, x, 0.0, max_t, dt, std::ref(obs));

            obs.fill_output(output[rep], determ_sys0,
                            static_cast<double>(rep) + 1.0);

        }

        return;
    }
};




// [[Rcpp::export]]
arma::mat plant_metacomm_stoch_cpp(const size_t& n_reps,
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
                                   const size_t& season_len,
                                   const double& season_surv,
                                   const double& q,
                                   const bool& open_sys,
                                   const double& dt,
                                   const double& max_t,
                                   const double& burnin,
                                   const size_t& save_every,
                                   const bool& begin_end,
                                   const int& summarize) {

    /*
     I can't just use 'stop()' because it causes a segfault (or similar).
     The system below is my workaround.
     ^^ Holdover from using OpenMP. Remove later if you want. ^^
     */
    bool err = lanscape_constF_arg_checks(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                                          L_0, u, X, Y0, B0, dt, max_t);
    min_val_check(err, n_sigma, "n_sigma", 0, false);
    min_val_check(err, season_surv, "season_surv", 0);
    max_val_check(err, season_surv, "season_surv", 1);
    min_val_check(err, q, "q", 0);
    if (err) return arma::mat(0,0);


    StochLandCFWorker worker(n_reps, m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                             L_0, u, X, Y0, B0, n_sigma,
                             season_len, season_surv, q, open_sys,
                             dt, max_t, burnin, save_every, begin_end, summarize);

    RcppParallel::parallelFor(0, n_reps, worker);

    arma::mat output = worker.make_output();

    return output;
}
