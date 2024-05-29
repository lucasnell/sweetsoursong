
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

#include "ode.h"
#include "landscape_constantF.h"

using namespace Rcpp;



class StochLandscapeStepper
{
public:

    typedef boost::numeric::odeint::stepper_tag stepper_category;

    static unsigned short order( void ) { return 1; }

    StochLandscapeStepper(const size_t& n_plants,
                          const size_t& n_states)
        : det(n_plants, n_states),
          stoch(n_plants, n_states) {}

    template< class System >
    void do_step(System system, MatType& x, double t, double dt) {
        system.first(x, det);
        system.second(x, stoch);
        double sqrt_dt = std::sqrt(dt);
        for (size_t i = 0 ; i < x.n_rows ; i++) {
            for (size_t j = 0 ; j < x.n_cols ; j++) {
                x(i,j) += dt * det(i,j) + sqrt_dt * stoch(i,j);
                // if (x(i,j) > 1) x(i,j) = 1;
                // if (x(i,j) < 0) x(i,j) = 0;
            }
        }
        return;
    }

private:
    MatType det;
    MatType stoch;
};




// struct LandscapeConstF
// {
//
//     double mu;
//
//     LandscapeConstF(const double& mu_) : mu(mu_) {}
//
//     void operator()(const MatType &x, MatType &dxdt) const {
//         double tau = 1;
//         dxdt[0] = -(x[0] - mu) / tau;
//         return;
//     }
// };



// Stochastic process of the stochastic landscape
struct StochLandscapeStochProcess
{
    boost::mt11213b &m_rng;
    boost::normal_distribution<> m_dist;
    std::vector<double> sigma;

    StochLandscapeStochProcess(boost::mt11213b &rng,
                               double sigma_y_,
                               double sigma_b_)
        : m_rng(rng),
          m_dist(0.0, 1.0),
          sigma({sigma_y_, sigma_b_}) {}
    StochLandscapeStochProcess(boost::mt11213b &rng,
                               std::vector<double> sigma_)
        : m_rng(rng),
          m_dist(0.0, 1.0),
          sigma(sigma_) {}

    void operator()(const MatType &x, MatType &dxdt) {

        double stdev;

        for (size_t i = 0 ; i < x.n_rows ; i++) {
            for (size_t j = 0 ; j < x.n_cols ; j++) {
                stdev = 4 * x(i,j) * (1 - x(i,j));
                // stdev = 1;
                dxdt(i,j) = stdev * sigma[j] * m_dist(m_rng);
            }
        }

        return;
    }
};








// // [[Rcpp::export]]
// NumericMatrix stoch_test2(const double& mu,
//                           const double& sigma,
//                           const double& dt = 0.1,
//                           const double& max_t = 10.0) {
//
//     int32_t seed = static_cast<int32_t>(R::runif(0, 2147483647));
//
//     boost::mt11213b rng;
//     rng.seed(seed);
//     MatType x = { 0.5 };
//     Observer<MatType> obs;
//     boost::numeric::odeint::integrate_const(
//         StochLandscapeStepper() ,
//         std::make_pair(LandscapeConstF(mu) ,
//                        StochLandscapeStochProcess(rng, sigma)),
//                        x, 0.0, max_t, dt, std::ref(obs));
//
//     size_t n_steps = obs.data.size();
//     NumericMatrix output(n_steps, x.size()+1U);
//     colnames(output) = CharacterVector::create("t", "x");
//     for (size_t i = 0; i < n_steps; i++) {
//         output(i,0) = obs.time[i];
//         for (size_t j = 0; j < x.size(); j++) {
//             output(i,j+1U) = obs.data[i][j];
//         }
//     }
//     return output;
//
// }

//' @export
// [[Rcpp::export]]
NumericMatrix landscape_constantF_stoch_ode(const std::vector<double>& m,
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
                                            const double& sigma_y,
                                            const double& sigma_b,
                                            const double& dt = 0.1,
                                            const double& max_t = 100.0) {

    size_t np = m.size();
    /*
     I can't just use 'stop()' because it causes a segfault (or similar).
     The system below is my workaround.
     */
    bool err = lanscape_constF_arg_checks(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                                          L_0, u, X, Y0, B0, dt, max_t);
    min_val_check(err, sigma_y, "sigma_y", 0);
    min_val_check(err, sigma_b, "sigma_b", 0);
    if (err) return NumericMatrix(0,0);

    size_t n_states = 2U;
    MatType x(np, n_states);
    for (size_t i = 0; i < np; i++) {
        x(i,0) = Y0[i];
        x(i,1) = B0[i];
    }


    int32_t seed = static_cast<int32_t>(R::unif_rand() * 2147483647.0);
    boost::mt11213b rng;
    rng.seed(seed);

    Observer<MatType> obs;

    LandscapeConstF determ_sys(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, u, X);

    boost::numeric::odeint::integrate_const(
        StochLandscapeStepper(np, n_states) ,
        std::make_pair(determ_sys,
                       StochLandscapeStochProcess(rng, sigma_y, sigma_b)),
                       x, 0.0, max_t, dt, std::ref(obs));

    size_t n_steps = obs.data.size();
    NumericMatrix output(n_steps * np, n_states+3U);
    colnames(output) = CharacterVector::create("t", "p", "Y", "B", "P");
    std::vector<double> wts(np);
    size_t i = 0;
    for (size_t t = 0; t < n_steps; t++) {
        determ_sys.make_weights(wts, obs.data[t]);
        for (size_t k = 0; k < np; k++) {
            output(i,0) = obs.time[t];
            output(i,1) = k;
            for (size_t j = 0; j < n_states; j++) {
                output(i,j+2U) = obs.data[t](k,j);
            }
            output(i, n_states+2U) = wts[k];
            i++;
        }

    }
    return output;
}
