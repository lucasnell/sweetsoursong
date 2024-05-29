
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

using namespace Rcpp;



// Geometric Brownian motion
class GeometricBrownian
{
public:

    typedef VecType state_type;
    typedef VecType deriv_type;
    typedef double value_type;
    typedef double time_type;
    typedef unsigned short order_type;

    typedef boost::numeric::odeint::stepper_tag stepper_category;

    static order_type order( void ) { return 1; }

    template< class System >
    void do_step( System system , state_type &x , time_type t , time_type dt ) const {
        deriv_type det(x.size()), stoch(x.size());
        system.first(x, det);
        system.second(x, stoch);
        for (size_t i = 0 ; i < x.size() ; i++) {
            x[i] += dt * det[i] + std::sqrt(dt) * stoch[i];
            // if (x[i] > 1) x[i] = 1;
            // if (x[i] < 0) x[i] = 0;
        }
    }
};
//]


typedef VecType state_type;

struct GeometricBrownian_det
{

    double mu;

    GeometricBrownian_det(const double& mu_) : mu(mu_) {}

    void operator()(const state_type &x, state_type &dxdt) const {
        double tau = 1;
        dxdt[0] = -(x[0] - mu) / tau;
        return;
    }
};

struct GeometricBrownian_stoch
{
    boost::mt11213b &m_rng;
    boost::normal_distribution<> m_dist;

    GeometricBrownian_stoch(boost::mt11213b &rng, double sigma)
        : m_rng(rng) , m_dist(0.0, sigma) {}

    void operator()(const state_type &x, state_type &dxdt) {
        // double stdev = std::sin(M_PI * x[0]);
        double stdev = 4 * x[0] * (1 - x[0]);
        // double stdev = 1 - 2 * std::abs(x[0] - 0.5);
        // double stdev = 1;
        dxdt[0] = stdev * m_dist(m_rng);
        return;
    }
};


 // [[Rcpp::export]]
 NumericMatrix stoch_test2(const double& mu,
                           const double& sigma,
                           const double& dt = 0.1,
                           const double& max_t = 10.0) {

     int32_t seed = static_cast<int32_t>(R::runif(0, 2147483647));

     boost::mt11213b rng;
     rng.seed(seed);
     state_type x = { 0.5 };
     Observer<state_type> obs;
     boost::numeric::odeint::integrate_const(
         GeometricBrownian() ,
         std::make_pair(GeometricBrownian_det(mu) ,
                        GeometricBrownian_stoch(rng, sigma)),
         x, 0.0, max_t, dt, std::ref(obs));

     size_t n_steps = obs.data.size();
     NumericMatrix output(n_steps, x.size()+1U);
     colnames(output) = CharacterVector::create("t", "x");
     for (size_t i = 0; i < n_steps; i++) {
         output(i,0) = obs.time[i];
         for (size_t j = 0; j < x.size(); j++) {
             output(i,j+1U) = obs.data[i][j];
         }
     }
     return output;

 }



// class DetermLandConstFSysFun
// {
// public:
//     std::vector<double> m;
//     std::vector<double> d_yp;
//     std::vector<double> d_b0;
//     std::vector<double> d_bp;
//     std::vector<double> g_yp;
//     std::vector<double> g_b0;
//     std::vector<double> g_bp;
//     std::vector<double> L_0;
//     double u;
//     double X;
//     size_t n_plants;
//
//
//     DetermLandConstFSysFun(const std::vector<double>& m_,
//                            const std::vector<double>& d_yp_,
//                            const std::vector<double>& d_b0_,
//                            const std::vector<double>& d_bp_,
//                            const std::vector<double>& g_yp_,
//                            const std::vector<double>& g_b0_,
//                            const std::vector<double>& g_bp_,
//                            const std::vector<double>& L_0_,
//                            const double& u_,
//                            const double& X_)
//         : m(m_),
//           d_yp(d_yp_),
//           d_b0(d_b0_),
//           d_bp(d_bp_),
//           g_yp(g_yp_),
//           g_b0(g_b0_),
//           g_bp(g_bp_),
//           L_0(L_0_),
//           u(u_),
//           X(X_),
//           n_plants(m_.size()),
//           weights(m_.size()) {};
//
//     void operator()(const MatType& x,
//                   MatType& dxdt,
//                   const double t) {
//
//         make_weights(this->weights, x);
//
//         for (size_t i = 0; i < n_plants; i++) {
//             one_plant(i, x, dxdt, t);
//         }
//
//         return;
//     }
//
//     void make_weights(std::vector<double>& wts_vec,
//                       const MatType& x) {
//
//         if (wts_vec.size() != n_plants) wts_vec.resize(n_plants);
//
//         double wt_sum = 0;
//         double YN; // proportion palatable nectar
//         for (size_t i = 0; i < n_plants; i++) {
//             YN = 1 - x(i,1);
//             wts_vec[i] = std::pow(YN, u);
//             wt_sum += wts_vec[i];
//         }
//
//         for (double& w : wts_vec) w /= (X + wt_sum);
//
//         return;
//     }
//
//
//
// private:
//
//     std::vector<double> weights;
//
//     void one_plant(const size_t& i,
//                    const MatType& x,
//                    MatType& dxdt,
//                    const double t) {
//
//         const double& Y(x(i,0));
//         const double& B(x(i,1));
//         double N = 1 - Y - B;
//
//         double& dYdt(dxdt(i,0));
//         double& dBdt(dxdt(i,1));
//
//         double P = weights[i];
//
//         double Lambda = P / (L_0[i] + P);
//
//         double gamma_y = g_yp[i] * Lambda;
//         double gamma_b = g_b0[i] + g_bp[i] * Lambda;
//
//         double delta_y = d_yp[i] * Lambda;
//         double delta_b = d_b0[i] + d_bp[i] * Lambda;
//
//         double disp_y = delta_y * Y + gamma_y;
//         double disp_b = delta_b * B + gamma_b;
//
//         dYdt = disp_y * N - m[i] * Y;
//         dBdt = disp_b * N - m[i] * B;
//
//         return;
//     }
//
// };
//
//
//
//
//
//
// //' @export
// // [[Rcpp::export]]
// NumericMatrix landscape_constantF_stoch_ode(const std::vector<double>& m,
//                                             const std::vector<double>& d_yp,
//                                             const std::vector<double>& d_b0,
//                                             const std::vector<double>& d_bp,
//                                             const std::vector<double>& g_yp,
//                                             const std::vector<double>& g_b0,
//                                             const std::vector<double>& g_bp,
//                                             const std::vector<double>& L_0,
//                                             const double& u,
//                                             const double& X,
//                                             const std::vector<double>& Y0,
//                                             const std::vector<double>& B0,
//                                             const double& sigma_y,
//                                             const double& sigma_b,
//                                             const double& dt = 0.1,
//                                             const double& max_t = 90.0) {
//
//     size_t np = m.size();
//     /*
//      I can't just use 'stop()' because it causes a segfault (or similar).
//      The system below is my workaround.
//      */
//     bool err = false;
//     len_check(err, d_yp, "d_yp", np);
//     len_check(err, d_b0, "d_b0", np);
//     len_check(err, d_bp, "d_bp", np);
//     len_check(err, g_yp, "g_yp", np);
//     len_check(err, g_b0, "g_b0", np);
//     len_check(err, g_bp, "g_bp", np);
//     len_check(err, L_0, "L_0", np);
//     len_check(err, Y0, "Y0", np);
//     len_check(err, B0, "B0", np);
//
//     if (err) return NumericMatrix(0,0);
//
//     size_t n_states = 2U;
//     MatType x(np, n_states);
//     for (size_t i = 0; i < np; i++) {
//         x(i,0) = Y0[i];
//         x(i,1) = B0[i];
//     }
//
//
//     Observer<MatType> obs;
//     DetermLandConstFSysFun system(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
//                                  L_0, u, X);
//
//     boost::numeric::odeint::integrate_const(
//         MatStepperType(), std::ref(system),
//         x, 0.0, max_t, dt, std::ref(obs));
//
//     size_t n_steps = obs.data.size();
//     NumericMatrix output(n_steps * np, n_states+3U);
//     colnames(output) = CharacterVector::create("t", "p", "Y", "B", "P");
//     std::vector<double> wts(np);
//     size_t i = 0;
//     for (size_t t = 0; t < n_steps; t++) {
//         system.make_weights(wts, obs.data[t]);
//         for (size_t k = 0; k < np; k++) {
//             output(i,0) = obs.time[t];
//             output(i,1) = k;
//             for (size_t j = 0; j < n_states; j++) {
//                 output(i,j+2U) = obs.data[t](k,j);
//             }
//             output(i, n_states+2U) = wts[k];
//             i++;
//         }
//
//     }
//     return output;
// }
