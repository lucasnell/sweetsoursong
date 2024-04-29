#define _USE_MATH_DEFINES


#include <RcppArmadillo.h>
#include <vector>
#include <cmath>

#include "ode.h"

using namespace Rcpp;



class LandscapeSeasonSystemFunction
{
public:
    std::vector<double> m;
    std::vector<double> d_yp;
    std::vector<double> d_b0;
    std::vector<double> d_bp;
    std::vector<double> g_yp;
    std::vector<double> g_b0;
    std::vector<double> g_bp;
    std::vector<double> L_0;
    std::vector<double> P_max;
    double u;
    double q;
    std::vector<double> W;
    std::vector<double> R_hat;
    std::vector<double> mu;
    std::vector<double> sigma;
    MatType exp_wz;
    size_t n_plants;


    LandscapeSeasonSystemFunction(const std::vector<double>& m_,
                            const std::vector<double>& d_yp_,
                            const std::vector<double>& d_b0_,
                            const std::vector<double>& d_bp_,
                            const std::vector<double>& g_yp_,
                            const std::vector<double>& g_b0_,
                            const std::vector<double>& g_bp_,
                            const std::vector<double>& L_0_,
                            const std::vector<double>& P_max_,
                            const double& u_,
                            const double& q_,
                            const std::vector<double>& W_,
                            const std::vector<double>& R_hat_,
                            const std::vector<double>& mu_,
                            const std::vector<double>& sigma_,
                            const double& w_,
                            const MatType& z_,
                            const std::vector<double>& Y0_,
                            const std::vector<double>& B0_,
                            const double& add_F_)
        : m(m_),
          d_yp(d_yp_),
          d_b0(d_b0_),
          d_bp(d_bp_),
          g_yp(g_yp_),
          g_b0(g_b0_),
          g_bp(g_bp_),
          L_0(L_0_),
          P_max(P_max_),
          u(u_),
          q(q_),
          W(W_),
          R_hat(R_hat_),
          mu(mu_),
          sigma(sigma_),
          exp_wz(z_.n_rows, z_.n_cols),
          n_plants(z_.n_rows),
          Y0(Y0_),
          B0(B0_),
          add_F(add_F_),
          YB_added(z_.n_rows, false),
          weights(z_.n_rows) {
        for (size_t i = 0; i < n_plants; i++) {
            for (size_t j = 0; j < n_plants; j++) {
                if (i == j) {
                    exp_wz(i,j) = 1;
                } else {
                    exp_wz(i,j) = std::exp(-w_ * z_(i, j));
                }
            }
        }
    };

    void operator()(const MatType& x,
                    MatType& dxdt,
                    const double t) {

        make_weights(this->weights, x);

        for (size_t i = 0; i < n_plants; i++) {
            one_plant(i, x, dxdt, t);
        }

        return;
    }

    void make_weights(std::vector<double>& wts_vec,
                      const MatType& x) {

        landscape_weights__(wts_vec, x, n_plants, u, q, W);

        return;
    }



private:

    std::vector<double> Y0;
    std::vector<double> B0;
    double add_F;
    std::vector<bool> YB_added;
    std::vector<double> weights;

   void one_plant(const size_t& i,
                   const MatType& x,
                   MatType& dxdt,
                   const double t) {

        const double& Y(x(i,0));
        const double& B(x(i,1));
        const double& N(x(i,2));

        double& dYdt(dxdt(i,0));
        double& dBdt(dxdt(i,1));
        double& dNdt(dxdt(i,2));

        double F = Y + B + N;

        double P = P_max[i] * weights[i];

        double PF = 0;
        double YF = 0;
        double BF = 0;
        if (F > 0) {
            PF = P / F;
            YF = Y / F;
            BF = B / F;
        }
        double Lambda = PF / (L_0[i] + PF);
        if (L_0[i] == 0 && F == 0) Lambda = 0;

        double gamma_y = g_yp[i] * Lambda;
        double gamma_b = g_b0[i] + g_bp[i] * Lambda;

        double delta_y = d_yp[i] * Lambda;
        double delta_b = d_b0[i] + d_bp[i] * Lambda;

        double disp_y = delta_y * YF + gamma_y;
        double disp_b = delta_b * BF + gamma_b;

        double R = (R_hat[i] / (sigma[i] * std::sqrt(2 * M_PI))) *
            std::exp(-0.5 * std::pow((t - mu[i]) / sigma[i], 2U));

        dYdt = disp_y * N - m[i] * Y;
        dBdt = disp_b * N - m[i] * B;
        dNdt = R - N * (m[i] + disp_y + disp_b);

        if (F >= add_F && ! YB_added[i]) {
            // note: using 'N + dNdt' below to avoid N going < 0
            double y0_i = (N + dNdt) * Y0[i];
            double b0_i = (N + dNdt) * B0[i];
            dYdt += y0_i;
            dBdt += b0_i;
            dNdt -= (y0_i + b0_i);
            YB_added[i] = true;
        }

        return;
    }

};





//' @export
// [[Rcpp::export]]
NumericMatrix landscape_season_ode(const std::vector<double>& m,
                                   const std::vector<double>& d_yp,
                                   const std::vector<double>& d_b0,
                                   const std::vector<double>& d_bp,
                                   const std::vector<double>& g_yp,
                                   const std::vector<double>& g_b0,
                                   const std::vector<double>& g_bp,
                                   const std::vector<double>& L_0,
                                   const std::vector<double>& P_max,
                                   const double& u,
                                   const double& q,
                                   const std::vector<double>& W,
                                   const std::vector<double>& R_hat,
                                   const std::vector<double>& mu,
                                   const std::vector<double>& sigma,
                                   const double& w,
                                   const arma::mat& z,
                                   const std::vector<double>& Y0,
                                   const std::vector<double>& B0,
                                   const double& add_F = 1.0,
                                   const double& dt = 0.1,
                                   const double& max_t = 90.0) {

    size_t np = z.n_rows;
    /*
     I can't just use 'stop()' because it causes a segfault (or similar).
     The system below is my workaround.
     */
    bool err = false;
    mat_dim_check(err, z, "z", np);
    len_check(err, m, "m", np);
    len_check(err, d_yp, "d_yp", np);
    len_check(err, d_b0, "d_b0", np);
    len_check(err, d_bp, "d_bp", np);
    len_check(err, g_yp, "g_yp", np);
    len_check(err, g_b0, "g_b0", np);
    len_check(err, g_bp, "g_bp", np);
    len_check(err, L_0, "L_0", np);
    len_check(err, P_max, "P_max", np);
    len_check(err, W, "W", np);
    len_check(err, R_hat, "R_hat", np);
    len_check(err, mu, "mu", np);
    len_check(err, sigma, "sigma", np);
    len_check(err, Y0, "Y0", np);
    len_check(err, B0, "B0", np);

    // Must be (or contain values) at least certain value (usually zero):
    min_val_check(err, z, "z", 0);
    min_val_check(err, m, "m", 0, false);
    min_val_check(err, d_yp, "d_yp", 0);
    min_val_check(err, d_b0, "d_b0", 0);
    min_val_check(err, d_bp, "d_bp", 0);
    min_val_check(err, g_yp, "g_yp", 0);
    min_val_check(err, g_b0, "g_b0", 0);
    min_val_check(err, g_bp, "g_bp", 0);
    min_val_check(err, L_0, "L_0", 0, false);
    min_val_check(err, P_max, "P_max", 0, false);
    min_val_check(err, W, "W", 0);
    min_val_check(err, R_hat, "R_hat", 0, false);
    min_val_check(err, mu, "mu", 0, false);
    min_val_check(err, sigma, "sigma", 0, false);
    min_val_check(err, Y0, "Y0", 0);
    min_val_check(err, B0, "B0", 0);
    min_val_check(err, u, "u", 0);
    min_val_check(err, q, "q", 0);
    min_val_check(err, w, "w", 0);
    min_val_check(err, add_F, "add_F", 0);
    min_val_check(err, dt, "dt", 0, false);
    min_val_check(err, max_t, "max_t", 1.0);
    for (size_t i = 0; i < std::min(B0.size(), Y0.size()); i++) {
        if ((Y0[i] + B0[i]) > 1) {
            Rcout << "Y0+B0 must always be <= 1." << std::endl;
            err = true;
            break;
        }
    }

    if (err) return NumericMatrix(0,0);

    size_t n_states = 3U;
    MatType x(np, n_states);
    for (size_t i = 0; i < np; i++) {
        x(i,0) = 0.0;  // Y0[i];
        x(i,1) = 0.0;  // B0[i];
        x(i,2) = 0.0;  // N0[i];
    }

    Observer<MatType> obs;
    LandscapeSeasonSystemFunction system(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                                   L_0, P_max, u, q, W, R_hat, mu, sigma, w, z,
                                   Y0, B0, add_F);

    boost::numeric::odeint::integrate_const(
        MatStepperType(), std::ref(system),
        x, 0.0, max_t, dt, std::ref(obs));

    size_t n_steps = obs.data.size();
    NumericMatrix output(n_steps * np, n_states+3U);
    colnames(output) = CharacterVector::create("t", "p", "Y", "B", "N", "P");
    std::vector<double> wts(np);
    size_t i = 0;
    for (size_t t = 0; t < n_steps; t++) {
        system.make_weights(wts, obs.data[t]);
        for (size_t k = 0; k < np; k++) {
            output(i,0) = obs.time[t];
            output(i,1) = k;
            for (size_t j = 0; j < n_states; j++) {
                output(i,j+2U) = obs.data[t](k,j);
            }
            output(i, n_states+2U) = P_max[k] * wts[k];
            i++;
        }

    }
    return output;
}
