#define _USE_MATH_DEFINES


#include <Rcpp.h>
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
    double S_0;
    double q;
    std::vector<double> X;
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
                            const double& S_0_,
                            const double& q_,
                            const std::vector<double>& X_,
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
          S_0(S_0_),
          q(q_),
          X(X_),
          R_hat(R_hat_),
          mu(mu_),
          sigma(sigma_),
          exp_wz(z_.size1(), z_.size2()),
          n_plants(z_.size1()),
          Y0(Y0_),
          B0(B0_),
          add_F(add_F_),
          YB_added(z_.size1(), false),
          weights(z_.size1()) {
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

        landscape_weights__(wts_vec, x, n_plants, S_0, q, X, exp_wz);

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


struct LandscapeSeasonObserver
{
    std::vector<MatType> data;
    std::vector<double> time;
    LandscapeSeasonObserver() : data(), time() {};

    void operator()(const MatType& x, const double& t) {
        data.push_back(x);
        time.push_back(t);
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
                                   const double& S_0,
                                   const double& q,
                                   const std::vector<double>& X,
                                   const std::vector<double>& R_hat,
                                   const std::vector<double>& mu,
                                   const std::vector<double>& sigma,
                                   const double& w,
                                   const NumericMatrix& z,
                                   const std::vector<double>& Y0,
                                   const std::vector<double>& B0,
                                   const double& add_F = 1.0,
                                   const double& dt = 0.1,
                                   const double& max_t = 90.0) {

    size_t np = z.nrow();
    /*
     I can't just use 'stop()' because it causes a segfault (or similar).
     The system below is my workaround.
     */
    bool err = false;
    if (z.ncol() != np) {
        Rcout << "z needs to be square!" << std::endl;
        err = true;
    }
    len_check<double>(err, m, "m", np);
    len_check<double>(err, d_yp, "d_yp", np);
    len_check<double>(err, d_b0, "d_b0", np);
    len_check<double>(err, d_bp, "d_bp", np);
    len_check<double>(err, g_yp, "g_yp", np);
    len_check<double>(err, g_b0, "g_b0", np);
    len_check<double>(err, g_bp, "g_bp", np);
    len_check<double>(err, L_0, "L_0", np);
    len_check<double>(err, P_max, "P_max", np);
    len_check<double>(err, X, "X", np);
    len_check<double>(err, R_hat, "R_hat", np);
    len_check<double>(err, mu, "mu", np);
    len_check<double>(err, sigma, "sigma", np);
    len_check<double>(err, Y0, "Y0", np);
    len_check<double>(err, B0, "B0", np);
    if ((*std::min_element(Y0.begin(), Y0.end())) < 0) {
        Rcout << "Y0 must only contain values >= 0." << std::endl;
        err = true;
    }
    if ((*std::min_element(B0.begin(), B0.end())) < 0) {
        Rcout << "B0 must only contain values >= 0." << std::endl;
        err = true;
    }
    if (B0.size() == Y0.size()) {
        for (size_t i = 0; i < Y0.size(); i++) {
            if ((Y0[i] + B0[i]) > 1) {
                Rcout << "Y0+B0 must always be <= 1." << std::endl;
                err = true;
                break;
            }
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

    MatType z_boost(np, np);
    for (size_t i = 0; i < np; i++) {
        for (size_t j = 0; j < np; j++) {
            z_boost(i,j) = z(i,j);
        }
    }

    LandscapeSeasonObserver obs;
    LandscapeSeasonSystemFunction system(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                                   L_0, P_max, S_0, q, X, R_hat, mu, sigma, w, z_boost,
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
