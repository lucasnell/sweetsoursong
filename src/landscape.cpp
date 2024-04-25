#include "ode.h"

#include <RcppArmadillo.h>
#include <vector>

using namespace Rcpp;



class LandscapeSystemFunction
{
public:
    std::vector<double> m;
    std::vector<double> R;
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
    MatType exp_wz;
    size_t n_plants;


    LandscapeSystemFunction(const std::vector<double>& m_,
                            const std::vector<double>& R_,
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
                            const double& w_,
                            const MatType& z_)
        : m(m_),
          R(R_),
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
          exp_wz(z_.n_rows, z_.n_cols),
          n_plants(z_.n_rows),
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

        landscape_weights__(wts_vec, x, n_plants, S_0, q, X, exp_wz);

        return;
    }



private:

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

        double PF = P / F;
        double Lambda = PF / (L_0[i] + PF);

        double gamma_y = g_yp[i] * Lambda;
        double gamma_b = g_b0[i] + g_bp[i] * Lambda;

        double delta_y = d_yp[i] * Lambda;
        double delta_b = d_b0[i] + d_bp[i] * Lambda;

        double disp_y = delta_y * Y / F + gamma_y;
        double disp_b = delta_b * B / F + gamma_b;

        dYdt = disp_y * N - m[i] * Y;
        dBdt = disp_b * N - m[i] * B;
        dNdt = R[i] - N * (m[i] + disp_y + disp_b);

        return;
    }

};


struct LandscapeObserver
{
    std::vector<MatType> data;
    std::vector<double> time;
    LandscapeObserver() : data(), time() {};

    void operator()(const MatType& x, const double& t) {
        data.push_back(x);
        time.push_back(t);
        return;
    }
};



//' @export
// [[Rcpp::export]]
NumericMatrix landscape_ode(const std::vector<double>& m,
                            const std::vector<double>& R,
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
                            const double& w,
                            const arma::mat& z,
                            const std::vector<double>& Y0,
                            const std::vector<double>& B0,
                            const std::vector<double>& N0,
                            const double& dt = 0.1,
                            const double& max_t = 90.0) {

    size_t np = z.n_rows;
    /*
     I can't just use 'stop()' because it causes a segfault (or similar).
     The system below is my workaround.
     */
    bool err = false;
    if (z.n_cols != np) {
        Rcout << "z needs to be square!" << std::endl;
        err = true;
    }
    len_check<double>(err, m, "m", np);
    len_check<double>(err, R, "R", np);
    len_check<double>(err, d_yp, "d_yp", np);
    len_check<double>(err, d_b0, "d_b0", np);
    len_check<double>(err, d_bp, "d_bp", np);
    len_check<double>(err, g_yp, "g_yp", np);
    len_check<double>(err, g_b0, "g_b0", np);
    len_check<double>(err, g_bp, "g_bp", np);
    len_check<double>(err, L_0, "L_0", np);
    len_check<double>(err, P_max, "P_max", np);
    len_check<double>(err, X, "X", np);
    len_check<double>(err, Y0, "Y0", np);
    len_check<double>(err, B0, "B0", np);
    len_check<double>(err, N0, "N0", np);

    if (err) return NumericMatrix(0,0);

    size_t n_states = 3U;
    MatType x(np, n_states);
    for (size_t i = 0; i < np; i++) {
        x(i,0) = Y0[i];
        x(i,1) = B0[i];
        x(i,2) = N0[i];
    }

    LandscapeObserver obs;
    LandscapeSystemFunction system(m, R, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                                   L_0, P_max, S_0, q, X, w, z);

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
