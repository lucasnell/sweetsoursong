#include <Rcpp.h>
#include <vector>

#include "ode.h"

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
          q(q_),
          X(X_),
          exp_wz(z_.size1(), z_.size2()),
          n_plants(z_.size1()),
          weights(z_.size1()) {
        for (size_t i = 0; i < n_plants; i++) {
            for (size_t j = 0; j < n_plants; j++) {
                exp_wz(i,j) = std::exp(-w_ * z_(i, j));
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

        std::vector<double> F;
        F.reserve(n_plants);
        for (size_t i = 0; i < n_plants; i++) {
            F.push_back(x(i,0) + x(i,1) + x(i,2));
        }

        double theta, wt_sum = 0;
        for (size_t i = 0; i < n_plants; i++) {
            theta = F[i];
            for (size_t j = 0; j < n_plants; j++) {
                if (j == i) continue;
                theta += (F[j] * exp_wz(i,j));
            }
            const double& B_i(x(i,1));
            wts_vec.push_back(std::exp(theta - q * B_i / F[i]));
            wt_sum += wts_vec.back();
        }

        for (size_t i = 0; i < n_plants; i++) {
            wts_vec[i] /= (wt_sum + X[i]);
        }

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
                            const double& q,
                            const std::vector<double>& X,
                            const double& w,
                            const NumericMatrix& z,
                            const std::vector<double>& Y0,
                            const std::vector<double>& B0,
                            const std::vector<double>& N0,
                            const double& dt = 0.1,
                            const double& max_t = 90.0) {

    size_t np = z.nrow();
    if (z.ncol() != np) stop("z needs to be square!");
    if(m.size() != np) stop("m is the wrong length!");
    if(R.size() != np) stop("R is the wrong length!");
    if(d_yp.size() != np) stop("d_yp is the wrong length!");
    if(d_b0.size() != np) stop("d_b0 is the wrong length!");
    if(d_bp.size() != np) stop("d_bp is the wrong length!");
    if(g_yp.size() != np) stop("g_yp is the wrong length!");
    if(g_b0.size() != np) stop("g_b0 is the wrong length!");
    if(g_bp.size() != np) stop("g_bp is the wrong length!");
    if(L_0.size() != np) stop("L_0 is the wrong length!");
    if(P_max.size() != np) stop("P_max is the wrong length!");
    if(X.size() != np) stop("X is the wrong length!");
    if(Y0.size() != np) stop("Y0 is the wrong length!");
    if(B0.size() != np) stop("B0 is the wrong length!");
    if(N0.size() != np) stop("N0 is the wrong length!");

    size_t n_states = 3U;
    MatType x(np, n_states);
    for (size_t i = 0; i < np; i++) {
        x(i,0) = Y0[i];
        x(i,1) = B0[i];
        x(i,2) = N0[i];
    }

    MatType z_boost(np, np);
    for (size_t i = 0; i < np; i++) {
        for (size_t j = 0; j < np; j++) {
            z_boost(i,j) = z(i,j);
        }
    }

    LandscapeObserver obs;
    LandscapeSystemFunction system(m, R, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                                   L_0, P_max, q, X, w, z_boost);

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
