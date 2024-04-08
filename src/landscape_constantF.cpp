
/*
 Lanscape of plants where the total number of flowers (F) at each plant is
 constant and where pollinators are not affected by flower densities.
 We don't even include F and instead model Y and B as proportions
 of total flowers with non-colonized N = 1 - Y - B.
 */

#include <Rcpp.h>
#include <vector>

#include "ode.h"

using namespace Rcpp;




class LandConstFSystemFunction
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
    double q;
    double X;
    size_t n_plants;


    LandConstFSystemFunction(const std::vector<double>& m_,
                            const std::vector<double>& d_yp_,
                            const std::vector<double>& d_b0_,
                            const std::vector<double>& d_bp_,
                            const std::vector<double>& g_yp_,
                            const std::vector<double>& g_b0_,
                            const std::vector<double>& g_bp_,
                            const std::vector<double>& L_0_,
                            const double& q_,
                            const double& X_)
        : m(m_),
          d_yp(d_yp_),
          d_b0(d_b0_),
          d_bp(d_bp_),
          g_yp(g_yp_),
          g_b0(g_b0_),
          g_bp(g_bp_),
          L_0(L_0_),
          q(q_),
          X(X_),
          n_plants(m_.size()),
          weights(m_.size()) {};

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

        if (wts_vec.size() != n_plants) wts_vec.resize(n_plants);

        double wt_sum = 0;
        for (size_t i = 0; i < n_plants; i++) {
            const double& B(x(i,1));
            wts_vec[i] = std::exp(-q * B);
            wt_sum += wts_vec[i];
        }

        for (double& w : wts_vec) w /= (X + wt_sum);

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
        double N = 1 - Y - B;

        double& dYdt(dxdt(i,0));
        double& dBdt(dxdt(i,1));

        double P = weights[i];

        double Lambda = P / (L_0[i] + P);

        double gamma_y = g_yp[i] * Lambda;
        double gamma_b = g_b0[i] + g_bp[i] * Lambda;

        double delta_y = d_yp[i] * Lambda;
        double delta_b = d_b0[i] + d_bp[i] * Lambda;

        double disp_y = delta_y * Y + gamma_y;
        double disp_b = delta_b * B + gamma_b;

        dYdt = disp_y * N - m[i] * Y;
        dBdt = disp_b * N - m[i] * B;

        return;
    }

};


struct LandConstFObserver
{
    std::vector<MatType> data;
    std::vector<double> time;
    LandConstFObserver() : data(), time() {};

    void operator()(const MatType& x, const double& t) {
        data.push_back(x);
        time.push_back(t);
        return;
    }
};





//' @export
// [[Rcpp::export]]
NumericMatrix landscape_constantF_ode(const std::vector<double>& m,
                                      const std::vector<double>& d_yp,
                                      const std::vector<double>& d_b0,
                                      const std::vector<double>& d_bp,
                                      const std::vector<double>& g_yp,
                                      const std::vector<double>& g_b0,
                                      const std::vector<double>& g_bp,
                                      const std::vector<double>& L_0,
                                      const double& q,
                                      const double& X,
                                      const std::vector<double>& Y0,
                                      const std::vector<double>& B0,
                                      const double& dt = 0.1,
                                      const double& max_t = 90.0) {

    size_t np = m.size();
    if(d_yp.size() != np) stop("d_yp is the wrong length!");
    if(d_b0.size() != np) stop("d_b0 is the wrong length!");
    if(d_bp.size() != np) stop("d_bp is the wrong length!");
    if(g_yp.size() != np) stop("g_yp is the wrong length!");
    if(g_b0.size() != np) stop("g_b0 is the wrong length!");
    if(g_bp.size() != np) stop("g_bp is the wrong length!");
    if(L_0.size() != np) stop("L_0 is the wrong length!");
    if(Y0.size() != np) stop("Y0 is the wrong length!");
    if(B0.size() != np) stop("B0 is the wrong length!");

    size_t n_states = 2U;
    MatType x(np, n_states);
    for (size_t i = 0; i < np; i++) {
        x(i,0) = Y0[i];
        x(i,1) = B0[i];
    }


    LandConstFObserver obs;
    LandConstFSystemFunction system(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                                    L_0, q, X);

    boost::numeric::odeint::integrate_const(
        MatStepperType(), std::ref(system),
        x, 0.0, max_t, dt, std::ref(obs));

    size_t n_steps = obs.data.size();
    NumericMatrix output(n_steps * np, n_states+3U);
    colnames(output) = CharacterVector::create("t", "p", "Y", "B", "P");
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
            output(i, n_states+2U) = wts[k];
            i++;
        }

    }
    return output;
}
