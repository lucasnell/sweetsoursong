# ifndef __SWEETSOURSONG_LANDSCAPE_CONSTANTF_H
# define __SWEETSOURSONG_LANDSCAPE_CONSTANTF_H



/*
 Lanscape of plants where the total number of flowers (F) at each plant is
 constant and where pollinators are not affected by flower densities.
 We don't even include F and instead model Y and B as proportions
 of total flowers with non-colonized N = 1 - Y - B.
 */

#include <RcppArmadillo.h>
#include <vector>

#include "ode.h"

using namespace Rcpp;



inline bool lanscape_constF_arg_checks(const std::vector<double>& m,
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
                                       const double& dt,
                                       const double& max_t) {

    size_t np = m.size();
    /*
     I can't just use 'stop()' because it causes a segfault (or similar).
     The system below is my workaround.
     */
    bool err = false;

    len_check(err, d_yp, "d_yp", np);
    len_check(err, d_b0, "d_b0", np);
    len_check(err, d_bp, "d_bp", np);
    len_check(err, g_yp, "g_yp", np);
    len_check(err, g_b0, "g_b0", np);
    len_check(err, g_bp, "g_bp", np);
    len_check(err, L_0, "L_0", np);
    len_check(err, Y0, "Y0", np);
    len_check(err, B0, "B0", np);

    min_val_check(err, d_yp, "d_yp", 0);
    min_val_check(err, d_b0, "d_b0", 0);
    min_val_check(err, d_bp, "d_bp", 0);
    min_val_check(err, g_yp, "g_yp", 0);
    min_val_check(err, g_b0, "g_b0", 0);
    min_val_check(err, g_bp, "g_bp", 0);
    min_val_check(err, L_0, "L_0", 0);
    min_val_check(err, u, "u", 0);
    min_val_check(err, X, "X", 0);

    min_val_check(err, Y0, "Y0", 0);
    min_val_check(err, B0, "B0", 0);
    max_val_check(err, Y0, "Y0", 1);
    max_val_check(err, B0, "B0", 1);

    min_val_check(err, dt, "dt", 0, false);
    min_val_check(err, max_t, "max_t", std::max(dt, 0.0), false);

    return err;

}






class LandscapeConstF
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
    double u;
    double X;
    size_t n_plants;
    bool open_sys;


    LandscapeConstF(const std::vector<double>& m_,
                    const std::vector<double>& d_yp_,
                    const std::vector<double>& d_b0_,
                    const std::vector<double>& d_bp_,
                    const std::vector<double>& g_yp_,
                    const std::vector<double>& g_b0_,
                    const std::vector<double>& g_bp_,
                    const std::vector<double>& L_0_,
                    const double& u_,
                    const double& X_,
                    const bool& open_sys_)
        : m(m_),
          d_yp(d_yp_),
          d_b0(d_b0_),
          d_bp(d_bp_),
          g_yp(g_yp_),
          g_b0(g_b0_),
          g_bp(g_bp_),
          L_0(L_0_),
          u(u_),
          X(X_),
          n_plants(m_.size()),
          open_sys(open_sys_),
          weights(m_.size()) {};


    LandscapeConstF(const LandscapeConstF& other)
        : m(other.m),
          d_yp(other.d_yp),
          d_b0(other.d_b0),
          d_bp(other.d_bp),
          g_yp(other.g_yp),
          g_b0(other.g_b0),
          g_bp(other.g_bp),
          L_0(other.L_0),
          u(other.u),
          X(other.X),
          n_plants(other.n_plants),
          open_sys(other.open_sys),
          weights(other.weights) {};

    LandscapeConstF& operator=(const LandscapeConstF& other) {
        m = other.m;
        d_yp = other.d_yp;
        d_b0 = other.d_b0;
        d_bp = other.d_bp;
        g_yp = other.g_yp;
        g_b0 = other.g_b0;
        g_bp = other.g_bp;
        L_0 = other.L_0;
        u = other.u;
        X = other.X;
        n_plants = other.n_plants;
        open_sys = other.open_sys;
        weights = other.weights;
        return *this;
    }


    void operator()(const MatType& x,
                    MatType& dxdt) {

        make_weights(this->weights, x);

        double Ybar = arma::mean(x.col(0));
        double Bbar = arma::mean(x.col(1));

        for (size_t i = 0; i < n_plants; i++) {
            one_plant(i, Ybar, Bbar, x, dxdt);
        }

        return;
    }

    void operator()(const MatType& x,
                  MatType& dxdt,
                  const double& t) {
        this->operator()(x, dxdt);
        return;
    }


    void make_weights(std::vector<double>& wts_vec,
                      const MatType& x) {

        if (wts_vec.size() != n_plants) wts_vec.resize(n_plants);

        double wt_sum = 0;
        double YN; // proportion palatable nectar
        for (size_t i = 0; i < n_plants; i++) {
            YN = 1 - x(i,1);
            wts_vec[i] = std::pow(YN, u);
            wt_sum += wts_vec[i];
        }

        double mult = static_cast<double>(n_plants) / (2 * (X + wt_sum));
        for (double& w : wts_vec) w *= mult;

        return;
    }



private:

    std::vector<double> weights;

   void one_plant(const size_t& i,
                  const double& Ybar,
                  const double& Bbar,
                  const MatType& x,
                  MatType& dxdt) {

        const double& Y(x(i,0));
        const double& B(x(i,1));
        double N = 1 - Y - B;

        double& dYdt(dxdt(i,0));
        double& dBdt(dxdt(i,1));

        double P = weights[i];

        double Lambda = P / (L_0[i] + P);

        double delta, gamma;

        if (open_sys) {
            delta = (g_yp[i] + d_yp[i] * Y) * Lambda;
            gamma = g_b0[i] + (d_b0[i] * B) + (g_bp[i] + d_bp[i] * B) * Lambda;
        } else {
            delta = (g_yp[i] * Ybar + d_yp[i] * Y) * Lambda;
            gamma = g_b0[i] * Bbar + (d_b0[i] * B) + (g_bp[i] * Bbar + d_bp[i] * B) * Lambda;
        }

        dYdt = delta * N - m[i] * Y;
        dBdt = gamma * N - m[i] * B;

        return;
    }

};





#endif
