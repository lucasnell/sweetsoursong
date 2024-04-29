# ifndef __SWEETSOURSONG_LANDSCAPE_H
# define __SWEETSOURSONG_LANDSCAPE_H


#include <RcppArmadillo.h>
#include <vector>

#include "ode.h"


using namespace Rcpp;


inline bool lanscape_arg_checks(const std::vector<double>& m,
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
                                const double& w,
                                const arma::mat& z,
                                const std::vector<double>& Y0,
                                const std::vector<double>& B0,
                                const double& dt,
                                const double& max_t) {
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
    len_check(err, Y0, "Y0", np);
    len_check(err, B0, "B0", np);

    // Must be (or contain values) at least certain value (usually zero):

    min_val_check(err, m, "m", 0, false);
    min_val_check(err, d_yp, "d_yp", 0);
    min_val_check(err, d_b0, "d_b0", 0);
    min_val_check(err, d_bp, "d_bp", 0);
    min_val_check(err, g_yp, "g_yp", 0);
    min_val_check(err, g_b0, "g_b0", 0);
    min_val_check(err, g_bp, "g_bp", 0);
    min_val_check(err, L_0, "L_0", 0, false);
    min_val_check(err, P_max, "P_max", 0, false);
    min_val_check(err, u, "u", 0);
    min_val_check(err, q, "q", 0);
    min_val_check(err, W, "W", 0);
    min_val_check(err, w, "w", 0);
    min_val_check(err, z, "z", 0);
    min_val_check(err, Y0, "Y0", 0);
    min_val_check(err, B0, "B0", 0);
    min_val_check(err, dt, "dt", 0, false);
    min_val_check(err, max_t, "max_t", 1.0);


    return err;
}






class LandscapeSystemFunction
{
public:
    arma::vec m;
    arma::vec d_yp;
    arma::vec d_b0;
    arma::vec d_bp;
    arma::vec g_yp;
    arma::vec g_b0;
    arma::vec g_bp;
    arma::vec L_0;
    arma::vec P_max;
    double u;
    double q;
    arma::vec W;
    arma::mat Phi;
    size_t n_plants;



    LandscapeSystemFunction(const std::vector<double>& m_,
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
                            const double& w_,
                            const arma::mat& z_)
        : m(arma::conv_to<arma::vec>::from(m_)),
          d_yp(arma::conv_to<arma::vec>::from(d_yp_)),
          d_b0(arma::conv_to<arma::vec>::from(d_b0_)),
          d_bp(arma::conv_to<arma::vec>::from(d_bp_)),
          g_yp(arma::conv_to<arma::vec>::from(g_yp_)),
          g_b0(arma::conv_to<arma::vec>::from(g_b0_)),
          g_bp(arma::conv_to<arma::vec>::from(g_bp_)),
          L_0(arma::conv_to<arma::vec>::from(L_0_)),
          P_max(arma::conv_to<arma::vec>::from(P_max_)),
          u(u_),
          q(q_),
          W(arma::conv_to<arma::vec>::from(W_)),
          Phi(z_.n_rows, z_.n_cols),
          n_plants(z_.n_rows),
          weights(z_.n_rows),
          F(z_.n_rows),
          R(z_.n_rows) {
        fill_Phi__(w_, z_);
    };


    void make_weights(arma::vec& wts_vec,
                      const MatType& x) {

        if (F.n_elem != n_plants) F.set_size(n_plants);
        if (wts_vec.n_elem != n_plants) wts_vec.set_size(n_plants);

        for (size_t i = 0; i < n_plants; i++) {
            F(i) = x(i,0) + x(i,1) + x(i,2);
        }

        double wt_sum = 0;
        double YN_i;
        for (size_t i = 0; i < n_plants; i++) {
            wts_vec(i) = std::pow(F(i), q);
            YN_i = x(i,0) + x(i,2);
            if (F(i) > 0) YN_i /= F(i);
            wts_vec(i) *= std::pow(YN_i, u);
            wt_sum += wts_vec(i);
        }

        for (size_t i = 0; i < n_plants; i++) {
            if (wt_sum > 0 || W[i] > 0) {
                wts_vec[i] /= (wt_sum + W[i]);
            }
        }

        return;

    }



protected:


    arma::vec weights;
    arma::vec F;
    arma::vec R;

    /*
      Everything but calculating R, which differs by derived class.
     */
    void all_but_R(const MatType& x,
                   MatType& dxdt,
                   const double t) {

        const arma::vec Y(x.unsafe_col(0));
        const arma::vec B(x.unsafe_col(1));
        const arma::vec N(x.unsafe_col(2));

        arma::vec dYdt = dxdt.unsafe_col(0);
        arma::vec dBdt = dxdt.unsafe_col(1);
        arma::vec dNdt = dxdt.unsafe_col(2);

        // this should already be calculated inside make_weights:
        // this->F = Y + B + N;

        // to avoid dividing by zeros:
        arma::uvec non_zeros = arma::find(F > 0);

        arma::vec P = P_max % weights;

        arma::vec PF(n_plants, arma::fill::zeros);
        arma::vec YF(n_plants, arma::fill::zeros);
        arma::vec BF(n_plants, arma::fill::zeros);
        PF(non_zeros) = P(non_zeros) / F(non_zeros);
        YF(non_zeros) = Y(non_zeros) / F(non_zeros);
        BF(non_zeros) = B(non_zeros) / F(non_zeros);

        // L_0 is not allowed to be exactly zero, so this should always be okay:
        arma::vec Lambda = PF / (L_0 + PF);

        arma::vec gamma_y = g_yp % Lambda;
        arma::vec gamma_b = g_b0 + g_bp % Lambda;

        arma::vec delta_y = d_yp % Lambda;
        arma::vec delta_b = d_b0 + d_bp % Lambda;

        arma::vec growth_y = Phi * ((delta_y % YF + gamma_y) % N);
        arma::vec growth_b = Phi * ((delta_b % BF + gamma_b) % N);

        dYdt = growth_y - m % Y;
        dBdt = growth_b - m % B;
        dNdt = R - m % N - growth_y - growth_b;

        return;
    }

    // fill Phi matrix.
    // `z` should be n_plants x n_plants in size
    void fill_Phi__(const double& w_, const arma::mat& z_) {

        Phi.set_size(n_plants, n_plants);

        double col_sum;
        for (size_t j = 0; j < n_plants; j++) {
            col_sum = 0;
            for (size_t i = 0; i < n_plants; i++) {
                if (i == j) {
                    Phi(i,j) = 1;
                } else {
                    Phi(i,j) = std::exp(-w_ * z_(i, j));
                }
                col_sum += Phi(i,j);
            }
            for (size_t i = 0; i < n_plants; i++) Phi(i,j) /= col_sum;
        }

    }



};






#endif
