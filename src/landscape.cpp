

#include <RcppArmadillo.h>
#include <vector>

#include "landscape.h"


using namespace Rcpp;





class NonSeasonalLandscape : public LandscapeSystemFunction
{
public:

    NonSeasonalLandscape(const std::vector<double>& m_,
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
                         const double& a_,
                         const arma::mat& z_,
                         const double& min_F_for_P_,
                         const std::vector<double>& R_)
        : LandscapeSystemFunction(m_, d_yp_, d_b0_, d_bp_, g_yp_, g_b0_, g_bp_,
                                  L_0_, P_max_, u_, q_, W_, a_, z_,
                                  min_F_for_P_) {

        this->R = arma::conv_to<arma::vec>::from(R_);

    };


    void operator()(const MatType& x,
                    MatType& dxdt,
                    const double t) {

        LandscapeSystemFunction::make_weights(this->weights, x);
        LandscapeSystemFunction::all_but_R(x, dxdt, t);
        return;

    }

    void make_weights(arma::vec& wts_vec,
                      const MatType& x) {
        LandscapeSystemFunction::make_weights(wts_vec, x);
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
                            const double& u,
                            const double& q,
                            const std::vector<double>& W,
                            const double& a,
                            const arma::mat& z,
                            const double& min_F_for_P,
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
    bool err = lanscape_arg_checks(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                                   L_0, P_max, u, q, W, a, z, min_F_for_P,
                                   Y0, B0, dt, max_t);
    len_check(err, R, "R", np);
    len_check(err, N0, "N0", np);
    min_val_check(err, R, "R", 0);
    min_val_check(err, N0, "N0", 0);
    if (err) return NumericMatrix(0,0);

    MatType x(np, 3);
    for (size_t i = 0; i < np; i++) {
        x(i,0) = Y0[i];
        x(i,1) = B0[i];
        x(i,2) = N0[i];
    }

    Observer<MatType> obs;

    NonSeasonalLandscape system(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                                L_0, P_max, u, q, W, a, z, min_F_for_P, R);

    boost::numeric::odeint::integrate_const(
        MatStepperType(), std::ref(system),
        x, 0.0, max_t, dt, std::ref(obs));

    size_t n_steps = obs.data.size();
    NumericMatrix output(n_steps * np, 6);
    colnames(output) = CharacterVector::create("t", "p", "Y", "B", "N", "P");
    arma::vec wts(np);
    size_t i = 0;
    for (size_t t = 0; t < n_steps; t++) {
        system.make_weights(wts, obs.data[t]);
        for (size_t k = 0; k < np; k++) {
            output(i,0) = obs.time[t];
            output(i,1) = k;
            output(i,2) = obs.data[t](k, 0);
            output(i,3) = obs.data[t](k, 1);
            output(i,4) = obs.data[t](k, 2);
            output(i,5) = P_max[k] * wts(k);
            i++;
        }

    }
    return output;
}
