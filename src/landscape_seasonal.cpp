#define _USE_MATH_DEFINES


#include <RcppArmadillo.h>
#include <vector>
#include <cmath>

#include "landscape.h"

using namespace Rcpp;





class SeasonalLandscape : public LandscapeSystemFunction
{
public:
    arma::vec R_hat;
    arma::vec par1;
    arma::vec par2;
    std::vector<char> distr_types;

    SeasonalLandscape(const std::vector<double>& m_,
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
                      const arma::mat& z_,
                      const double& min_F_for_P_,
                      const std::vector<double>& R_hat_,
                      const std::vector<double>& par1_,
                      const std::vector<double>& par2_,
                      const std::vector<char>& distr_types_,
                      const std::vector<double>& Y0_,
                      const std::vector<double>& B0_,
                      const double& add_F_)
        : LandscapeSystemFunction(m_, d_yp_, d_b0_, d_bp_, g_yp_, g_b0_, g_bp_,
                                  L_0_, P_max_, u_, q_, W_, w_, z_, min_F_for_P_),
          R_hat(arma::conv_to<arma::vec>::from(R_hat_)),
          par1(arma::conv_to<arma::vec>::from(par1_)),
          par2(arma::conv_to<arma::vec>::from(par2_)),
          distr_types(distr_types_),
          Y0(arma::conv_to<arma::vec>::from(Y0_)),
          B0(arma::conv_to<arma::vec>::from(B0_)),
          add_F(add_F_),
          YB_added(z_.n_rows, false) {

        for (size_t i = 0; i < n_plants; i++) {
            // No reason to add these if these are set to zero:
            if ((Y0(i) + B0(i)) == 0) YB_added[i] = true;
        }

    };


    void operator()(const MatType& x,
                    MatType& dxdt,
                    const double t) {

        LandscapeSystemFunction::make_weights(this->weights, x);
        make_R(t);
        LandscapeSystemFunction::all_but_R(x, dxdt, t);

        for (size_t i = 0; i < n_plants; i++) {
            if (F(i) >= add_F && ! YB_added[i]) {
                const double& N(x(i,2));
                double& dYdt = dxdt(i,0);
                double& dBdt = dxdt(i,1);
                double& dNdt = dxdt(i,2);
                // note: using 'N + dNdt' below to avoid N going < 0
                double y0_i = (N + dNdt) * Y0(i);
                double b0_i = (N + dNdt) * B0(i);
                dYdt += y0_i;
                dBdt += b0_i;
                dNdt -= (y0_i + b0_i);
                YB_added[i] = true;
            }
        }


        return;
    }

    void make_weights(arma::vec& wts_vec,
                      const MatType& x) {
        LandscapeSystemFunction::make_weights(wts_vec, x);
    }



private:
    arma::vec Y0;
    arma::vec B0;
    double add_F;
    std::vector<bool> YB_added;
    const double sqrt_2pi = std::sqrt(2 * M_PI);


    void make_R(const double& t) {

        for (size_t i = 0; i < n_plants; i++) {
            const char& dtype(distr_types[i]);
            switch (dtype) {
            case 'N':
                // normal
                normal_R(t, i);
                break;
            case 'W':
                // weibull
                weibull_R(t, i);
                break;
            case 'L':
                // lognormal
                lognormal_R(t, i);
                break;
            default:
                // revert to normal
                normal_R(t, i);
                break;
            }
        }
    }


    inline void normal_R(const double& t, const size_t& i) {
        const double& mu(par1(i));
        const double& sigma(par2(i));
        double tmp = (t - mu) / sigma;
        R(i) = (R_hat(i) / (sigma * sqrt_2pi)) *
            std::exp(-0.5 * (tmp*tmp));
        return;
    }

    inline void weibull_R(const double& t, const size_t& i) {
        const double& lambda(par1(i));
        const double& k(par2(i));
        R(i) = R_hat(i) * (k / lambda) * std::pow(t / lambda, k-1) *
            std::exp(- std::pow(t / lambda, k));
        return;
    }

    inline void lognormal_R(const double& t, const size_t& i) {
        const double& mu(par1(i));
        const double& sigma(par2(i));
        double tmp = std::log(t) - mu;
        R(i) = (R_hat(i) / (t * sigma * sqrt_2pi)) *
            std::exp(- (tmp*tmp) / (2 * sigma * sigma));
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
                                   const std::vector<double>& par1,
                                   const std::vector<double>& par2,
                                   const StringVector& distr_types,
                                   const double& w,
                                   const arma::mat& z,
                                   const double& min_F_for_P,
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
    bool err = lanscape_arg_checks(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp,
                                   L_0, P_max, u, q, W, w, z, min_F_for_P,
                                   Y0, B0, dt, max_t);
    len_check(err, R_hat, "R_hat", np);
    len_check(err, par1, "par1", np);
    len_check(err, par2, "par2", np);
    len_check(err, distr_types, "distr_types", np);
    min_val_check(err, R_hat, "R_hat", 0, false);
    min_val_check(err, par1, "par1", 0, false);
    min_val_check(err, par2, "par2", 0, false);
    std::vector<char> distr_types_char;
    distr_types_char.reserve(np);
    std::string d;
    for (size_t i = 0; i < distr_types.size(); i++) {
        d = distr_types(i);
        if (d != "N" && d != "W" && d != "L") {
            Rcout << "distr_types must only contain 'N', 'W', or 'L'. ";
            Rcout << "Yours contains at least one '" << d << "'." << std::endl;
            err = true;
            break;
        }
        distr_types_char.push_back(d[0]);
    }
    min_val_check(err, add_F, "add_F", 0, false);
    for (size_t i = 0; i < std::min(B0.size(), Y0.size()); i++) {
        if ((Y0[i] + B0[i]) > add_F) {
            Rcout << "Y0+B0 must always be <= `add_F`." << std::endl;
            err = true;
            break;
        }
    }
    if (err) return NumericMatrix(0,0);


    MatType x(np, 3);
    for (size_t i = 0; i < np; i++) {
        x(i,0) = 0.0;  // Y
        x(i,1) = 0.0;  // B
        x(i,2) = 0.0;  // N
    }


    Observer<MatType> obs;
    SeasonalLandscape system(m, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, P_max,
                             u, q, W, w, z, min_F_for_P,
                             R_hat, par1, par2, distr_types_char,
                             Y0, B0, add_F);

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
