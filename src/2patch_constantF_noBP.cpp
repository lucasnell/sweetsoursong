
/*
 Two plants where the total number of flowers (F) at each plant is
 constant and where pollinators are not affected by flower densities.
 We don't even include F and instead model Y and B as proportions
 of total flowers with non-colonized N = 1 - Y - B.
 Also, there is no effect of bacteria on pollinator visits and no immigration.
 */

#include <RcppArmadillo.h>
#include <vector>

#include "ode.h"

using namespace Rcpp;




class LandConstFnoBPSystemFunction
{
public:
    std::vector<double> m;
    std::vector<double> d_yp;
    std::vector<double> d_b0;
    std::vector<double> d_bp;
    std::vector<double> L_0;
    const size_t n_plants = 2U;


    LandConstFnoBPSystemFunction(const std::vector<double>& m_,
                            const std::vector<double>& d_yp_,
                            const std::vector<double>& d_b0_,
                            const std::vector<double>& d_bp_,
                            const std::vector<double>& L_0_)
        : m(m_),
          d_yp(d_yp_),
          d_b0(d_b0_),
          d_bp(d_bp_),
          L_0(L_0_) {};

    void operator()(const MatType& x,
                    MatType& dxdt,
                    const double t) {

        for (size_t i = 0; i < n_plants; i++) {

            const double& Y(x(i,0));
            const double& B(x(i,1));
            double N = 1 - Y - B;

            double& dYdt(dxdt(i,0));
            double& dBdt(dxdt(i,1));

            double P = 0.5;

            double Lambda = P / (L_0[i] + P);

            dYdt = d_yp[i] * Lambda * Y * N - m[i] * Y;
            dBdt = (d_b0[i] + d_bp[i] * Lambda) * B * N - m[i] * B;
        }

        return;
    }




};






//' @export
// [[Rcpp::export]]
NumericMatrix twopatch_constantF_noBP_ode(const std::vector<double>& m,
                                          const std::vector<double>& d_yp,
                                          const std::vector<double>& d_b0,
                                          const std::vector<double>& d_bp,
                                          const std::vector<double>& L_0,
                                          const std::vector<double>& Y0,
                                          const std::vector<double>& B0,
                                          const double& dt = 0.1,
                                          const double& max_t = 90.0) {

    size_t np = 2U;
    /*
     I can't just use 'stop()' because it causes a segfault (or similar).
     The system below is my workaround.
     */
    bool err = false;
    len_check(err, m, "m", np);
    len_check(err, d_yp, "d_yp", np);
    len_check(err, d_b0, "d_b0", np);
    len_check(err, d_bp, "d_bp", np);
    len_check(err, L_0, "L_0", np);
    len_check(err, Y0, "Y0", np);
    len_check(err, B0, "B0", np);

    if (err) return NumericMatrix(0,0);

    size_t n_states = 2U;
    MatType x(np, n_states);
    for (size_t i = 0; i < np; i++) {
        x(i,0) = Y0[i];
        x(i,1) = B0[i];
    }


    Observer<MatType> obs;
    LandConstFnoBPSystemFunction system(m, d_yp, d_b0, d_bp, L_0);

    boost::numeric::odeint::integrate_const(
        MatStepperType(), std::ref(system),
        x, 0.0, max_t, dt, std::ref(obs));

    size_t n_steps = obs.data.size();
    NumericMatrix output(n_steps * np, n_states+2U);
    colnames(output) = CharacterVector::create("t", "p", "Y", "B");
    std::vector<double> wts(np);
    size_t i = 0;
    for (size_t t = 0; t < n_steps; t++) {
        for (size_t k = 0; k < np; k++) {
            output(i,0) = obs.time[t];
            output(i,1) = k;
            for (size_t j = 0; j < n_states; j++) {
                output(i,j+2U) = obs.data[t](k,j);
            }
            i++;
        }

    }
    return output;
}
