#include <RcppArmadillo.h>
#include <vector>

#include "ode.h"

using namespace Rcpp;




class LocalSystemFunction
{
public:
    size_t n_Y;
    size_t n_B;
    double Y_delay;
    double Y0;
    double B_delay;
    double B0;
    double D;
    double A_0;
    std::vector<double> r_Y;
    std::vector<double> m_Y;
    std::vector<double> q_Y;
    std::vector<double> c_Y;
    std::vector<double> h_Y;
    std::vector<double> r_B;
    std::vector<double> m_B;
    std::vector<double> e_B;
    std::vector<double> q_B;
    std::vector<double> c_B;
    std::vector<double> h_B;

    LocalSystemFunction(const size_t& n_Y_,
                   const size_t& n_B_,
                   const double& Y_delay_,
                   const double& Y0_,
                   const double& B_delay_,
                   const double& B0_,
                   const double& D_,
                   const double& A_0_,
                   const double& r_Y_,
                   const double& m_Y_,
                   const double& q_Y_,
                   const double& c_Y_,
                   const double& h_Y_,
                   const double& r_B_,
                   const double& m_B_,
                   const double& e_B_,
                   const double& q_B_,
                   const double& c_B_,
                   const double& h_B_)
        : n_Y(n_Y_),
          n_B(n_B_),
          Y_delay(Y_delay_),
          Y0(Y0_),
          B_delay(B_delay_),
          B0(B0_),
          D(D_),
          A_0(A_0_),
          r_Y(n_Y, r_Y_),
          m_Y(n_Y, m_Y_),
          q_Y(n_Y, q_Y_),
          c_Y(n_Y, c_Y_),
          h_Y(n_Y, h_Y_),
          r_B(n_B, r_B_),
          m_B(n_B, m_B_),
          e_B(n_B, e_B_),
          q_B(n_B, q_B_),
          c_B(n_B, c_B_),
          h_B(n_B, h_B_) {};
    LocalSystemFunction(const size_t& n_Y_,
                   const size_t& n_B_,
                   const double& Y_delay_,
                   const double& Y0_,
                   const double& B_delay_,
                   const double& B0_,
                   const double& D_,
                   const double& A_0_,
                   const std::vector<double>& r_Y_,
                   const std::vector<double>& m_Y_,
                   const std::vector<double>& q_Y_,
                   const std::vector<double>& c_Y_,
                   const std::vector<double>& h_Y_,
                   const std::vector<double>& r_B_,
                   const std::vector<double>& m_B_,
                   const std::vector<double>& e_B_,
                   const std::vector<double>& q_B_,
                   const std::vector<double>& c_B_,
                   const std::vector<double>& h_B_)
        : n_Y(n_Y_),
          n_B(n_B_),
          Y_delay(Y_delay_),
          Y0(Y0_),
          B_delay(B_delay_),
          B0(B0_),
          D(D_),
          A_0(A_0_),
          r_Y(r_Y_),
          m_Y(m_Y_),
          q_Y(q_Y_),
          c_Y(c_Y_),
          h_Y(h_Y_),
          r_B(r_B_),
          m_B(m_B_),
          e_B(e_B_),
          q_B(q_B_),
          c_B(c_B_),
          h_B(h_B_) {
        if (r_Y.size() != n_Y) Rcpp::stop("r_Y.size() != n_Y");
        if (m_Y.size() != n_Y) Rcpp::stop("m_Y.size() != n_Y");
        if (q_Y.size() != n_Y) Rcpp::stop("q_Y.size() != n_Y");
        if (c_Y.size() != n_Y) Rcpp::stop("c_Y.size() != n_Y");
        if (h_Y.size() != n_Y) Rcpp::stop("h_Y.size() != n_Y");
        if (r_B.size() != n_B) Rcpp::stop("r_B.size() != n_B");
        if (m_B.size() != n_B) Rcpp::stop("m_B.size() != n_B");
        if (e_B.size() != n_B) Rcpp::stop("e_B.size() != n_B");
        if (q_B.size() != n_B) Rcpp::stop("q_B.size() != n_B");
        if (c_B.size() != n_B) Rcpp::stop("c_B.size() != n_B");
        if (h_B.size() != n_B) Rcpp::stop("h_B.size() != n_B");
    };
    void operator()(const VecType& x, VecType& dxdt, const double t) {
        size_t n_YB = n_Y+n_B;
        const double& A(x[n_YB]);
        const double& H(x[n_YB+1U]);
        double& dAdt(dxdt[n_YB]);
        double& dHdt(dxdt[n_YB+1U]);
        dAdt = D * (A_0 - A);
        dHdt = - D * H;
        double Ryi, Rbj;

        bool add_B = B_delay > 0 && t >= B_delay;
        if (add_B) B_delay = -1;
        bool add_Y = Y_delay > 0 && t >= Y_delay;
        if (add_Y) Y_delay = -1;

        for (size_t i = 0; i < n_Y; i++) {
            const double& Yi(x[i]);
            if (Yi == 0) {
                if (add_Y) dxdt[i] = Y0;
                continue;
            }
            Ryi = r_Y[i] * A / ((c_Y[i] + A) * (1 + H / h_Y[i]));
            // Ryi = r_Y[i] * A * A / ((c_Y[i] + A * A) * (1 + H / h_Y[i]));
            dxdt[i] = Yi * (Ryi - m_Y[i]);
            dAdt -= (q_Y[i] * Yi * Ryi);
        }
        for (size_t j = 0; j < n_B; j++) {
            const double& Bj(x[j+n_Y]);
            if (Bj == 0) {
                if (add_B) dxdt[j+n_Y] = B0;
                continue;
            }
            Rbj = r_B[j] * A / ((c_B[j] + A) * (1 + H / h_B[j]));
            // Rbj = r_B[j] * A * A / (c_B[j] + A * A);
            dxdt[j+n_Y] = Bj * (Rbj - m_B[j]);
            dAdt -= (q_B[j] * Bj * Rbj);
            dHdt += (e_B[j] * Bj * Rbj);
        }
        return;
    }
};



//' @export
 // [[Rcpp::export]]
 NumericMatrix run_ode_cpp(const double& dt = 0.01,
                           const double& max_t = 36.0,
                           const double& Y_delay = 0,
                           const double& B_delay = 0,
                           const double& Y0 = 1.0,
                           const double& B0 = 1.0,
                           const double& A0 = 1.46,
                           const double& H0 = 0.0,
                           const double& D = 0.214,
                           double A_0 = -999,
                           const double& r_Y = 0.44,
                           const double& r_B = 0.264,
                           const double& m_Y = 0.01,
                           const double& m_B = 0.01,
                           const double& e_B = 0.84,
                           const double& q_Y = 0.022,
                           const double& q_B = 0.0,
                           const double& c_Y = 0.152,
                           const double& c_B = 1,
                           const double& h_B = 0.124,
                           const double& h_Y = 0.044) {

     if (A_0 == -999) A_0 = A0;

     VecType x(4U, 0.0);
     if (Y_delay <= 0) x[0] = Y0;
     if (B_delay <= 0) x[1] = B0;
     x[2] = A0;
     x[3] = H0;

     Observer<VecType> obs;
     LocalSystemFunction system(1U, 1U, Y_delay, Y0, B_delay, B0, D, A_0, r_Y, m_Y, q_Y,
                           c_Y, h_Y, r_B, m_B, e_B, q_B, c_B, h_B);

     boost::numeric::odeint::integrate_const(
         VecStepperType(), std::ref(system),
         x, 0.0, max_t, dt, std::ref(obs));

     size_t n_steps = obs.data.size();
     NumericMatrix output(n_steps, x.size()+1U);
     colnames(output) = CharacterVector::create("t", "Y", "B", "A", "H");
     for (size_t i = 0; i < n_steps; i++) {
         output(i,0) = obs.time[i];
         for (size_t j = 0; j < x.size(); j++) {
             output(i,j+1U) = obs.data[i][j];
         }
     }
     return output;
 }
