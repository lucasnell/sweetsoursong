# ifndef __SWEETSOURSONG_ODE_H
# define __SWEETSOURSONG_ODE_H


#include <Rcpp.h>
#include <vector>


// To avoid many warnings from BOOST
#pragma clang diagnostic ignored "-Wlanguage-extension-token"
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#pragma clang diagnostic warning "-Wlanguage-extension-token"


using namespace Rcpp;


typedef std::vector<double> VecType;

typedef boost::numeric::odeint::runge_kutta_dopri5<VecType> VecStepperType;

// Comment if you want to use boost matrix types:
typedef boost::numeric::ublas::matrix<double> MatType;

typedef boost::numeric::odeint::runge_kutta_dopri5<MatType> MatStepperType;

template <class C>
void len_check(bool& err,
               const std::vector<C>& vec,
               const std::string& vec_name,
               const size_t& req_len) {
    if (vec.size() != req_len) {
        Rcout << vec_name << " is length " << std::to_string(vec.size());
        Rcout << " but should be " << std::to_string(req_len) << "!" << std::endl;
        err = true;
    }
    return;
}



/*
 Weight pollinator movement when F is NOT constant.
 */
inline void landscape_weights__(std::vector<double>& wts_vec,
                                const MatType& x,
                                const size_t& n_plants,
                                const double& S_0,
                                const double& q,
                                const std::vector<double>& X,
                                const MatType& exp_wz) {

    std::vector<double> F;
    F.reserve(n_plants);
    for (size_t i = 0; i < n_plants; i++) {
        F.push_back(x(i,0) + x(i,1) + x(i,2));
    }

    if (wts_vec.size() != n_plants) wts_vec.resize(n_plants);

    double wt_F, wt_sum = 0;
    int power = 3;
    double pow_S_0 = std::pow(S_0, power);
    for (size_t i = 0; i < n_plants; i++) {
        wt_F = 0;
        for (size_t j = 0; j < n_plants; j++) {
            // Note: exp_wz(i,i) == 1
            wt_F += (F[j] * exp_wz(i,j));
        }
        double B_i = x(i,1);
        if (F[i] > 0) B_i /= F[i];
        wts_vec[i] = std::pow(wt_F, q) * (pow_S_0 / (pow_S_0 + std::pow(B_i, power)));
        wt_sum += wts_vec[i];
    }

    for (size_t i = 0; i < n_plants; i++) {
        if (wt_sum > 0 || X[i] > 0) {
            wts_vec[i] /= (wt_sum + std::pow(X[i], q));
        }
    }

    return;
}



#endif
