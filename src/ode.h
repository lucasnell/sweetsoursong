# ifndef __SWEETSOURSONG_ODE_H
# define __SWEETSOURSONG_ODE_H


#include <RcppArmadillo.h>
#include <vector>
#include <string>


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
// typedef boost::numeric::ublas::matrix<double> MatType;
typedef arma::mat MatType;

typedef boost::numeric::odeint::runge_kutta_dopri5<MatType> MatStepperType;


template< class C >
struct Observer
{
    std::vector<C> data;
    std::vector<double> time;
    Observer() : data(), time() {};

    void operator()(const C& x, const double& t) {
        data.push_back(x);
        time.push_back(t);
        return;
    }
};




/*
 ==============================================================================
 ==============================================================================
 Checking functions
 ==============================================================================
 ==============================================================================
 */


inline void len_check(bool& err,
                      const std::vector<double>& vec,
                      const std::string& vec_name,
                      const size_t& req_len) {
    if (vec.size() != req_len) {
        Rcout << vec_name << " is length " << std::to_string(vec.size());
        Rcout << " but should be " << std::to_string(req_len) << "!" << std::endl;
        err = true;
    }
    return;
}
inline void len_check(bool& err,
                      const StringVector& vec,
                      const std::string& vec_name,
                      const size_t& req_len) {
    if (vec.size() != req_len) {
        Rcout << vec_name << " is length " << std::to_string(vec.size());
        Rcout << " but should be " << std::to_string(req_len) << "!" << std::endl;
        err = true;
    }
    return;
}

// check square matrix dimensions:
inline void mat_dim_check(bool& err,
                          const MatType& mat,
                          const std::string& mat_name,
                          const size_t& req_rows_cols) {
    if (mat.n_rows != req_rows_cols || mat.n_cols != req_rows_cols) {
        std::string r = std::to_string(req_rows_cols);
        Rcout << mat_name << " is " << std::to_string(mat.n_rows) << " x ";
        Rcout << std::to_string(mat.n_cols) << " (rows x columns), but ";
        Rcout << "should be " << r << " x " << r << "!" << std::endl;
        err = true;
    }
    return;
}

// check values:
inline void min_val_check(bool& err,
                          const double& obj,
                          const std::string& obj_name,
                          const double& min_val,
                          const bool& allow_equal = true) {
    if (allow_equal) {
        if (obj < min_val) {
            Rcout << obj_name << " is " << std::to_string(obj);
            Rcout << " but should be >= " << std::to_string(min_val) << "!" << std::endl;
            err = true;
        }
    } else {
        if (obj <= min_val) {
            Rcout << obj_name << " is " << std::to_string(obj);
            Rcout << " but should be > " << std::to_string(min_val) << "!" << std::endl;
            err = true;
        }
    }
    return;
}
// overloaded for vector
inline void min_val_check(bool& err,
                          const std::vector<double>& obj,
                          const std::string& obj_name,
                          const double& min_val,
                          const bool& allow_equal = true) {
    double obj_min = *std::min_element(obj.begin(), obj.end());
    if (allow_equal) {
        if (obj_min < min_val) {
            Rcout << obj_name << " has a minimum value of " << std::to_string(obj_min);
            Rcout << " but should be >= " << std::to_string(min_val) << "!" << std::endl;
            err = true;
        }
    } else {
        if (obj_min <= min_val) {
            Rcout << obj_name << " has a minimum value of " << std::to_string(obj_min);
            Rcout << " but should be > " << std::to_string(min_val) << "!" << std::endl;
            err = true;
        }
    }
    return;
}
// overloaded for matrix
inline void min_val_check(bool& err,
                          const arma::mat& obj,
                          const std::string& obj_name,
                          const double& min_val,
                          const bool& allow_equal = true) {
    double obj_min = obj.min();
    if (allow_equal) {
        if (obj_min < min_val) {
            Rcout << obj_name << " has a minimum value of " << std::to_string(obj_min);
            Rcout << " but should be >= " << std::to_string(min_val) << "!" << std::endl;
            err = true;
        }
    } else {
        if (obj_min <= min_val) {
            Rcout << obj_name << " has a minimum value of " << std::to_string(obj_min);
            Rcout << " but should be > " << std::to_string(min_val) << "!" << std::endl;
            err = true;
        }
    }
    return;
}







#endif
