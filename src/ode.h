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
#include <boost/random.hpp>
#pragma clang diagnostic warning "-Wlanguage-extension-token"


using namespace Rcpp;


typedef std::vector<double> VecType;

typedef boost::numeric::odeint::runge_kutta_dopri5<VecType> VecStepperType;

// Comment if you want to use boost matrix types:
// typedef boost::numeric::ublas::matrix<double> MatType;
typedef arma::mat MatType;

typedef boost::numeric::odeint::runge_kutta_dopri5<MatType> MatStepperType;


// From https://stackoverflow.com/a/41820991
namespace boost { namespace numeric { namespace odeint {

template <>
struct is_resizeable<arma::vec>
{
    typedef boost::true_type type;
    const static bool value = type::value;
};

template <>
struct same_size_impl<arma::vec, arma::vec>
{
    static bool same_size(const arma::vec& x, const arma::vec& y)
    {
        return x.n_elem == y.n_elem;
    }
};

template<>
struct resize_impl<arma::vec, arma::vec>
{
    static void resize(arma::vec &v1, const arma::vec& v2)
    {
        v1.resize(v2.n_elem);
    }
};

template <>
struct is_resizeable<arma::mat>
{
    typedef boost::true_type type;
    const static bool value = type::value;
};

template <>
struct same_size_impl<arma::mat, arma::mat>
{
    static bool same_size(const arma::mat& x, const arma::mat& y)
    {
        return x.n_rows == y.n_rows && x.n_cols == y.n_cols;
    }
};

template<>
struct resize_impl<arma::mat, arma::mat>
{
    static void resize(arma::mat &v1, const arma::mat& v2)
    {
        v1.resize(v2.n_rows, v2.n_cols);
    }
};

} } } // namespace boost::numeric::odeint


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


// ------------------------------------------------------------
// ------------------------------------------------------------
// Same but for max
inline void max_val_check(bool& err,
                          const double& obj,
                          const std::string& obj_name,
                          const double& max_val,
                          const bool& allow_equal = true) {
    if (allow_equal) {
        if (obj > max_val) {
            Rcout << obj_name << " is " << std::to_string(obj);
            Rcout << " but should be <= " << std::to_string(max_val) << "!" << std::endl;
            err = true;
        }
    } else {
        if (obj >= max_val) {
            Rcout << obj_name << " is " << std::to_string(obj);
            Rcout << " but should be < " << std::to_string(max_val) << "!" << std::endl;
            err = true;
        }
    }
    return;
}
// overloaded for vector
inline void max_val_check(bool& err,
                          const std::vector<double>& obj,
                          const std::string& obj_name,
                          const double& max_val,
                          const bool& allow_equal = true) {
    double obj_max = *std::max_element(obj.begin(), obj.end());
    if (allow_equal) {
        if (obj_max > max_val) {
            Rcout << obj_name << " has a maximum value of " << std::to_string(obj_max);
            Rcout << " but should be <= " << std::to_string(max_val) << "!" << std::endl;
            err = true;
        }
    } else {
        if (obj_max >= max_val) {
            Rcout << obj_name << " has a maximum value of " << std::to_string(obj_max);
            Rcout << " but should be < " << std::to_string(max_val) << "!" << std::endl;
            err = true;
        }
    }
    return;
}
// overloaded for matrix
inline void max_val_check(bool& err,
                          const arma::mat& obj,
                          const std::string& obj_name,
                          const double& max_val,
                          const bool& allow_equal = true) {
    double obj_max = obj.max();
    if (allow_equal) {
        if (obj_max > max_val) {
            Rcout << obj_name << " has a maximum value of " << std::to_string(obj_max);
            Rcout << " but should be <= " << std::to_string(max_val) << "!" << std::endl;
            err = true;
        }
    } else {
        if (obj_max >= max_val) {
            Rcout << obj_name << " has a maximum value of " << std::to_string(obj_max);
            Rcout << " but should be < " << std::to_string(max_val) << "!" << std::endl;
            err = true;
        }
    }
    return;
}







#endif
