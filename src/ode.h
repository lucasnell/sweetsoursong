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

/*
 This function checks for whether the remainder of numer / denom matches
 one of the remainders inside the rmd vector.
 The use of 1e-10 is to avoid rounding issues with std::fmod.
 This function only works if numer > denom.
 */
inline bool remainder_equals(const double& numer,
                             const double& denom,
                             const std::vector<double>& rmd) {
    double nd_rmd = std::fmod(numer, denom);
    for (const double& r : rmd) {
        if (std::abs(nd_rmd - r) < 1e-10) return true;
    }
    return false;
}

/*
 This is similar to above but only checks for zero remainder
 The use of 1e-10 is to avoid rounding issues with std::remainder.
 This function works when numer < denom.
 */
inline bool zero_remainder(const double& numer, const double& denom) {
    return std::abs(std::remainder(numer, denom)) < 1e-10;
}




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






/*
 ==============================================================================
 ==============================================================================
 Observer template classes
 ==============================================================================
 ==============================================================================
 */


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

    void clear() {
        data.clear();
        time.clear();
    }
};

/*
 Same as above but for simulations with a burn-in period and not saving every
 time step.
 This includes just saving begin and end of each season.
 */
template< class C >
struct ObserverBurnEvery
{
    std::vector<C> data;
    std::vector<double> time;
    double burnin;
    size_t save_every;
    std::vector<double> remainders;
    ObserverBurnEvery(const double& burnin_, const size_t& save_every_,
                      const size_t& season_len_, const bool& begin_end_)
        : data(), time(), burnin(burnin_), save_every(save_every_),
          season_len(season_len_), begin_end(begin_end_),
          iters(save_every_) {};
    // note: setting `iters` to `save_every_` so that the first time step
    // is always included

    void operator()(const C& x, const double& t) {
        if (t >= burnin) {
            if (begin_end) {
                if (iters >= season_len) {
                    push_back__(x, t);
                    iters = 0;
                } else if (iters == 1U) {
                    push_back__(x, t);
                }
            } else if (iters >= save_every) {
                push_back__(x, t);
                iters = 0;
            }
            // Don't start iterating until after burn-in:
            iters++;
        }
        return;
    }

    void clear() {
        data.clear();
        time.clear();
    }

protected:

    size_t season_len;
    bool begin_end;
    size_t iters;

private:

    void push_back__(const C& x, const double& t) {
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
