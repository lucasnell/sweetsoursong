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


#endif
