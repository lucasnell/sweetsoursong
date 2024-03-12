# ifndef __SWEETSOURSONG_MATH_H
# define __SWEETSOURSONG_MATH_H


#include <Rcpp.h>
#include <cmath>


using namespace Rcpp;


/*
 ---------
 Different parameterizations of the PDF for the Gamma distribution:
 ---------
 */
inline double gamma_pdf__(double x, double shape, double rate) {
    double f = std::pow(x, shape-1) * std::exp(- rate * x);
    f *= std::pow(rate, shape);
    f /= std::tgamma(shape);
    return f;
}
inline double gamma_pdf2__(double x, double mean, double variance) {
    double shape = mean * mean / variance;
    double rate = mean / variance;
    double f = std::pow(x, shape-1) * std::exp(- rate * x);
    f *= std::pow(rate, shape);
    f /= std::tgamma(shape);
    return f;
}
inline double gamma_pdf3__(double x, double mean, double skew) {
    double shape = 4 / (skew * skew);
    double rate = 4 / (mean * skew * skew);
    double f = std::pow(x, shape-1) * std::exp(- rate * x);
    f *= std::pow(rate, shape);
    f /= std::tgamma(shape);
    return f;
}
inline double gamma_pdf4__(double x, double mode, double skew) {
    double shape = 4 / (skew * skew);
    double rate = (shape - 1) / mode;
    double f = std::pow(x, shape-1) * std::exp(- rate * x);
    f *= std::pow(rate, shape);
    f /= std::tgamma(shape);
    return f;
}


#endif


