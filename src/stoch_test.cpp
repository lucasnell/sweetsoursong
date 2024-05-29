#include <RcppArmadillo.h>


using namespace Rcpp;

/*
 Testing stochastic simulations.
 */


/*
 libs/numeric/odeint/examples/stochastic_euler.hpp

 Copyright 2012 Karsten Ahnert
 Copyright 2012 Mario Mulansky

 Stochastic euler stepper example and Ornstein-Uhlenbeck process

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#include <vector>
// To avoid many warnings from BOOST
#pragma clang diagnostic ignored "-Wlanguage-extension-token"
#include <boost/random.hpp>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#pragma clang diagnostic warning "-Wlanguage-extension-token"



// //[ stochastic_euler_ornstein_uhlenbeck_def
// const static size_t N = 1;
// typedef boost::array< double , N > state_type;


typedef std::vector<double> state_type;



//[ stochastic_euler_class
// template< size_t N >
class stochastic_euler
{
public:

    // typedef boost::array< double , N > state_type;
    // typedef boost::array< double , N > deriv_type;
    // typedef double value_type;
    // typedef double time_type;
    // typedef unsigned short order_type;

    typedef boost::numeric::odeint::stepper_tag stepper_category;

    static unsigned short order( void ) { return 1; }

    template< class System >
    void do_step( System system , state_type &x , double t , double dt ) const
    {
        state_type det(x.size()) , stoch(x.size()) ;
        system.first( x , det );
        system.second( x , stoch );
        for( size_t i=0 ; i<x.size() ; ++i )
            x[i] += dt * det[i] + sqrt( dt ) * stoch[i];
    }
};
//]





struct ornstein_det
{
    void operator()( const state_type &x , state_type &dxdt ) const
    {
        dxdt[0] = -x[0];
    }
};

struct ornstein_stoch
{
    boost::mt11213b &m_rng;
    boost::normal_distribution<> m_dist;

    ornstein_stoch( boost::mt11213b &rng , double sigma ) : m_rng( rng ) , m_dist( 0.0 , sigma ) { }

    void operator()( const state_type &x , state_type &dxdt )
    {
        dxdt[0] = m_dist( m_rng );
    }
};
//]

// struct streaming_observer
// {
//     template< class State >
//     void operator()( const State &x , double t ) const
//     {
//         std::cout << t << "\t" << x[0] << "\n";
//     }
// };

struct vec_observer
{
    std::vector<state_type> data;
    std::vector<double> time;
    vec_observer() : data(), time() {};

    void operator()(const state_type& x, const double& t) {
        data.push_back(x);
        time.push_back(t);
        return;
    }
};



// int main( int argc , char **argv )
// {
//     using namespace std;
//     using namespace boost::numeric::odeint;
//
//     //[ ornstein_uhlenbeck_main
//     boost::mt11213b rng;
//     double dt = 0.1;
//     state_type x = {{ 1.0 }};
//     integrate_const( stochastic_euler< N >() ,
//                      std::make_pair( ornstein_det() , ornstein_stoch( rng , 1.0 ) ),
//                      x , 0.0 , 10.0 , dt , streaming_observer() );
//     //]
//     return 0;
// }

//' @export
// [[Rcpp::export]]
NumericMatrix stoch_test() {

    int32_t seed = static_cast<int32_t>(R::runif(0, 2147483647));

    boost::mt11213b rng;
    rng.seed(seed);
    double dt = 0.1;
    double max_t = 10.0;
    state_type x {1.0};
    vec_observer obs;
    boost::numeric::odeint::integrate_const(
        stochastic_euler() ,
        std::make_pair( ornstein_det() , ornstein_stoch( rng , 1.0 ) ),
        x , 0.0 , max_t , dt , std::ref(obs) );

    size_t n_steps = obs.data.size();
    NumericMatrix output(n_steps, x.size()+1U);
    colnames(output) = CharacterVector::create("t", "x");
    for (size_t i = 0; i < n_steps; i++) {
        output(i,0) = obs.time[i];
        for (size_t j = 0; j < x.size(); j++) {
            output(i,j+1U) = obs.data[i][j];
        }
    }
    return output;

}


