#include <iostream>
#include <vector>
#include <utility>
#include <tuple>
#include <random>
#include <algorithm>

#include <vexcl/vexcl.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/external/vexcl/vexcl_resize.hpp>

namespace odeint = boost::numeric::odeint;

typedef vex::vector< double > vector_type;
typedef vex::multivector< double, 2 > state_type;

struct ensemble
{
    const vector_type &m_eps;
    const vector_type &m_omega;
    double m_mu;

    ensemble( const vector_type &eps , const vector_type &omega , double mu )
        : m_eps( eps ) , m_omega( omega ) , m_mu( mu )
    {
    }

    void operator()( const state_type &x , state_type &dxdt , double t )
    {
        dxdt = std::make_tuple( x(1) , - sin( x(0) ) - m_mu * x(1) + m_eps * sin( m_omega * t ) );
    }
};


const double dt = 0.01;
const double t_max = 100.0;

int main( int argc , char **argv )
{
    const size_t n_eps = 132;
    const size_t n_omega = 132;
    const size_t n = n_eps * n_omega;

    std::vector<double> x( 2 * n );
    std::fill( x.begin() , x.end() , 1.0 );

    std::vector< double > eps( n ) , omega( n );
    double eps_start = 0.0 , eps_end = 5.0;
    double omega_start = 0.5 , omega_end = 1.5;
    double e = eps_start , de = ( eps_end - eps_start ) / double( n_eps - 1 );
    for( size_t i=0 ; i<n_eps ; ++i,e+=de )
    {
        double o = omega_start , dom = ( omega_end - omega_start ) / double( n_omega - 1 );
        for( size_t j=0 ; j<n_omega ; ++j,o+=dom )
        {
            eps  [ i * n_omega + j ] = e;
            omega[ i * n_omega + j ] = o;
        }
    }

    vex::Context ctx( vex::Filter::Exclusive( vex::Filter::DoublePrecision && vex::Filter::Env ) );
    std::cout << ctx << std::endl;

    

    state_type X(ctx.queue(), n);
    vex::copy( x.begin() , x.begin() + n, X(0).begin() );
    vex::copy( x.begin() + n, x.end() , X(1).begin() );

    vector_type Eps( ctx.queue() , n );
    vector_type Omega( ctx.queue() , n );
    vex::copy( eps.begin() , eps.end() , Eps.begin() );
    vex::copy( omega.begin() , omega.end() , Omega.begin() );


    odeint::runge_kutta4<
        state_type , double , state_type , double ,
        odeint::vector_space_algebra , odeint::default_operations
        > stepper;

    odeint::integrate_const( stepper , ensemble( Eps , Omega , 0.1 )
                             , X , double(0.0) , t_max , dt );

    std::vector< double > res( 2 * n );
    vex::copy( X(0) , res );

    for( size_t i=0 ; i<n ; ++i )
        std::cout << eps[i] << "\t" << omega[i] << "\t" << res[i] << "\n";

}
