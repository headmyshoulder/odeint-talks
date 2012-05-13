/*
 * pendulum1.cpp
 *
 *  Created on: Apr 15, 2012
 *      Author: karsten
 */

#include <iostream>
#include <array>

#include <boost/numeric/odeint.hpp>

namespace odeint = boost::numeric::odeint;

typedef std::array< double , 2 > state_type;

struct pendulum
{
    double m_mu , m_omega , m_epsilon;

    pendulum( double mu , double omega , double epsilon )
    : m_mu( mu ) , m_omega( omega ) , m_epsilon( epsilon ) { }

    void operator()( const state_type &x , state_type &dxdt , double t ) const
    {
        dxdt[0] = x[1];
        dxdt[1] = - sin( x[0] ) - m_mu * x[1] + m_epsilon * sin( m_omega * t );
    }
};

int main( int argc , char **argv )
{
    odeint::runge_kutta4< state_type > rk4;

    // 0.0 , 0.0 , 0.0
    // 0.1 , 0.0 , 0.0
    // 0.1 , 1.05 , 1.5
    pendulum p( 0.1 , 1.05 , 1.5 );

    state_type x = {{ 1.0 , 0.0 }};
    double t = 0.0;
    const double dt = 0.025;

    std::cout << "unset key" << "\n";
    std::cout << "set size square" << std::endl;

    for( size_t i=0 ; i<2000 ; ++i )
    {
        std::cout << "p [-1.1:1.1][-1.1:1.1]  '-' w lp pt 7 ps 2 lw 3" << std::endl;
        std::cout << "0.0 0.0" << "\n";
        std::cout << -sin(x[0]) << " " << -cos(x[0]) << "\n";
        std::cout << "e" << std::endl;

        rk4.do_step( p , x , t , dt );
        t += dt;
    }


    return 0;
}
