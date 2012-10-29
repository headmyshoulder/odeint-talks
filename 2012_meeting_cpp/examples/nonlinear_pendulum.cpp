/*
 * duffing.cpp
 *
 *  Created on: Apr 15, 2012
 *      Author: karsten
 */

#include <iostream>
#include <array>

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix.hpp>

typedef std::array< double , 2 > state_type;

class pendulum
{
public:

    pendulum( double mu , double omega , double epsilon )
    : m_mu( mu ) , m_omega( omega ) , m_epsilon( epsilon ) { }

    void operator()( const state_type &x , state_type &dxdt , double t ) const
    {
        dxdt[0] = x[1];
        dxdt[1] = - sin( x[0] ) - m_mu * x[1] + m_epsilon * sin( m_omega * t );
    }

private:

    double m_mu;
    double m_omega;
    double m_epsilon;
};

struct write_for_gnuplot
{
    void operator()( const state_type &x , double t ) const
    {
        using namespace std;
        cout << "unset key" << endl;
        cout << "set size square" << endl;
        cout << "p [-1.1:1.1][-1.1:1.1]  '-' w lp pt 7 ps 2 lw 3" << endl;
        cout << "0.0 0.0" << "\n";
        cout << -sin(x[0]) << " " << -cos(x[0]) << "\n";
        cout << "e" << endl;
    }
};

int main( int argc , char **argv )
{
    using namespace boost::numeric::odeint;
    using namespace boost::phoenix::arg_names;
    state_type x = {{ 1.0 , 0.0 }};
    pendulum d( 0.1 , 1.05 , 0.0 );
//    duffing d( 0.1 , 1.05 , 1.5 );
//    integrate_const( runge_kutta4< state_type >() , d , x , 0.0 , 1000.0 , 0.1 ,
//            std::cout << arg2 << "\t" << arg1[0] << "\t" << arg1[1] << "\n" );

    integrate_const( runge_kutta4< state_type >() , d , x , 0.0 , 1000.0 , 0.1 , write_for_gnuplot() );




    return 0;
}
