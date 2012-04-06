/*
 * harmonic_oscillator.cpp
 *
 *  Created on: Apr 1, 2012
 *      Author: karsten
 */

#include <iostream>
#include <array>

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix.hpp>

typedef std::array< double , 2 > state_type;

class harmonic_oscillator
{
public:

    harmonic_oscillator( double omega ) : m_omega( omega ) { }

    void operator()( const state_type &x , state_type &dxdt , double t ) const
    {
        dxdt[0] = x[1];
        dxdt[1] = - m_omega * x[0];
    }

private:

    double m_omega;
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
    state_type x = {{ 0.1 , 0.0 }};
//    integrate_const( runge_kutta4< state_type >() , harmonic_oscillator( 1.0 ) , x , 0.0 , 10.0 , 0.01 ,
//            std::cout << arg2 << "\t" << arg1[0] << "\t" << arg1[1] << "\n" );

    integrate_const( runge_kutta4< state_type >() , harmonic_oscillator( 1.0 ) , x , 0.0 , 100.0 , 0.1 , write_for_gnuplot() );




    return 0;
}
