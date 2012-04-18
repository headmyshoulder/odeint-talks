/*
 * phase_compacton_lattice.cpp
 *
 *  Created on: Apr 15, 2012
 *      Author: karsten
 */

#include <iostream>
#include <vector>
#include <cmath>

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix.hpp>

namespace odeint = boost::numeric::odeint;

typedef std::vector< double > state_type;

struct phase_compacton_lattice
{
    void operator()( const state_type &x , state_type &dxdt , double t ) const
    {
        size_t N = x.size();
        double tmp = cos( x[0] );
        dxdt[N-1] = tmp;
        dxdt[1] = - tmp;
        dxdt[0] = 0.0;
        for( size_t i=1 ; i<N-2 ; ++i )
        {
            double tmp = cos( x[i] );
            dxdt[i+1] = -tmp;
            dxdt[i-1] += tmp;
        }
        tmp = cos( x[N-2] );
        dxdt[N-1] -= tmp;
        dxdt[N-3] += tmp;
        tmp = cos( x[N-1] );
        dxdt[0] -= tmp;
        dxdt[N-2] += tmp;
    }
};

void initial_condition( state_type &x , size_t start , size_t end , double height )
{
    double hw = double( end - start - 1 );
    for( size_t i = start ; i<end ; ++i )
    {
        x[i] = sin( double( i - start ) / hw * M_PI ) * height;
    }
}

double mod2pi(double x)
{
    return x-double(int(x/2.0/M_PI))*2.0*M_PI+ ((x<0.0)?2.0*M_PI : 0.0);
}

double mod2pi2(double x)
{
    return mod2pi(x+M_PI)-M_PI;
}




int main( int argc , char **argv )
{
    const double dt = 0.01;
    const size_t N = 64;

    odeint::runge_kutta4< state_type > rk4;
    state_type x( N , 0.0 );
    initial_condition( x , 10 , 30 , 2.0 );


    double t = 0.0;
    for( size_t oi=0 ; oi<3000 ; ++oi )
    {
        std::cout << "p [][-pi:pi] '-' w lp pt 7" << "\n";
        for( size_t i=0 ; i<N ; ++i )
            std::cout << i << " " << mod2pi2( x[i] ) << "\n";
        std::cout << "e" << std::endl;

        for( size_t ii=0 ; ii<20 ; ++ii , t+=dt )
            rk4.do_step( phase_compacton_lattice() , x , t , dt );
    }
    return 0;
}
