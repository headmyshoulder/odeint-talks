/*
 *  van_der_pol_stiff.cpp
 *  
 *  Created on: Dec 12, 2011
 *      Author: rajeev
 */

#include <iostream>
#include <fstream>
#include <utility>
#include <array>

#include <boost/numeric/odeint.hpp>

#include <boost/spirit/include/phoenix.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef std::array< double , 2 > state_type;

struct vdp
{
    void operator()( const state_type &x , state_type &dxdt , double t )
    {
    	const double mu = 1.0;
    	dxdt[0] = x[1];
    	dxdt[1] = -x[0] - mu * x[1] * (x[0]*x[0]-1.0);
    }
};

struct writer
{
	std::ostream &m_out;
	writer( std::ostream &out ) : m_out( out ) { }

	template< class State , class Time >
	void operator()( const State &x , Time t ) const
	{
		m_out << t << " " << x[0] << " " << x[1] << "\n";
	}
};


int main( int argc , char **argv )
{
	using namespace boost::phoenix::arg_names;

	state_type x = {{ 1.0 , 1.0 }};

	ofstream fout1( "vdp_constant.dat" );
	ofstream fout2( "vdp_adaptive.dat" );

	integrate_const( make_dense_output( 1.0e-3 , 1.0e-3 , runge_kutta_dopri5< state_type >() ) ,
	            vdp() , x , 0.0 , 102.0 , 0.1 );

	state_type y = x;
    integrate_const( make_dense_output( 1.0e-4 , 1.0e-4 , runge_kutta_dopri5< state_type >() ) ,
            vdp() , x , 0.0 , 10.0 , 0.1 , writer( fout1 ) );


    x = {{ 1.0 , 1.0 }};

    integrate_adaptive( make_dense_output( 1.0e-4 , 1.0e-4 , runge_kutta_dopri5< state_type >() ) ,
                vdp() , y , 0.0 , 10.0 , 0.1 , writer( fout2 ) );

    return 0;
}
