/*
 * pendulum1.cpp
 *
 *  Created on: Apr 15, 2012
 *      Author: karsten
 */

#include <iostream>
#include <array>

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix.hpp>

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
    pendulum p( 0.0 , 0.0 , 0.0 );

    state_type x = {{ 1.0 , 0.0 }};
    double t = 0.0;
    const double dt = 0.025;

    std::cout << "unset key" << "\n";
    std::cout << "set size square" << std::endl;

    for( size_t i=0 ; i<1000 ; ++i )
    {
        std::cout << "p [-1.1:1.1][-1.1:1.1]  '-' w lp pt 7 ps 2 lw 3" << std::endl;
        std::cout << "0.0 0.0" << "\n";
        std::cout << -sin(x[0]) << " " << -cos(x[0]) << "\n";
        std::cout << "e" << std::endl;

        rk4.do_step( p , x , t , dt );
        t += dt;
    }

    {
        odeint::runge_kutta_fehlberg78< state_type > stepper;
    }

    {
        odeint::runge_kutta_dopri5< state_type > stepper;
    }

    {
        double dt = 0.025;
        auto stepper = make_controlled( 1.0e-6 , 1.0e6 ,  odeint::runge_kutta_fehlberg78< state_type >() );
        odeint::controlled_step_result res = stepper.try_step( p , x , t , dt );

        double t_end = 1000.0;
        while( t < t_end )
        {
            odeint::controlled_step_result res = stepper.try_step( p , x , t , dt );
            while( res != odeint::success )
            {
                res = stepper.try_step( p , x , t , dt );
            }
        }

        double t_start = 0.0;
        integrate_adaptive( stepper , p , x , t_start , t_end , dt );

        namespace phoenix = boost::phoenix;
        using namespace phoenix::arg_names;
        integrate_adaptive( stepper , p , x , t_start , t_end , dt ,
                std::cout << arg2 << "\t" << arg1[0] << "\t" << arg1[1] << "\n" );
    }




//    {
//        odeint::runge_kutta4< state_type > rk4;
//        pendulum p( 0.1 , 1.05 , 1.5 );
//
//        state_type x = {{ 1.0 , 0.0 }};
//        double t = 0.0;
//        const double dt = 0.01;
//
//        rk4.do_step( p , x , t , dt );
//        t += dt;
//    }


//    odeint::runge_kutta4< state_type > rk4;
//    pendulum p( 0.1 , 1.05 , 1.5 );
//
//    state_type x = {{ 1.0 , 0.0 }};
//    double t = 0.0;
//    const double dt = 0.01;
//
//    std::cout << t << " " << x[0] << " " << x[1] << "\n";
//    for( size_t i=0 ; i<10 ; ++i )
//    {
//      rk4.do_step( p , x , t , dt );
//      t += dt;
//      std::cout << t << " " << x[0] << " " << x[1] << "\n";
//    }

    return 0;
}
