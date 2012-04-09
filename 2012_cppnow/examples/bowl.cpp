/*
 * bowl.cpp
 *
 *  Created on: Apr 6, 2012
 *      Author: karsten
 */


#include <iostream>
#include <array>

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix.hpp>

const size_t num_of_particles = 2;
typedef std::array< double , 4 * num_of_particles > state_type;
// rx0, ry0, rx1, ry1, ..., vx0, vy0, vx1, vy1, ...

class bowl
{
public:

    bowl( void )
    : m_mass( 1.0 ) , m_a( 1.0 ) , m_g( 9.81 ) , m_radius( 0.25 ) ,
      m_interaction_fac( 100.0 ) , m_mu( 0.25 ) { }

    void operator()( const state_type &x , state_type &dxdt , double t ) const
    {
        for( size_t i=0 ; i<num_of_particles ; ++i )
        {
            dxdt[ 2 * i     ] = x[ 2 * num_of_particles + 2 * i     ];
            dxdt[ 2 * i + 1 ] = x[ 2 * num_of_particles + 2 * i + 1 ];
            local_func( dxdt[ 2 * num_of_particles + 2 * i ] , dxdt[ 2 * num_of_particles + 2 * i + 1 ] ,
                    x[ 2 * i ] , x[ 2 * i + 1 ] , x[ 2 * num_of_particles + 2 * i ] , x[ 2 * num_of_particles + 2 * i + 1 ] );
        }

        for( size_t i=1 ; i<num_of_particles ; ++i )
        {
            for( size_t j=0 ; j<i ; ++j )
            {
                interaction_func( dxdt[ 2 * num_of_particles + 2 * i ] , dxdt[ 2 * num_of_particles + 2 * i + 1 ] ,
                        dxdt[ 2 * num_of_particles + 2 * j ] , dxdt[ 2 * num_of_particles + 2 * j + 1 ] ,
                        x[ 2 * i ] , x[ 2 * i + 1 ] , x[ 2 * j ] , x[ 2 * j + 1 ] );
            }
        }
    }

private:

    void local_func( double &dvx , double &dvy , double x , double y , double vx , double vy ) const
    {
        double r = sqrt( x * x + y * y );
        double R = 2.0 * m_a * r;
        dvx = - x / r * m_mass * R / sqrt( 1.0 + R * R ) - m_mu * vx;
        dvy = - y / r * m_mass * R / sqrt( 1.0 + R * R ) - m_mu * vy ;
    }

    void interaction_func( double &dvx1 , double &dvy1 , double &dvx2 , double &dvy2 , double x1 , double y1 , double x2 , double y2 ) const
    {
        double dx = ( x1 - x2 ) , dy = ( y1 - y2 );
        double dist = sqrt( dx * dx + dy * dy );
        if( dist < 2.0 * m_radius )
        {
            double fac = m_interaction_fac * pow( 2.0 * m_radius - dist , 1.5 );
            dvx1 = dx * fac;
            dvx2 = - dx * fac;
            dvy1 = dy * fac;
            dvy2 = - dy * fac;
        }
    }

    double m_mass ;
    double m_a ;
    double m_g ;
    double m_radius;
    double m_interaction_fac;
    double m_mu;
};

struct write_for_gnuplot
{
    void operator()( const state_type &x , double t ) const
    {
        using namespace std;
        cout << "unset key" << endl;
        cout << "sp [-2:2][-2:2][0:10] x*x+y*y,'-' w p pt 7 ps 3 lt 3 " << endl;
        for( size_t i=0 ; i<num_of_particles ; ++i )
        {
            double rx = x[2*i] , ry= x[2*i+1];
            double r2 = rx * rx + ry * ry;
            std::cout << rx << "\t" << ry << "\t" << r2 << endl;
        }
        cout << "e" << endl;
    }
};

struct write1
{
    void operator()( const state_type &x , double t ) const
    {
        std::cout << t;
        for( size_t i=0 ; i<num_of_particles ; ++i )
        {
            double rx = x[2*i] , ry= x[2*i+1];
            double r2 = rx * rx + ry * ry;
            std::cout << "\t" << rx << "\t" << ry ;
            std::cout << "\t" << r2;
        }
        std::cout << "\n";
    }
};

int main( int argc , char **argv )
{
    using namespace boost::numeric::odeint;
    using namespace boost::phoenix::arg_names;
    state_type x = {{ 1.0 , 1.0 , -1.0 , 1.0 , 0.0 , -0.1 , 0.0 , 0.01 }};

//    integrate_const( runge_kutta4< state_type >() , bowl() , x , 0.0 , 100.0 , 0.1 ,
//            std::cout << arg2 << "\t" << arg1[0] << "\t" << arg1[1] << "\t" << arg1[2] << "\t" << arg1[3] << "\n" );

//    integrate_const( runge_kutta4< state_type >() , bowl() , x , 0.0 , 100.0 , 0.1 , write1() );
    integrate_const( runge_kutta4< state_type >() , bowl() , x , 0.0 , 1000.0 , 0.01 , write_for_gnuplot() );




    return 0;
}
