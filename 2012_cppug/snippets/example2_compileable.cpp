#include <vector>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

typedef std::vector<double> state_type;

struct fpu {
    double m_beta;
    fpu(double beta) : m_beta(beta) { }
    void operator()(const state_type &x, state_type &dxdt, double t) const {
        // ...
    }
};

struct statistics_observer {
    void operator()( const state_type &x , double t ) {
        // write the statistics
    }
};

int main(int argc, char **argv) {
    state_type x(256);
    // initialize x
    integrate_const(runge_kutta4<state_type>(), fpu(1.0), x, 0.0, 10.0, 0.01, statistics_observer());
    return 0;
}