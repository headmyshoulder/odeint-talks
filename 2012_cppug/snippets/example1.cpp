#include <boost/numeric/odeint.hpp>
#include <tr1/array>

using namespace boost::numeric::odeint;

typedef std::tr1::array<double,3> state_type;

void lorenz(const state_type &x, state_type &dxdt, double t) {
    // ...
}

int main(int argc, char **argv) {
    state_type x = {{10.0, 10.0, 10.0}};
    typedef dense_output_runge_kutta<
        controlled_runge_kutta<
            runge_kutta_dopri5<state_type> > > stepper_type;
    integrate_const(stepper_type(), lorenz, x, 0.0, 10.0, 0.01);
    return 0;
}