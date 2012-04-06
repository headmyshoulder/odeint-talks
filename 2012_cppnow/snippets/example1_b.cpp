// ...

int main(int argc, char **argv)
{
    state_type x = {{ 10.0 , 1.0 , 1.0 }};
    typedef dense_output_runge_kutta<
        controlled_runge_kutta<
            runge_kutta_dopri5<state_type> > > stepper_type;
    integrate_const(stepper_type(), lorenz, x, 0.0, 10.0, 0.01);
    return 0;
}