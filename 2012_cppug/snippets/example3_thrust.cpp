typedef thrust::device_vector<double> state_type;
typedef runge_kutta4<state_type ,value_type ,state_type ,value_type ,thrust_algebra ,thrust_operations > stepper_type;

struct lorenz_system {
    
    lorenz_system(size_t N ,const state_type &beta)
    : m_N(N) , m_beta(beta) {}

    void operator()(  const state_type &x , state_type &dxdt , double t ){
        // ..
    }

    size_t m_N;
    const state_type &m_beta;
};

int main( int arc , char* argv[] )
{
    const size_t N = 1024;

    vector<value_type> beta_host(N);
    for( size_t i=0 ; i<N ; ++i )
        beta_host[i] = 56.0 + value_type( i ) * ( 56.0 ) / value_type( N - 1 );
    state_type beta = beta_host;
    state_type x( 3 * N , 10.0 );
    integrate_const( stepper_type() , lorenz(N,beta) , x , 0.0 , 10.0 , 0.01 );

    return 0;
}
