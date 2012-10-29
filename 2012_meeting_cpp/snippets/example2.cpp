typedef std::vector<double> state_type;

struct fpu {
    double m_beta;
    fpu(double beta) : m_beta(beta) { }
    void operator()(const state_type &q, state_type &dpdt) const {
        // ...
    }
};

void statistics_observer ( const state_type &x , double t ) {
    // write the statistics
}

int main(int argc, char **argv) {
    state_type q(256),p(256);
    // initialize q,p
    integrate_const(symplectic_rkn_sb3a_mclachlan<state_type>(), fpu(1.0), make_pair(q,p), 0.0, 10.0, 0.01, statistics_observer());
    return 0;
}