#include "sa_solver.hpp"

int main(int argc, char* argv[]) {
// Runs SA with many inputfiles
// Input: beta0 is the starting temperature of SA (use 0.1 for bimodal instances)
//        beta1 is the ending temperature of SA (use 3.0 for bimodal instances)
//        Ns is the number of Monte Carlo Sweeps within each run of SA
//        num_rep is the number of repetitions of SA (with different seeds)

  hamiltonian_type& H = *(hamiltonian_type*)argv[1];
  property_type& P = *(property_type*)argv[2];

  const unsigned Ns = P.get_param<unsigned>("nsweeps");
  const double beta0 = P.get_param<double>("b0");
  const double beta1 = P.get_param<double>("b1");
  const unsigned seed = P.get_seed();

  auto res = solve(H,beta0,beta1,Ns,seed);
  P.set_cfg(res.spins_,res.E_);

  return 0;
}
