#ifndef SA_SOLVER_HPP
#define SA_SOLVER_HPP

#include <random>
#include "wdb.hpp"
typedef wdb::entities::ising::model hamiltonian_type;
typedef wdb::entities::ising::property property_type;
// FIXME: duplicate symbols when using .cpp and linking

#include "result.hpp"

class sa_solver
// Simulated annealing algorithm to find ground state of spin glass
// Input: Hamiltonian of the spin glass (when constructing a class object)
// Ouput: Energy and Spin configuration (accessible via get_config())
//
// run(...) initialized and runs SA
{
public:

  //initialise sa solver with hamiltonian
  sa_solver(const hamiltonian_type& H)
  : N_(H.num_nodes()),
    H_(H)
  {
  }


  //single run of sa from random initial state on hamiltonian H_
  result run(const double beta0,
             const double beta1,
             const std::size_t Ns,
             const std::size_t seed)
  {
    // instantiate a generator for random numbers
    std::minstd_rand0 lcg;
    lcg.seed(seed);

    // binary and real distributions
    std::uniform_int_distribution<int> binaries(0, 1);
    std::uniform_real_distribution<double> realnums(0.0, 1.0);

    // stop STL from reallocating space as vector grows
    spins_.reserve(N_);

    // generate N random binary states
    for(unsigned i = 0; i < N_; ++i)
      spins_.push_back(binaries(lcg));

    // compute initial energy
    double E = H_.total_energy(spins_);

    // perform annealing
    for(unsigned s = 0; s < Ns; ++s){
      const double beta(beta0 + (beta1-beta0)/(Ns-1)*s);
      for(unsigned i = 0; i < N_; ++i){
        const double dE(H_.delta_energy(spins_,i));
        if(dE <= 0.0 || realnums(lcg) < std::exp(-beta*dE)){
          spins_[i] = spins_[i] ^ 1;
          E += dE;
        }
      }
    }

    // return result
    result res;
    res.E_=E;
    res.spins_.swap(spins_);
    return res;

  }

private:

  const std::size_t N_;
  const hamiltonian_type& H_;

  std::vector<int> spins_;
};

// SA solver
// Input: H is the Hamiltonian of the spin glass
//        beta0 is the starting inverse temperature (use 0.1)
//        beta1 is the ending inverse temperature (use 3.0)
//        Ns is the number of Monte Carlo sweeps per run of SA
//        seed initialized the random number generator, use a different one for each repetition
// Output: Configuration of spins and energy
result solve(const hamiltonian_type& H,
             const double beta0,
             const double beta1,
             const std::size_t Ns,
             const std::size_t seed)
{
  sa_solver solver(H);
  return solver.run(beta0, beta1, Ns, seed);
}

#endif
