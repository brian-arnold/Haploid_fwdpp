#ifndef __SAMPLING_FUNCTIONS_HPP__
#define __SAMPLING_FUNCTIONS_HPP__

#include <utility>
#include <vector>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/insertion_policies.hpp>

using namespace fwdpp ;

/*! \brief Sample the next generation of dipliods in an individual-based
  simulation.  Constant population size case.
  \param r GSL random number generator
  \param gametes Gametes currently in population
  \param diploids Vector of parents from which we sample offspring
  \param mutations Mutations currently in population
  \param mcounts Vector of integers corresponding to counts of each element
  in mutations
  \param N_curr The population size
  \param mu The total mutation rate per gamete
  \param mmodel Mutation model policy
  \param rec_pol Recombination model policy
  \param ff Policy calculating the fitness of a diploid
  \param neutral
  \param selected
  \param f Probability that a mating is a selfing event
  \param mp Policy determining how whether or not to remove fixed variants
  from the gametes.
  \param gpolicy_mut Policy determining how new gametes are added to
  population after a mutation event

  \note diploids will be updated to reflect the new diploid genotypes
  post-sampling (the descedants).  Gametes will be changed by mutation,
  recombination, and sampling.  Mutations will be changed by mutation and
  sampling.
  \return The mean fitness of the parental generation
  \example diploid_ind.cc
  \example pfix.cc
  \example diploid_fixed_sh_ind_lambda.cc
*/
//
// single deme, constant N
//
template <typename gamete_type, typename gamete_cont_type_allocator,
          typename mutation_type, typename mutation_cont_type_allocator,
          typename haploid_fitness_function, typename mutation_model,
          typename recombination_policy,
          template <typename, typename> class gamete_cont_type,
          template <typename, typename> class mutation_cont_type,
          typename mutation_removal_policy = std::true_type>
double sample_haploid(
        // 13 args
        const gsl_rng *r,
        gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
        mutation_cont_type<mutation_type, mutation_cont_type_allocator> &mutations,
        std::vector<uint_t> &mcounts,
        const uint_t &N_curr,
        const double &mu,
        const mutation_model &mmodel,
        const recombination_policy &rec_pol,
        const haploid_fitness_function &ff,
        typename gamete_type::mutation_container &neutral,
        typename gamete_type::mutation_container &selected,
        const double f = 0.,
        const mutation_removal_policy mp = mutation_removal_policy());

/*! \brief Sample the next generation of dipliods in an individual-based
  simulation.  Changing population size case.
  \param r GSL random number generator
  \param gametes Gametes currently in population
  \param diploids Vector of parents from which we sample offspring
  \param mutations Mutations currently in population
  \param mcounts Vector of integers corresponding to counts of each element
  in mutations
  \param N_curr The current population size
  \param N_next The population size after sampling
  \param mu The total mutation rate per gamete
  \param mmodel Mutation model policy
  \param rec_pol Recombination model policy
  \param ff Policy calculating the fitness of a diploid
  \param neutral
  \param selected
  \param f Probability that a mating is a selfing event
  \param mp Policy determining how whether or not to remove fixed variants
  from the gametes.
  \param gpolicy_mut Policy determining how new gametes are added to
  population after a mutation event

  \note diploids will be updated to reflect the new diploid genotypes
  post-sampling (the descedants).  Gametes will be changed by mutation,
  recombination, and sampling.  Mutations will be changed by mutation and
  sampling.
  \return The mean fitness of the parental generation
  \example bneck_selection_ind.cc
*/
//
// single deme, N changing
//
template <typename gamete_type, typename gamete_cont_type_allocator,
          typename mutation_type, typename mutation_cont_type_allocator,
          typename haploid_fitness_function, typename mutation_model,
          typename recombination_policy,
          template <typename, typename> class gamete_cont_type,
          template <typename, typename> class mutation_cont_type,
          typename mutation_removal_policy = std::true_type>
double sample_haploid(
        // 14 args
        const gsl_rng *r,
        gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
        mutation_cont_type<mutation_type, mutation_cont_type_allocator> &mutations,
        std::vector<uint_t> &mcounts,
        const uint_t &N_curr,
        const uint_t &N_next,
        const double &mu,
        const mutation_model &mmodel,
        const recombination_policy &rec_pol,
        const haploid_fitness_function &ff,
        typename gamete_type::mutation_container &neutral,
        typename gamete_type::mutation_container &selected,
        const double f = 0.,
        const mutation_removal_policy mp = mutation_removal_policy());

// this needs to be included after, o.w. compiler doesnt like adding defaults to
//func tmeplate already declared
//#include "sample_haploid_singlepop.cpp"

#endif
