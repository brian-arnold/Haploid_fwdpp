//  -*- C++ -*-

/*!
  \file sample_haploid.tcc

  \brief Definitions of functions for evolving populations of diploids.
*/

#ifndef __SAMPLE_HAPLOID_STRUCTURE_CPP__
#define __SAMPLE_HAPLOID_STRUCTURE_CPP__

#include <iostream>
#include <fwdpp/debug.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpp/internal/recycling.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/internal/gamete_cleaner.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>
#include "misc_functions.hpp"
#include "mutate_recombine_haploid.hpp"
#include "sample_haploid_structure.hpp"

using namespace fwdpp ;

//
// For structured deme, constant N, includes haploid vector
//
template <typename gamete_type, typename gamete_cont_type_allocator,
            typename mutation_type, typename mutation_cont_type_allocator,
            typename haploid_geno_t, typename haploid_vector_type_allocator,
            typename haploid_fitness_function, typename mutation_model,
            typename recombination_policy,
            template <typename, typename> class gamete_cont_type,
            template <typename, typename> class mutation_cont_type,
            template <typename, typename> class haploid_vector_type,
            typename mutation_removal_policy = std::true_type>
double
sample_haploid_structure(
    //14 args
    const gsl_rng *r,
    gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
    haploid_vector_type<haploid_geno_t, haploid_vector_type_allocator> &haploids,
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
    const mutation_removal_policy mp = mutation_removal_policy())
{
    // run changing N version with N_next == N_curr
    return sample_haploid_structure(r, gametes, haploids, mutations, mcounts, N_curr,
                          N_curr, mu, mmodel, rec_pol, ff, neutral,
                          selected, f, mp);
}

//
// For structured deme, N changing, includes haploid vector
//
template <typename gamete_type, typename gamete_cont_type_allocator,
            typename mutation_type, typename mutation_cont_type_allocator,
            typename haploid_geno_t, typename haploid_vector_type_allocator,
            typename haploid_fitness_function, typename mutation_model,
            typename recombination_policy,
            template <typename, typename> class gamete_cont_type,
            template <typename, typename> class mutation_cont_type,
            template <typename, typename> class haploid_vector_type,
            typename mutation_removal_policy = std::true_type>
double
sample_haploid_structure(
    //15 args
   const gsl_rng *r,
   gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
   haploid_vector_type<haploid_geno_t, haploid_vector_type_allocator> &haploids,
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
   const mutation_removal_policy mp = mutation_removal_policy())
{
    /*
      The main part of fwdpp does not throw exceptions.
      Rather, testing is performed via C's assert macro.
      This macro should be disabled in "production" builds via
      -DNEBUG as is standard practice.  It is the developer's
      responsibility to properly set up a build system to distinguish
      'debug' from 'production' builds.

      More complex debugging blocks will be wrapped in #ifndef
      NDEBUG/#endif
      blocks as needed.

      Compiling in a 'debug' mode slows simulations down several-fold.
    */

    // DEBUG
    assert(happopdata_sane(haploids, gametes, mutations, mcounts));
    assert(mcounts.size() == mutations.size());
    assert(check_sum(gametes, N_curr)); // BA: uses Kevin's check_sum func to cycle through gametes, counting n
    assert(mcounts.size() == mutations.size());

    /*
      The mutation and gamete containers contain both extinct and extant
      objects.
      The former are useful b/c the represent already-allocated memory.
      The library
      uses these extinct objects to 'recycle' them into new objects.  The
      function calls
      below create FIFO queues of where extinct objects are.  These queues
      are passed to
      mutation and recombination functions and used to decide if recyling
      is possible or
      if a new object needs to be 'emplace-back'-ed into a container.

      The type of the FIFO queue is abstracted with the name
      fwdpp::fwdpp_internal::recycling_bin_t,
      which is a C++11 template alias.

      The details of recycling are implemented in
      fwdpp/internal/recycling.hpp
    */
    auto mut_recycling_bin = fwdpp_internal::make_mut_queue(mcounts);
    // this tags gametes with n=0 for recycling; b4 they're marked as 0 during fitness calc
    auto gam_recycling_bin = fwdpp_internal::make_gamete_queue(gametes);
   
    // vect of fitnesses
    std::vector<double> fitnesses(haploids.size());
    double wbar = 0.; // pop'n mean fitness
    for (uint_t i = 0; i < N_curr; ++i)
        {
            //from fitness_models_haploid.cpp, ff returns double
            fitnesses[i] = ff(gametes[haploids[i]], mutations);
            gametes[haploids[i]].n = 0;
            wbar += fitnesses[i];
        }
    // divide by popsize, which may != ngametes since gametes may have n>1
    wbar /= double(N_curr);

    // DEBUG
    for (const auto &g : gametes)
        assert(!g.n);

    
    /*
      This is a lookup table for rapid sampling of diploids proportional to
      their fitnesses.
      This is a unique_ptr wrapper around an object from the GNU Scientific
      Library.  A custom deleter
      is required to make this work, which is why there is no cleanup call
      down below.
    */
    fwdpp_internal::gsl_ran_discrete_t_ptr lookup_sel(
        gsl_ran_discrete_preproc(N_curr, fitnesses.data()));
    const auto parents(haploids); // Copy the parents, which is trivally
    
    // Change the population size
    if (haploids.size() != N_next)
    {
        haploids.resize(N_next);
    }
    assert(haploids.size() == N_next);
    // Fill in the next generation!
    for (auto &hap : haploids)
        {
            // Choose parent 1 based on fitness
            auto p1 = gsl_ran_discrete(r, lookup_sel.get());
            // Donor for recombination
            // RANDOMLY sample for recombinant
            auto p2 = gsl_rng_uniform_int(r, N_curr) ;
            
            assert(p1 < parents.size());
            assert(p2 < parents.size());

            auto h1 = parents[p1] ; // need to cp from parents,
            auto h2 = parents[p2] ; // haploids being modified

            mutate_recombine_update_haploid(
                r, gametes, mutations,
                std::make_tuple(h1, h2), rec_pol, mmodel,
                mu, gam_recycling_bin, mut_recycling_bin, hap,
                neutral, selected);
    }

    // DEBUG
    assert(check_sum(gametes, N_next));
    for (const auto &hap : haploids)
        {
            assert(gametes[hap].n > 0);
            assert(gametes[hap].n <= N_next);
        }
     
    /*
      At the end of the above loop, we have a bunch of new diploids
      that are all recombined and mutated sampling of the parental
      generation.

      Our problem is that we no longer know how many times each mutation is
      present, which
      is corrected by the following call.

      Although the implementation of process_gametes is super-trivial, it
      is actually the
      most computationally-expensive part of a simulation once mutation
      rates are large.

      The implementation is in fwdpp/internal/sample_diploid_helpers.hpp
     */
    fwdpp_internal::process_gametes(gametes, mutations, mcounts);
    assert(mcounts.size() == mutations.size());

    // DEBUG
    for (const auto &mc : mcounts)
        {
            assert(mc <= N_next);
        }
    assert(happopdata_sane(haploids, gametes, mutations, mcounts));

    
    /*
      The last thing to do is handle fixations.  In many contexts, we
      neither want nor need
      to keep indexes to fixed variants in our gametes.  Such decisions are
      implemented via
      simple policies, which are in the variable 'mp'.

      The implementation is in fwdpp/internal/gamete_cleaner.hpp.
    */
    fwdpp_internal::gamete_cleaner(gametes, mutations, mcounts, N_next,
                                   mp);
    return wbar;
}



#endif
