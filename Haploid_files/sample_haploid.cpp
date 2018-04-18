//  -*- C++ -*-

/*!
  \file sample_haploid.tcc

  \brief Definitions of functions for evolving populations of diploids.
*/

#ifndef __FWDPP_SAMPLE_HAPLOID_CPP__
#define __FWDPP_SAMPLE_HAPLOID_CPP__

#include <iostream>
#include <fwdpp/debug.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpp/internal/recycling.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/internal/gamete_cleaner.hpp>
#include <fwdpp/internal/multilocus_rec.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>
#include <misc_functions.hpp>

using namespace fwdpp ;

// single deme, constant N
template <typename gamete_type, typename gamete_cont_type_allocator,
          typename mutation_type, typename mutation_cont_type_allocator,
          typename haploid_fitness_function, typename mutation_model,
          typename recombination_policy,
          template <typename, typename> class gamete_cont_type,
          template <typename, typename> class mutation_cont_type,
          typename mutation_removal_policy>
double
sample_haploid(
    //13 args
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
    const double f,
    const mutation_removal_policy mp)
{
    // run changing N version with N_next == N_curr
    return sample_haploid(r, gametes, mutations, mcounts, N_curr,
                          N_curr, mu, mmodel, rec_pol, ff, neutral,
                          selected, f, mp);
}

// single deme, N changing
template <typename gamete_type, typename gamete_cont_type_allocator,
          typename mutation_type, typename mutation_cont_type_allocator,
          typename haploid_fitness_function, typename mutation_model,
          typename recombination_policy,
          template <typename, typename> class gamete_cont_type,
          template <typename, typename> class mutation_cont_type,
          typename mutation_removal_policy>
double
sample_haploid(
    //14 args
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
const double f,
    const mutation_removal_policy mp)
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

    // test preconditions in debugging mode
    for (const auto &g : gametes)
    {
        assert(gamete_data_sane(g, mutations, mcounts)) ;
    }
    std::cout << "in sample_haploid\n" ;
    assert(mcounts.size() == mutations.size());
    assert(check_sum(gametes, N_curr)); // BA: uses Kevin's check_sum func to cycle through gametes, counting n
    std::cout << "init N: " << get_sum(gametes) << "\n" ;
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
    auto gam_recycling_bin = fwdpp_internal::make_gamete_queue(gametes);


    std::size_t ngametes = gametes.size();
    // vect of fitnesses
    std::vector<double> fitnesses(ngametes);
    double wbar = 0.; // pop'n mean fitness
    // taken from ./fwdpp/io/gamete.hpp:
    for (std::size_t i = 0; i < ngametes; ++i)
        {
            std::cout << "i:" <<i << "N:" << get_sum(gametes) << "\n" ;
          //from fitness_models_haploid.cpp, ff returns double
            fitnesses[i] = ff(gametes[i], mutations);
            std::cout << fitnesses[i] << "\t" ;
            // (gamete fitness)*(pop freq) = total sampling probability
            fitnesses[i] = fitnesses[i]*static_cast<double>(gametes[i].n) ;
            std::cout << fitnesses[i] << "\n" ;
            gametes[i].n = 0;
            wbar += fitnesses[i];
        }
    // divide by popsize, which may != ngametes since gametes may have n>1
    wbar /= double(N_curr);
#ifndef NDEBUG
    for (const auto &g : gametes)
        assert(!g.n);
#endif

    /*
      This is a lookup table for rapid sampling of diploids proportional to
      their fitnesses.
      This is a unique_ptr wrapper around an object from the GNU Scientific
      Library.  A custom deleter
      is required to make this work, which is why there is no cleanup call
      down below.
    */
    fwdpp_internal::gsl_ran_discrete_t_ptr lookup(
        gsl_ran_discrete_preproc(ngametes, fitnesses.data()));
    // Fill in the next generation!
    std::cout << "ngametes: " << ngametes << "\n" ;
    for (std::size_t i = 0; i < N_next; ++i)
        {
            // Choose parent 1 based on fitness
            auto p1 = gsl_ran_discrete(r, lookup.get());
            std::cout << p1 << "\n" ;
            gametes[p1].n++ ;
            //For Recombination
            /*
            auto p2 = p1 ;
            // select different individual p2 as DNA donor
            while(p2 == p1){
                p2 = gsl_ran_discrete(r, lookup.get());
            }
            assert(p1 < ngametes);
            assert(p2 < ngametes);
            */

            /*
              These are the gametes from each parent.
              This is a trivial assignment if keys.
            */
            //auto p1g1 = parents[p1].first;
            //auto p1g2 = parents[p1].second;
            //auto p2g1 = parents[p2].first;
            //auto p2g2 = parents[p2].second;

    ////////
    // make overloaded versions of functions in fwdpp/mutate_recombine.hpp?
    // first passes can have no rec... but have to modify this func in order to mutate!
           /*
            mutate_recombine_update(
                r, gametes, mutations,
                std::make_tuple(p1g1, p1g2, p2g1, p2g2), rec_pol, mmodel,
                mu, gam_recycling_bin, mut_recycling_bin, dip, neutral,
                selected);
            */
    ////////


    }
    std::cout << "N aft: " << get_sum(gametes) << "\n" ;
    assert(check_sum(gametes, N_next));
#ifndef NDEBUG
    /*
for (const auto &dip : diploids)
        {
            assert(gametes[dip.first].n > 0);
            assert(gametes[dip.first].n <= 2 * N_next);
            assert(gametes[dip.second].n > 0);
            assert(gametes[dip.second].n <= 2 * N_next);
        }
*/    
#endif
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

      Further, the function is hard to optimize. Recall that gametes store
      mutations in order
      according to position.  Thus, when we go from a gamete to a position
      in mcounts, we are
      accessing the latter container out of order with respect to location
      in memory.  process_gametes
      is thus the "scatter" part of a "scatter-gather" idiom.  Modern x86
      CPU have little available
      for vectorizing such cases.  I've experimented with CPU intrinsics to
      attempt memory prefetches,
      but never saw any performance improvement, and the code got complex,
      and possibly less portable.

      The implementation is in fwdpp/internal/sample_diploid_helpers.hpp
     */
    fwdpp_internal::process_gametes(gametes, mutations, mcounts);
    assert(mcounts.size() == mutations.size());
#ifndef NDEBUG
    for (const auto &mc : mcounts)
        {
            assert(mc <= 2 * N_next);
        }
#endif
    
    for (const auto &g : gametes)
    {
        assert(gamete_data_sane(g, mutations, mcounts)) ;
    }

    /*
      The last thing to do is handle fixations.  In many contexts, we
      neither want nor need
      to keep indexes to fixed variants in our gametes.  Such decisions are
      implemented via
      simple policies, which are in the variable 'mp'.

      The implementation is in fwdpp/internal/gamete_cleaner.hpp.

      The implementation is the "erase/remove idiom" (Effective STL, Item
      32), but with a twist
      that the function will exit early if there are no fixations present
      in the population at
      the moment.

      Example policies are fwdpp::remove_nothing and fwdpp::remove_neutral,
      both found
      in fwdpp/fwd_functional.hpp.  If mp is std::true_type, then all
      fixations (e.g., neutral
      and selected)  will be removed from all gametes.
    */
    fwdpp_internal::gamete_cleaner(gametes, mutations, mcounts, N_next,
                                   mp);
    return wbar;
}



#endif
