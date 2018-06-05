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
// For structured deme, whole individuals migrate, includes haploid vector
//
template <typename structpoptype,
            typename haploid_fitness_function, typename mutation_model,
            typename recombination_policy,
            typename mutation_removal_policy = std::true_type>
void
sample_haploid_struct_indmig(
    //12 args
   const gsl_rng *r,
   structpoptype &pop,
   //gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
   //haploid_vector_type<haploid_geno_t, haploid_vector_type_allocator> &haploids,
   //mutation_cont_type<mutation_type, mutation_cont_type_allocator> &mutations,
   //std::vector<uint_t> &mcounts,
   const uint_t &N1,
   const uint_t &N2,
   const double m12,
   const double m21,
   const double &mu,
   const mutation_model &mmodel,
   const recombination_policy &rec_pol,
   const haploid_fitness_function &ff,
   //typename gamete_type::mutation_container &neutral,
   //typename gamete_type::mutation_container &selected,
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
    /*
    assert(happopdata_sane(pop.haploids, pop.gametes, pop.mutations, pop.mcounts));
    assert(pop.mcounts.size() == pop.mutations.size());
    assert(check_sum(pop.gametes, (N1 + N2)));
    */
     
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
    auto mut_recycling_bin = fwdpp_internal::make_mut_queue(pop.mcounts);
    // this tags gametes with n=0 for recycling; b4 they're marked as 0 during fitness calc
    auto gam_recycling_bin = fwdpp_internal::make_gamete_queue(pop.gametes);
   
    // Migration and build pop-specific lookup tables:
    auto lookups = migrate_and_calc_fitness(r, pop, ff, N1, N2, m12, m21);
    
    // DEBUG
    /*
    for (const auto &g : pop.gametes)
        assert(!g.n);
     */

    const auto parents(pop.haploids); // Copy the parents

    // Fill in the next generation!
    // We generate the offspring for deme 1 first, and then for deme 2.
    // because of migration, the number of individuals for N1(N2) increases
    // in lookup table, but we only sample N1(N2) ind's from it,
    // conserving pop sizes.
    // Lookup tables from migrate_and_calc_fitness allow individuals
    // from N2 to get included while sampling for N1
    for (uint_t i = 0; i < N1 + N2; ++i)
    {
        std::size_t p1 = std::numeric_limits<std::size_t>::max();
        std::size_t p2 = std::numeric_limits<std::size_t>::max();
        if (i < N1) // pick parents from pop 1
        {
            // some of these may be indices from N2 if mig occurred
            p1 = lookups.parents1[gsl_ran_discrete(r, lookups.seln1.get())];
            p2 = lookups.parents1[gsl_ran_discrete(r, lookups.rand1.get())];
        }
        else // pick parents from pop 2
        {
            p1 = lookups.parents2[gsl_ran_discrete(r, lookups.seln2.get())];
            p2 = lookups.parents2[gsl_ran_discrete(r, lookups.rand2.get())];
        }
        assert(p1 < parents.size());
        assert(p2 < parents.size());
        
        auto h1 = parents[p1] ; // need to cp from parents,
        auto h2 = parents[p2] ; // haploids being modified
        
        // mutant/recombinant made from parent haplotypes h1+h2, assigned to
        // pop.haploids[i], transferring haplotypes associated with indices
        // from N2 to indices associated with N1
        mutate_recombine_update_haploid(r, pop.gametes, pop.mutations,
                                        std::make_tuple(h1, h2), rec_pol, mmodel,
                                        mu, gam_recycling_bin, mut_recycling_bin,
                                        pop.haploids[i], pop.neutral, pop.selected);
        
    }
    
    // DEBUG
    /*
    assert(check_sum(pop.gametes, (N1 + N2)));
    for (const auto &hap : pop.haploids)
        {
            assert(pop.gametes[hap].n > 0);
            assert(pop.gametes[hap].n <= (N1 + N2));
        }
     */
     
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
    fwdpp_internal::process_gametes(pop.gametes, pop.mutations, pop.mcounts);
    assert(pop.mcounts.size() == pop.mutations.size());

    // DEBUG
    /*
    for (const auto &mc : pop.mcounts)
        {
            assert(mc <= (N1 + N2));
        }
    assert(happopdata_sane(pop.haploids, pop.gametes, pop.mutations, pop.mcounts));
     */
    
    /*
      The last thing to do is handle fixations.  In many contexts, we
      neither want nor need
      to keep indexes to fixed variants in our gametes.  Such decisions are
      implemented via
      simple policies, which are in the variable 'mp'.

      The implementation is in fwdpp/internal/gamete_cleaner.hpp.
    */
    fwdpp_internal::gamete_cleaner(pop.gametes, pop.mutations, pop.mcounts, (N1 + N2),
                                   mp);
}

//
// For structured deme, whole individuals migrate, includes haploid vector
//
template <typename structpoptype,
typename haploid_fitness_function, typename mutation_model,
typename recombination_policy,
typename mutation_removal_policy = std::true_type>
void
sample_haploid_struct_disrupted_indmig(
                             //14 args
                             const gsl_rng *r,
                             structpoptype &pop,
                             //gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
                             //haploid_vector_type<haploid_geno_t, haploid_vector_type_allocator> &haploids,
                             //mutation_cont_type<mutation_type, mutation_cont_type_allocator> &mutations,
                             //std::vector<uint_t> &mcounts,
                             const uint_t &N1,
                             const uint_t &N2,
                             double m12,
                             double m21,
                             const double &mu,
                             const mutation_model &mmodel,
                             const recombination_policy &rec_pol,
                             const unsigned &gen,
                             const unsigned &ngenMig,
                             const haploid_fitness_function &ff,
                             //typename gamete_type::mutation_container &neutral,
                             //typename gamete_type::mutation_container &selected,
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
    /*
     assert(happopdata_sane(pop.haploids, pop.gametes, pop.mutations, pop.mcounts));
     assert(pop.mcounts.size() == pop.mutations.size());
     assert(check_sum(pop.gametes, (N1 + N2)));
     */
    
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
    auto mut_recycling_bin = fwdpp_internal::make_mut_queue(pop.mcounts);
    // this tags gametes with n=0 for recycling; b4 they're marked as 0 during fitness calc
    auto gam_recycling_bin = fwdpp_internal::make_gamete_queue(pop.gametes);
    
    
    // Disrupt migration after ngenMig generations
    if( gen > ngenMig ){
        m12 = 0.0 ;
        m21 = 0.0 ;
    }
    auto lookups = migrate_and_calc_fitness(r, pop, ff, N1, N2, m12, m21);
    
    // DEBUG
    /*
     for (const auto &g : pop.gametes)
     assert(!g.n);
     */
    
    const auto parents(pop.haploids); // Copy the parents
    
    // Fill in the next generation!
    // We generate the offspring for deme 1 first, and then for deme 2.
    // because of migration, the number of individuals for N1(N2) increases
    // in lookup table, but we only sample N1(N2) ind's from it,
    // conserving pop sizes.
    // Lookup tables from migrate_and_calc_fitness allow individuals
    // from N2 to get included while sampling for N1
    for (uint_t i = 0; i < N1 + N2; ++i)
    {
        std::size_t p1 = std::numeric_limits<std::size_t>::max();
        std::size_t p2 = std::numeric_limits<std::size_t>::max();
        if (i < N1) // pick parents from pop 1
        {
            // some of these may be indices from N2 if mig occurred
            p1 = lookups.parents1[gsl_ran_discrete(r, lookups.seln1.get())];
            p2 = lookups.parents1[gsl_ran_discrete(r, lookups.rand1.get())];
        }
        else // pick parents from pop 2
        {
            p1 = lookups.parents2[gsl_ran_discrete(r, lookups.seln2.get())];
            p2 = lookups.parents2[gsl_ran_discrete(r, lookups.rand2.get())];
        }
        assert(p1 < parents.size());
        assert(p2 < parents.size());
        
        auto h1 = parents[p1] ; // need to cp from parents,
        auto h2 = parents[p2] ; // haploids being modified
        
        // mutant/recombinant made from parent haplotypes h1+h2, assigned to
        // pop.haploids[i], transferring haplotypes associated with indices
        // from N2 to indices associated with N1
        mutate_recombine_update_haploid(r, pop.gametes, pop.mutations,
                                        std::make_tuple(h1, h2), rec_pol, mmodel,
                                        mu, gam_recycling_bin, mut_recycling_bin,
                                        pop.haploids[i], pop.neutral, pop.selected);
        
    }
    
    // DEBUG
    /*
     assert(check_sum(pop.gametes, (N1 + N2)));
     for (const auto &hap : pop.haploids)
     {
     assert(pop.gametes[hap].n > 0);
     assert(pop.gametes[hap].n <= (N1 + N2));
     }
     */
    
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
    fwdpp_internal::process_gametes(pop.gametes, pop.mutations, pop.mcounts);
    assert(pop.mcounts.size() == pop.mutations.size());
    
    // DEBUG
    /*
     for (const auto &mc : pop.mcounts)
     {
     assert(mc <= (N1 + N2));
     }
     assert(happopdata_sane(pop.haploids, pop.gametes, pop.mutations, pop.mcounts));
     */
    
    /*
     The last thing to do is handle fixations.  In many contexts, we
     neither want nor need
     to keep indexes to fixed variants in our gametes.  Such decisions are
     implemented via
     simple policies, which are in the variable 'mp'.
     
     The implementation is in fwdpp/internal/gamete_cleaner.hpp.
     */
    fwdpp_internal::gamete_cleaner(pop.gametes, pop.mutations, pop.mcounts, (N1 + N2),
                                   mp);
}







//
// For structured deme, migration only from interpopulation recombination
//
template <typename structpoptype,
typename haploid_fitness_function, typename mutation_model,
typename recombination_policy,
typename mutation_removal_policy = std::true_type>
void
sample_haploid_struct_recmig(
                             //12 args
                             const gsl_rng *r,
                             structpoptype &pop,
                             //gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
                             //haploid_vector_type<haploid_geno_t, haploid_vector_type_allocator> &haploids,
                             //mutation_cont_type<mutation_type, mutation_cont_type_allocator> &mutations,
                             //std::vector<uint_t> &mcounts,
                             const uint_t &N1,
                             const uint_t &N2,
                             const double m12,
                             const double m21,
                             const double &mu,
                             const mutation_model &mmodel,
                             const recombination_policy &rec_pol,
                             const haploid_fitness_function &ff,
                             //typename gamete_type::mutation_container &neutral,
                             //typename gamete_type::mutation_container &selected,
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
    /*
     assert(happopdata_sane(pop.haploids, pop.gametes, pop.mutations, pop.mcounts));
     assert(pop.mcounts.size() == pop.mutations.size());
     assert(check_sum(pop.gametes, (N1 + N2)));
     */
    
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
    auto mut_recycling_bin = fwdpp_internal::make_mut_queue(pop.mcounts);
    // this tags gametes with n=0 for recycling; b4 they're marked as 0 during fitness calc
    auto gam_recycling_bin = fwdpp_internal::make_gamete_queue(pop.gametes);
    
    // Migration and build pop-specific lookup tables:
    auto lookups = calc_fitness_twopop(r, pop, ff, N1, N2);
    
    // DEBUG
    /*
     for (const auto &g : pop.gametes)
     assert(!g.n);
     */
    
    const auto parents(pop.haploids); // Copy the parents
    
    // Fill in the next generation!
    // We generate the offspring for deme 1 first, and then for deme 2.
    // because of migration, the number of individuals for N1(N2) increases
    // in lookup table, but we only sample N1(N2) ind's from it,
    // conserving pop sizes.
    // Lookup tables from migrate_and_calc_fitness allow individuals
    // from N2 to get included while sampling for N1
    for (uint_t i = 0; i < N1 + N2; ++i)
    {
        std::size_t p1 = std::numeric_limits<std::size_t>::max();
        std::size_t p2 = std::numeric_limits<std::size_t>::max();
        if (i < N1) // pick parents from pop 1
        {
            // some of these may be indices from N2 if mig occurred
            p1 = lookups.parents1[gsl_ran_discrete(r, lookups.seln1.get())];
            if(gsl_rng_uniform(r) < m21){
                // randomly pick rec. donor from N2 pop
                p2 = lookups.parents2[gsl_ran_discrete(r, lookups.rand2.get())];
            }else{
                // randomly pick rec. donor from same (N1) pop
                p2 = lookups.parents1[gsl_ran_discrete(r, lookups.rand1.get())];
            }
        }
        else // pick parents from pop 2
        {
            p1 = lookups.parents2[gsl_ran_discrete(r, lookups.seln2.get())];
            if(gsl_rng_uniform(r) < m12){
                p2 = lookups.parents1[gsl_ran_discrete(r, lookups.rand1.get())];
            }else{
                p2 = lookups.parents2[gsl_ran_discrete(r, lookups.rand2.get())];
            }
        }
        assert(p1 < parents.size());
        assert(p2 < parents.size());
        
        auto h1 = parents[p1] ; // need to cp from parents,
        auto h2 = parents[p2] ; // haploids being modified
        
        // mutant/recombinant made from parent haplotypes h1+h2, assigned to
        // pop.haploids[i], transferring haplotypes associated with indices
        // from N2 to indices associated with N1
        mutate_recombine_update_haploid(r, pop.gametes, pop.mutations,
                                        std::make_tuple(h1, h2), rec_pol, mmodel,
                                        mu, gam_recycling_bin, mut_recycling_bin,
                                        pop.haploids[i], pop.neutral, pop.selected);
        
    }
    
    // DEBUG
    /*
     assert(check_sum(pop.gametes, (N1 + N2)));
     for (const auto &hap : pop.haploids)
     {
     assert(pop.gametes[hap].n > 0);
     assert(pop.gametes[hap].n <= (N1 + N2));
     }
     */
    
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
    fwdpp_internal::process_gametes(pop.gametes, pop.mutations, pop.mcounts);
    assert(pop.mcounts.size() == pop.mutations.size());
    
    // DEBUG
    /*
     for (const auto &mc : pop.mcounts)
     {
     assert(mc <= (N1 + N2));
     }
     assert(happopdata_sane(pop.haploids, pop.gametes, pop.mutations, pop.mcounts));
     */
    
    /*
     The last thing to do is handle fixations.  In many contexts, we
     neither want nor need
     to keep indexes to fixed variants in our gametes.  Such decisions are
     implemented via
     simple policies, which are in the variable 'mp'.
     
     The implementation is in fwdpp/internal/gamete_cleaner.hpp.
     */
    fwdpp_internal::gamete_cleaner(pop.gametes, pop.mutations, pop.mcounts, (N1 + N2),
                                   mp);
}

//
// For structured deme, migration only from interpopulation recombination
//
template <typename structpoptype,
typename haploid_fitness_function, typename mutation_model,
typename recombination_policy,
typename mutation_removal_policy = std::true_type>
void
sample_haploid_struct_delayed_recmig(
                             //14 args
                             const gsl_rng *r,
                             structpoptype &pop,
                             //gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
                             //haploid_vector_type<haploid_geno_t, haploid_vector_type_allocator> &haploids,
                             //mutation_cont_type<mutation_type, mutation_cont_type_allocator> &mutations,
                             //std::vector<uint_t> &mcounts,
                             const uint_t &N1,
                             const uint_t &N2,
                             const double m12,
                             const double m21,
                             const double &mu,
                             const mutation_model &mmodel,
                             const recombination_policy &rec_pol,
                             const unsigned &gen,
                             const unsigned &ngenMig,
                             const haploid_fitness_function &ff,
                             //typename gamete_type::mutation_container &neutral,
                             //typename gamete_type::mutation_container &selected,
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
    /*
     assert(happopdata_sane(pop.haploids, pop.gametes, pop.mutations, pop.mcounts));
     assert(pop.mcounts.size() == pop.mutations.size());
     assert(check_sum(pop.gametes, (N1 + N2)));
     */
    
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
    auto mut_recycling_bin = fwdpp_internal::make_mut_queue(pop.mcounts);
    // this tags gametes with n=0 for recycling; b4 they're marked as 0 during fitness calc
    auto gam_recycling_bin = fwdpp_internal::make_gamete_queue(pop.gametes);
    
    // Migration and build pop-specific lookup tables:
    auto lookups = calc_fitness_twopop(r, pop, ff, N1, N2);
    
    // DEBUG
    /*
     for (const auto &g : pop.gametes)
     assert(!g.n);
     */
    
    const auto parents(pop.haploids); // Copy the parents
    
    // Fill in the next generation!
    // We generate the offspring for deme 1 first, and then for deme 2.
    // because of migration, the number of individuals for N1(N2) increases
    // in lookup table, but we only sample N1(N2) ind's from it,
    // conserving pop sizes.
    // Lookup tables from migrate_and_calc_fitness allow individuals
    // from N2 to get included while sampling for N1
    for (uint_t i = 0; i < N1 + N2; ++i)
    {
        std::size_t p1 = std::numeric_limits<std::size_t>::max();
        std::size_t p2 = std::numeric_limits<std::size_t>::max();
        if (i < N1) // pick parents from pop 1
        {
            // some of these may be indices from N2 if mig occurred
            p1 = lookups.parents1[gsl_ran_discrete(r, lookups.seln1.get())];
            
            if( gen > ngenMig && gsl_rng_uniform(r) < m21){
                // randomly pick rec. donor from N2 pop
                p2 = lookups.parents2[gsl_ran_discrete(r, lookups.rand2.get())];
            }else{
                // randomly pick rec. donor from same (N1) pop
                p2 = lookups.parents1[gsl_ran_discrete(r, lookups.rand1.get())];
            }
        }
        else // pick parents from pop 2
        {
            p1 = lookups.parents2[gsl_ran_discrete(r, lookups.seln2.get())];
            if( gen > ngenMig && gsl_rng_uniform(r) < m12){
                p2 = lookups.parents1[gsl_ran_discrete(r, lookups.rand1.get())];
            }else{
                p2 = lookups.parents2[gsl_ran_discrete(r, lookups.rand2.get())];
            }
        }
        assert(p1 < parents.size());
        assert(p2 < parents.size());
        
        auto h1 = parents[p1] ; // need to cp from parents,
        auto h2 = parents[p2] ; // haploids being modified
        
        // mutant/recombinant made from parent haplotypes h1+h2, assigned to
        // pop.haploids[i], transferring haplotypes associated with indices
        // from N2 to indices associated with N1
        mutate_recombine_update_haploid(r, pop.gametes, pop.mutations,
                                        std::make_tuple(h1, h2), rec_pol, mmodel,
                                        mu, gam_recycling_bin, mut_recycling_bin,
                                        pop.haploids[i], pop.neutral, pop.selected);
        
    }
    
    // DEBUG
    /*
     assert(check_sum(pop.gametes, (N1 + N2)));
     for (const auto &hap : pop.haploids)
     {
     assert(pop.gametes[hap].n > 0);
     assert(pop.gametes[hap].n <= (N1 + N2));
     }
     */
    
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
    fwdpp_internal::process_gametes(pop.gametes, pop.mutations, pop.mcounts);
    assert(pop.mcounts.size() == pop.mutations.size());
    
    // DEBUG
    /*
     for (const auto &mc : pop.mcounts)
     {
     assert(mc <= (N1 + N2));
     }
     assert(happopdata_sane(pop.haploids, pop.gametes, pop.mutations, pop.mcounts));
     */
    
    /*
     The last thing to do is handle fixations.  In many contexts, we
     neither want nor need
     to keep indexes to fixed variants in our gametes.  Such decisions are
     implemented via
     simple policies, which are in the variable 'mp'.
     
     The implementation is in fwdpp/internal/gamete_cleaner.hpp.
     */
    fwdpp_internal::gamete_cleaner(pop.gametes, pop.mutations, pop.mcounts, (N1 + N2),
                                   mp);
}


//
// For structured deme, migration only from interpopulation recombination
//
template <typename structpoptype,
typename haploid_fitness_function, typename mutation_model,
typename recombination_policy,
typename mutation_removal_policy = std::true_type>
void
sample_haploid_struct_disrupted_recmig(
                                     //14 args
                                     const gsl_rng *r,
                                     structpoptype &pop,
                                     //gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes,
                                     //haploid_vector_type<haploid_geno_t, haploid_vector_type_allocator> &haploids,
                                     //mutation_cont_type<mutation_type, mutation_cont_type_allocator> &mutations,
                                     //std::vector<uint_t> &mcounts,
                                     const uint_t &N1,
                                     const uint_t &N2,
                                     const double m12,
                                     const double m21,
                                     const double &mu,
                                     const mutation_model &mmodel,
                                     const recombination_policy &rec_pol,
                                     const unsigned &gen,
                                     const unsigned &ngenMig,
                                     const haploid_fitness_function &ff,
                                     //typename gamete_type::mutation_container &neutral,
                                     //typename gamete_type::mutation_container &selected,
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
    /*
     assert(happopdata_sane(pop.haploids, pop.gametes, pop.mutations, pop.mcounts));
     assert(pop.mcounts.size() == pop.mutations.size());
     assert(check_sum(pop.gametes, (N1 + N2)));
     */
    
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
    auto mut_recycling_bin = fwdpp_internal::make_mut_queue(pop.mcounts);
    // this tags gametes with n=0 for recycling; b4 they're marked as 0 during fitness calc
    auto gam_recycling_bin = fwdpp_internal::make_gamete_queue(pop.gametes);
    
    // Migration and build pop-specific lookup tables:
    auto lookups = calc_fitness_twopop(r, pop, ff, N1, N2);
    
    // DEBUG
    /*
     for (const auto &g : pop.gametes)
     assert(!g.n);
     */
    
    const auto parents(pop.haploids); // Copy the parents
    
    // Fill in the next generation!
    // We generate the offspring for deme 1 first, and then for deme 2.
    // because of migration, the number of individuals for N1(N2) increases
    // in lookup table, but we only sample N1(N2) ind's from it,
    // conserving pop sizes.
    // Lookup tables from migrate_and_calc_fitness allow individuals
    // from N2 to get included while sampling for N1
    for (uint_t i = 0; i < N1 + N2; ++i)
    {
        std::size_t p1 = std::numeric_limits<std::size_t>::max();
        std::size_t p2 = std::numeric_limits<std::size_t>::max();
        if (i < N1) // pick parents from pop 1
        {
            // some of these may be indices from N2 if mig occurred
            p1 = lookups.parents1[gsl_ran_discrete(r, lookups.seln1.get())];
            
            if( gen < ngenMig && gsl_rng_uniform(r) < m21){
                // randomly pick rec. donor from N2 pop
                p2 = lookups.parents2[gsl_ran_discrete(r, lookups.rand2.get())];
            }else{
                // randomly pick rec. donor from same (N1) pop
                p2 = lookups.parents1[gsl_ran_discrete(r, lookups.rand1.get())];
            }
        }
        else // pick parents from pop 2
        {
            p1 = lookups.parents2[gsl_ran_discrete(r, lookups.seln2.get())];
            if( gen < ngenMig && gsl_rng_uniform(r) < m12){
                p2 = lookups.parents1[gsl_ran_discrete(r, lookups.rand1.get())];
            }else{
                p2 = lookups.parents2[gsl_ran_discrete(r, lookups.rand2.get())];
            }
        }
        assert(p1 < parents.size());
        assert(p2 < parents.size());
        
        auto h1 = parents[p1] ; // need to cp from parents,
        auto h2 = parents[p2] ; // haploids being modified
        
        // mutant/recombinant made from parent haplotypes h1+h2, assigned to
        // pop.haploids[i], transferring haplotypes associated with indices
        // from N2 to indices associated with N1
        mutate_recombine_update_haploid(r, pop.gametes, pop.mutations,
                                        std::make_tuple(h1, h2), rec_pol, mmodel,
                                        mu, gam_recycling_bin, mut_recycling_bin,
                                        pop.haploids[i], pop.neutral, pop.selected);
        
    }
    
    // DEBUG
    /*
     assert(check_sum(pop.gametes, (N1 + N2)));
     for (const auto &hap : pop.haploids)
     {
     assert(pop.gametes[hap].n > 0);
     assert(pop.gametes[hap].n <= (N1 + N2));
     }
     */
    
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
    fwdpp_internal::process_gametes(pop.gametes, pop.mutations, pop.mcounts);
    assert(pop.mcounts.size() == pop.mutations.size());
    
    // DEBUG
    /*
     for (const auto &mc : pop.mcounts)
     {
     assert(mc <= (N1 + N2));
     }
     assert(happopdata_sane(pop.haploids, pop.gametes, pop.mutations, pop.mcounts));
     */
    
    /*
     The last thing to do is handle fixations.  In many contexts, we
     neither want nor need
     to keep indexes to fixed variants in our gametes.  Such decisions are
     implemented via
     simple policies, which are in the variable 'mp'.
     
     The implementation is in fwdpp/internal/gamete_cleaner.hpp.
     */
    fwdpp_internal::gamete_cleaner(pop.gametes, pop.mutations, pop.mcounts, (N1 + N2),
                                   mp);
}


#endif
