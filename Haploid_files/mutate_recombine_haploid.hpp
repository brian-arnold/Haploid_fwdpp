/*!
  \file mutate_recombine.hpp

  \brief Handling of recombination and mutation in one step.

  \note Introduced in fwdpp 0.5.7
*/
#ifndef _MUTATE_RECOMBINE_HAPLOID_HPP__
#define _MUTATE_RECOMBINE_HAPLOID_HPP__

#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <tuple>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fwdpp/type_traits.hpp>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/internal/mutation_internal.hpp>
#include <fwdpp/internal/rec_gamete_updater.hpp>
#include <fwdpp/mutate_recombine.hpp>

using namespace fwdpp ;

template <typename recombination_policy, typename gcont_t,
            typename mcont_t>
std::vector<double>
generate_breakpoints_haploid(const std::size_t g1, const std::size_t g2,
                     const gcont_t &gametes, const mcont_t &mutations,
                     const recombination_policy &rec_pol)
/// Generate vector of recombination breakpoints
///
/// \param diploid A single-locus diploid.
/// \param g1 Index of gamete 1
/// \param g2 Index of gamete 2
/// \param gametes Vector of gametes
/// \param mutations Vector of mutations
/// \param rec_pol Function to generate breakpoints
///
/// \return std::vector<double> containing sorted breakpoints
///
/// \note An empty return value means no breakpoints.  Otherwise,
/// the breakpoints are returned and are terminated by
/// std::numeric_limits<double>::max()
{
    auto nm1
        = gametes[g1].mutations.size() + gametes[g1].smutations.size();
    auto nm2
        = gametes[g2].mutations.size() + gametes[g2].smutations.size();
    if ((std::min(nm1, nm2) == 0 && std::max(nm1, nm2) == 1)
        || gametes[g1] == gametes[g2])
        {
            return {};
        }
    return rec_pol() ;
}




template <typename queue_type, typename mutation_model,
          typename gcont_t, typename mcont_t>
std::vector<uint_t>
generate_new_mutations_haploid(queue_type &recycling_bin,
                       const gsl_rng *r,
                       const double &mu,
                       gcont_t &gametes,
                       mcont_t &mutations,
                       const std::size_t g,
                       const mutation_model &mmodel)
///
/// Return a vector of keys to new mutations.  The keys
/// will be sorted according to mutation postition.
///
/// \param recycling_bin The queue for recycling mutations
/// \param r A random number generator
/// \param mu The total mutation rate
/// \param dip A single-locus diploid
/// \param gametes Vector of gametes
/// \param mutations Vector of mutations
/// \param g index of gamete to mutate
/// \param mmodel The mutation policy
///
/// \return Vector of mutation keys, sorted according to position
///
{
    unsigned nm = gsl_ran_poisson(r, mu);
    std::vector<uint_t> rv;
    rv.reserve(nm);
    for (unsigned i = 0; i < nm; ++i)
        {
            //BA: mmodel creates position for new mmutation, feeds to
            //BA: fwdpp_internal::recycle_mutation_helper, which returns
            //BA: index to new mut in mutation cont
            rv.emplace_back( mmodel(recycling_bin, mutations) );
        }
    std::sort(rv.begin(), rv.end(),
              [&mutations](const uint_t a, const uint_t b) {
                  return mutations[a].pos < mutations[b].pos;
              });
    return rv;
}


template <typename gcont_t, typename mcont_t,
          typename recmodel, typename mutmodel>
std::tuple<std::size_t, std::size_t>
mutate_recombine_update_haploid(
    //11 args
    const gsl_rng *r,
    gcont_t &gametes,
    mcont_t &mutations,
    std::tuple<std::size_t, std::size_t> parental_gametes,
    const recmodel &rec_pol,
    const mutmodel &mmodel,
    const double mu,
    typename traits::recycling_bin_t<gcont_t> &gamete_recycling_bin,
    typename traits::recycling_bin_t<mcont_t> &mutation_recycling_bin,
    typename gcont_t::value_type::mutation_container &neutral,
    typename gcont_t::value_type::mutation_container &selected)
///
/// "Convenience" function for generating offspring gametes.
///
/// This function calls fwdpp::generate_breakpoints,
/// fwdpp::generate_new_mutations, and fwdpp::mutate_recombine,
/// resulting in two offspring gametes.
///
/// \param r A gsl_rng *
/// \param gametes Vector of gametes in population
/// \param mutations Vector of mutations in population
/// \param parental_gametes Tuple of gamete keys for each parent
/// \param rec_pol Policy to generate recombination breakpoints
/// \param mmodel Policy to generate new mutations
/// \param mu Total mutation rate (per gamete).
/// \param gamete_recycling_bin FIFO queue for gamete recycling
/// \param mutation_recycling_bin FIFO queue for mutation recycling
/// \param neutral Temporary container for updating neutral mutations
/// \param selected Temporary container for updating selected mutations
///
/// \return Number of recombination breakpoints and mutations for each
/// gamete.
///
/// \version
/// Added in fwdpp 0.5.7.
{
    auto p1 = std::get<0>(parental_gametes);
    auto p2 = std::get<1>(parental_gametes);
    // Now, we generate the crossover breakpoints for
    // both parents,as well as the new mutations that we'll place
    // onto each gamete.  The specific order of operations below
    // is done to ensure the exact same output as fwdpp 0.5.6 and
    // earlier.
    // The breakpoints are of type std::vector<double>, and
    // the new_mutations are std::vector<fwdpp::uint_t>, with
    // the integers representing the locations of the new mutations
    // in "mutations".
    
    auto breakpoints
            = generate_breakpoints_haploid(p1, p2, gametes, mutations,
                                           rec_pol);
    auto new_mutations
        = generate_new_mutations_haploid(mutation_recycling_bin, r, mu,
                                 gametes, mutations, p1, mmodel);

    // Pass the breakpoints and new mutation keys on to
    // fwdpp::mutate_recombine (defined in
    // fwdpp/mutate_recombine.hpp),
    // which splices together the offspring gamete and returns its
    // location in gametes.  The location of the offspring gamete
    // is
    // either reycled from an extinct gamete or it is the location
    // of a
    // new gamete emplace_back'd onto the end.
    
    
    auto p1newIdx = fwdpp::mutate_recombine(new_mutations, breakpoints, p1, p2,
                                 gametes, mutations, gamete_recycling_bin,
                                 neutral, selected);
    gametes[p1newIdx].n++;
    assert(gametes[p1newIdx].n);
    
    return std::make_tuple( breakpoints.size(), new_mutations.size() );

    
}



// Overloaded function that also takes container of haploids
// containing keys to gametes, which aren't exchangeable with structure
template <typename haploid_t, typename gcont_t, typename mcont_t,
typename recmodel, typename mutmodel>
std::tuple<std::size_t, std::size_t>
mutate_recombine_update_haploid(
                                //11 args
                                const gsl_rng *r,
                                gcont_t &gametes,
                                mcont_t &mutations,
                                std::tuple<std::size_t, std::size_t> parental_gametes,
                                const recmodel &rec_pol,
                                const mutmodel &mmodel,
                                const double mu,
                                typename traits::recycling_bin_t<gcont_t> &gamete_recycling_bin,
                                typename traits::recycling_bin_t<mcont_t> &mutation_recycling_bin,
                                haploid_t &hap,
                                typename gcont_t::value_type::mutation_container &neutral,
                                typename gcont_t::value_type::mutation_container &selected)
///
/// "Convenience" function for generating offspring gametes.
///
/// This function calls fwdpp::generate_breakpoints,
/// fwdpp::generate_new_mutations, and fwdpp::mutate_recombine,
/// resulting in two offspring gametes.
///
/// \param r A gsl_rng *
/// \param gametes Vector of gametes in population
/// \param mutations Vector of mutations in population
/// \param parental_gametes Tuple of gamete keys for each parent
/// \param rec_pol Policy to generate recombination breakpoints
/// \param mmodel Policy to generate new mutations
/// \param mu Total mutation rate (per gamete).
/// \param gamete_recycling_bin FIFO queue for gamete recycling
/// \param mutation_recycling_bin FIFO queue for mutation recycling
/// \param hap The offspring
/// \param neutral Temporary container for updating neutral mutations
/// \param selected Temporary container for updating selected mutations
///
/// \return Number of recombination breakpoints and mutations for each
/// gamete.
///
///
/// \version
/// Added in fwdpp 0.5.7.
{
    auto p1 = std::get<0>(parental_gametes);
    auto p2 = std::get<1>(parental_gametes);
    // The breakpoints are of type std::vector<double>, and
    // the new_mutations are std::vector<fwdpp::uint_t>, with
    // the integers representing the locations of the new mutations
    // in "mutations".
    
    auto breakpoints
    = generate_breakpoints_haploid(p1, p2, gametes, mutations,
                                   rec_pol);
    auto new_mutations
    = generate_new_mutations_haploid(mutation_recycling_bin, r, mu,
                                     gametes, mutations, p1, mmodel);
    
    // Pass the breakpoints and new mutation keys on to
    // fwdpp::mutate_recombine (defined in
    // fwdpp/mutate_recombine.hpp),
    // which splices together the offspring gamete and returns its
    // location in gametes.  The location of the offspring gamete
    // is
    // either reycled from an extinct gamete or it is the location
    // of a
    // new gamete emplace_back'd onto the end.
    
    
    hap = fwdpp::mutate_recombine(new_mutations, breakpoints, p1, p2,
                                            gametes, mutations, gamete_recycling_bin,
                                            neutral, selected);
    gametes[hap].n++;
    assert(gametes[hap].n);
    
    return std::make_tuple( breakpoints.size(), new_mutations.size() );
    
    
}



#endif
