#ifndef _MISC_FUNCTIONS_HPP_
#define _MISC_FUNCTIONS_HPP_

#include <vector>
#include <algorithm>
#include <fwdpp/forward_types.hpp>
#include <numeric>
#include <functional>
#include "misc_functions.cpp"

using namespace fwdpp ;

uint_t get_sum(const std::vector<gamete> &gametes) ;


template <typename hapcont_t, typename gcont_t,
            typename mcont_t>
bool
happopdata_sane(const hapcont_t &haploids, const gcont_t &gametes,
             const mcont_t &mutations,
             const std::vector<uint_t> &mutcounts) ;


template <typename gamete_type,
            typename gamete_cont_type_allocator,
            template <typename, typename> class gamete_cont_type>
std::vector< std::pair<std::size_t, std::size_t> >
group_haps_into_dips(const gsl_rng *r,
                     gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes) ;


#endif


