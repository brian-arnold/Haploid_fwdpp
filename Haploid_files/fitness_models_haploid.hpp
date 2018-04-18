#ifndef _FITNESS_MODELS_HAPLOID_HPP_
#define _FITNESS_MODELS_HAPLOID_HPP_

#include <vector>
#include <algorithm>
#include <fwdpp/forward_types.hpp>
#include <numeric>
#include <functional>

using namespace fwdpp ;

struct multiplicative_haploid
{
    
    template<typename gamete_type, typename mtype>
    // function has to be const! in template declaration, functor labelled as const
    // thus, any attempt to change a member variable or call non-const member function
    // results in compiler error
    inline double
    operator()(const gamete_type &g, const std::vector<mtype> &mutations) const noexcept
    {
        double product = 1.0 ;
        for(const std::uint32_t &key : g.smutations){
            product *= (1.0 + mutations[key].s) ;
        }
        return std::max(0., product);
    }
};

#endif
