#ifndef _FITNESS_MODELS_HAPLOID_CPP_
#define _FITNESS_MODELS_HAPLOID_CPP_

#include <vector>
#include <algorithm>
#include <fwdpp/forward_types.hpp>
#include <numeric>
#include <functional>

using namespace fwdpp ;

// gamete defined in <fwdpp/forward_types.hpp>
/*
template <typename mtype>
double multiplicative_haploid(const gamete &g, const std::vector<mtype> &mutations)
{
     double product = 1.0 ;
     for(const std::uint32_t &key : g.smutations){
        product *= (1.0 + mutations[key].s) ;
     }
    return std::max(0., product);
}
*/

/*
struct multiplicative_haploid
{
    double product = 1.0 ;
    
    template<typename gamete_type, typename mtype>
    //template<typename mtype>
   inline double
    operator()(const gamete_type &g, const std::vector<mtype> &mutations) noexcept
    {
        for(const std::uint32_t &key : g.smutations){
            product *= (1.0 + mutations[key].s) ;
        }
        return std::max(0., product);
    }
};
*/

#endif


