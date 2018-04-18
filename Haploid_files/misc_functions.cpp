#ifndef _MISC_FUNCTIONS_CPP_
#define _MISC_FUNCTIONS_CPP_

#include <vector>
#include <algorithm>
#include <fwdpp/forward_types.hpp>
#include <numeric>
#include <functional>

using namespace fwdpp ;

//template<typename gcont>
uint_t get_sum(const std::vector<gamete> &gametes)
{
    uint_t sum = 0 ;
    //uint_t init = 0 ;
 
    for (const auto &g : gametes)
    {
        sum += g.n ;
    }
    
    /*
    sum = std::accumulate(gametes.cbegin(), gametes.cend(), init,
                          [](uint_t init,
                             const typename gcont::value_type &__g) {
                              return init + __g.n;
                          });
    */
    
    return sum ;
}

#endif


