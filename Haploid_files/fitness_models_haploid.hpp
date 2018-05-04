#ifndef _FITNESS_MODELS_HAPLOID_HPP_
#define _FITNESS_MODELS_HAPLOID_HPP_

#include <vector>
#include <algorithm>
#include <fwdpp/forward_types.hpp>
#include <numeric>
#include <functional>
#include <math.h>
#include <iostream>

using namespace fwdpp ;

struct multiplicative_negseln_haploid
{
    
    template<typename gamete_type, typename mtype>
    // function has to be const! in template declaration, functor labelled as const
    // thus, any attempt to change a member variable or call non-const member function
    // results in compiler error
    inline double
    operator()(const gamete_type &g, const std::vector<mtype> &mutations, const std::vector<uint_t> &mcounts) const noexcept
    {
        double product = 1.0 ;
        for(const std::uint32_t &key : g.smutations){
            product *= (1.0 - mutations[key].s) ;
        }
        return std::max(0., product);
    }
};


class nfds
{
  private:
    double N ; // populaiton size, double b/c
    double eq ; // equilibrium frequency of selected mutations
  public:
    
    // Constructor
    explicit nfds(const uint_t &popsize, const double equilibrium_freq) :
        N(static_cast<double>(popsize)), eq(equilibrium_freq)
    {        
    }
    
    template<typename gamete_type, typename mtype>
    // function has to be const! in template declaration, functor labelled as const
    // thus, any attempt to change a member variable or call non-const member function
    // results in compiler error
    inline double
    operator()(const gamete_type &g, const std::vector<mtype> &mutations, const std::vector<uint_t> &mcounts) const noexcept
    {
        double product = 1.0 ;
        for(const std::uint32_t &key : g.smutations){
            product *= (1.0 + mutations[key].s*((eq - static_cast<double>(mcounts[key])/N)/eq)) ;
            //product *= pow( (1.0 + mutations[key].s), (eq - static_cast<double>(mcounts[key])/N)/eq ) ;
        }
        return std::max(0., product);
    }
};



struct pop_sign_seln_multiplicative
{
    
    template<typename gamete_type, typename mtype>
    // function has to be const! in template declaration, functor labelled as const
    // thus, any attempt to change a member variable or call non-const member function
    // results in compiler error
    inline double
    operator()(const gamete_type &g, const std::vector<mtype> &mutations, const int deme) const noexcept
    {
        double product = 1.0 ;
        for(const std::uint32_t &key : g.smutations){
            if(deme){   // neg seln in N2
                product *= (1.0 - mutations[key].s) ;
            }else{      // pos seln in N1
                product *= (1.0 + mutations[key].s) ;
            }
        }
        return std::max(0., product);
    }
};

#endif
