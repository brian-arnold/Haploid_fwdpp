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

//
// FITNESS FUNCTIONS FOR SINGLE POPULATION
//
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

struct multiplicative_posseln_haploid
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
            product *= (1.0 + mutations[key].s) ;
        }
        return std::max(0., product);
    }
};

class quadratic_synepistasis_negseln
{
private:
    double b ; // pairwise epistasis coefficient
    double s ; // selection coefficient
public:
    // Constructor
    explicit quadratic_synepistasis_negseln(const double sel, const double beta) :
            s(sel), b(beta)
    {
    }
    
    template<typename gamete_type, typename mtype>
    // function has to be const! in template declaration, functor labelled as const
    // thus, any attempt to change a member variable or call non-const member function
    // results in compiler error
    inline double
    operator()(const gamete_type &g, const std::vector<mtype> &mutations, const std::vector<uint_t> &mcounts) const noexcept
    {
        double muts = 0.0 ;
        double fitness ;
        for(const std::uint32_t &key : g.smutations){
            muts += 1.0 ;
        }
        //fitness = 1.0 - (s*muts) - b*(muts*(muts-1.0)/2.0) ;
        fitness = exp(-(s*muts) - b*(muts*(muts-1.0)/2.0)) ;
        //std::cout << fitness << "\n" ;
        return std::max(0., fitness);
    }
};

class antagonistic_epistasis_negseln
{
private:
    double b ; // pairwise epistasis coefficient
    double s ; // selection coefficient
public:
    // Constructor
    explicit antagonistic_epistasis_negseln(const double sel, const double beta) :
            s(sel), b(beta)
    {
    }
    
    template<typename gamete_type, typename mtype>
    // function has to be const! in template declaration, functor labelled as const
    // thus, any attempt to change a member variable or call non-const member function
    // results in compiler error
    inline double
    operator()(const gamete_type &g, const std::vector<mtype> &mutations, const std::vector<uint_t> &mcounts) const noexcept
    {
        double muts = 0.0 ;
        double fitness ;
        for(const std::uint32_t &key : g.smutations){
            muts += 1.0 ;
        }
	// From Desai et al 2007, Genetics
	// form: e^(-s*muts^(1-b))
	// negative b (or minus a positive b) corresponds to antagonistic epistasis
        fitness = exp( pow(-(s*muts), (1.0-b)) ) ;
        //std::cout << fitness << "\n" ;
        return std::max(0., fitness);
    }
};

class synergistic_epistasis_posseln
{
private:
    double b ; // pairwise epistasis coefficient
    double s ; // selection coefficient
public:
    // Constructor
    explicit synergistic_epistasis_posseln(const double sel, const double beta) :
            s(sel), b(beta)
    {
    }
    
    template<typename gamete_type, typename mtype>
    // function has to be const! in template declaration, functor labelled as const
    // thus, any attempt to change a member variable or call non-const member function
    // results in compiler error
    inline double
    operator()(const gamete_type &g, const std::vector<mtype> &mutations, const std::vector<uint_t> &mcounts) const noexcept
    {
        double muts = 0.0 ;
        double fitness ;
        for(const std::uint32_t &key : g.smutations){
            muts += 1.0 ;
        }
	// From Desai et al 2007, Genetics
	// form: e^(s*muts^(1+b))
	// positive b corresponds to synergistic epistasis
        fitness = exp( pow((s*muts), (1.0+b)) ) ;
        //std::cout << fitness << "\n" ;
        return std::max(0., fitness);
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
            product *= (1.0 + mutations[key].s*((eq - static_cast<double>(mcounts[key])/N))) ;
            //product *= pow( (1.0 + mutations[key].s), (eq - static_cast<double>(mcounts[key])/N)/eq ) ;
        }
        return std::max(0., product);
    }
};



//
// FITNESS FUNCTIONS FOR 2 POPULATIONS
//
struct pop_sign_seln_multiplicative
{
    // demes have opposite selection pressures, positive in deme0, negative in deme1
    template<typename gamete_type, typename mtype>
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

class pop_sign_seln_epistasis
{
private:
    double b ; // pairwise epistasis coefficient
    double s ; // selection coefficient
public:
    // Constructor
    explicit pop_sign_seln_epistasis(const double sel, const double beta) :
    s(sel), b(beta)
    {
    }
    
    template<typename gamete_type, typename mtype>
    // function has to be const! in template declaration, functor labelled as const
    // thus, any attempt to change a member variable or call non-const member function
    // results in compiler error
    inline double
    operator()(const gamete_type &g, const std::vector<mtype> &mutations, const int deme) const noexcept
    {
        double muts = 0.0 ;
        double fitness ;
        for(const std::uint32_t &key : g.smutations){
            muts += 1.0 ;
        }
        //fitness = 1.0 - (s*muts) - b*(muts*(muts-1.0)/2.0) ;
	if(deme){
        	fitness = exp(-(s*muts) - b*(muts*(muts-1.0)/2.0)) ;
	}else{
        	fitness = exp((s*muts) + b*(muts*(muts-1.0)/2.0)) ;
	}
        //std::cout << fitness << "\n" ;
        return std::max(0., fitness);
    }
};



struct pop_specific_seln_multiplicative
{
    // deme0 under pos selection, deme1 neutral
    template<typename gamete_type, typename mtype>
    inline double
    operator()(const gamete_type &g, const std::vector<mtype> &mutations, const int deme) const noexcept
    {
        double product = 1.0 ;
        for(const std::uint32_t &key : g.smutations){
            if(deme){   // neutral in N2
                // neutral, do nothing
            }else{      // pos seln in N1
                product *= (1.0 + mutations[key].s) ;
            }
        }
        return std::max(0., product);
    }
};

struct pop_pos_seln_multiplicative
{
    // both demes under multiplicative positive selection
    template<typename gamete_type, typename mtype>
    inline double
    operator()(const gamete_type &g, const std::vector<mtype> &mutations, const int deme) const noexcept
    {
        double product = 1.0 ;
        for(const std::uint32_t &key : g.smutations){
            product *= (1.0 + mutations[key].s) ;
        }
        return std::max(0., product);
    }
};

struct pop_neg_seln_multiplicative
{
    // both demes under multiplicative negative selection
    template<typename gamete_type, typename mtype>
    inline double
    operator()(const gamete_type &g, const std::vector<mtype> &mutations, const int deme) const noexcept
    {
        double product = 1.0 ;
        for(const std::uint32_t &key : g.smutations){
            product *= (1.0 - mutations[key].s) ;
        }
        return std::max(0., product);
    }
};


class pop_pos_seln_synepistasis
{
private:
    double b ; // pairwise epistasis coefficient
    double s ; // selection coefficient
public:
    // Constructor
    explicit pop_pos_seln_synepistasis(const double sel, const double beta) :
    s(sel), b(beta)
    {
    }
    
    template<typename gamete_type, typename mtype>
    // function has to be const! in template declaration, functor labelled as const
    // thus, any attempt to change a member variable or call non-const member function
    // results in compiler error
    inline double
    operator()(const gamete_type &g, const std::vector<mtype> &mutations, const int deme) const noexcept
    {
        double muts = 0.0 ;
        double fitness ;
        for(const std::uint32_t &key : g.smutations){
            muts += 1.0 ;
        }
        //fitness = 1.0 - (s*muts) - b*(muts*(muts-1.0)/2.0) ;
        fitness = exp((s*muts) + b*(muts*(muts-1.0)/2.0)) ;
        //std::cout << fitness << "\n" ;
        return std::max(0., fitness);
    }
};


/*
// GOAL HERE IS TO DO A POP SPECIFIC NFDS, BUT YOU'D NEED TO CHANGE OTHER FUNCTIONS TO ALSO TAKE
 // MCOUNTS ARG UNLESS YOU WANT TO MAKE A SEPARATE SAMPLING FUNCTION JUST FOR NFDS
// BEFORE THIS CAN BE IMPLEMENTED, YOU NEED A STRUCTURE THAT KEEPS TRACK OF MCOUNTS PER POPULATION
class nfds_two_pop
{
private:
    double N1 ; // populaiton size, double b/c
    double N2 ;
    double eq ; // equilibrium frequency of selected mutations
public:
    
    // Constructor
    explicit nfds(const uint_t &popsize1, const uint_t &popsize2, const double equilibrium_freq) :
    N1(static_cast<double>(popsize1)), N2(static_cast<double>(popsize2)), eq(equilibrium_freq)
    {
    }
    
    template<typename gamete_type, typename mtype>
    // function has to be const! in template declaration, functor labelled as const
    // thus, any attempt to change a member variable or call non-const member function
    // results in compiler error
    inline double
    operator()(const gamete_type &g, const std::vector<mtype> &mutations, const std::vector<uint_t> &mcounts, const int deme) const noexcept
    {
        double product = 1.0 ;
        for(const std::uint32_t &key : g.smutations){
            if(deme==0){ // 
                product *= (1.0 + mutations[key].s*((eq - static_cast<double>(mcounts[key])/N1)/eq)) ;
            }else{
                product *= (1.0 + mutations[key].s*((eq - static_cast<double>(mcounts[key])/N2)/eq)) ;
            }
        }
        return std::max(0., product);
    }
};
*/

#endif
