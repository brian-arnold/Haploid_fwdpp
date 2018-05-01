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
    for (const auto &g : gametes)
    {
        sum += g.n ;
    }
    return sum ;
}


template <typename hapcont_t, typename gcont_t, typename mcont_t>
bool
happopdata_sane(const hapcont_t &haploids, const gcont_t &gametes,
             const mcont_t &mutations,
             const std::vector<uint_t> &mutcounts)
/*
 \brief Check that all haploids refer to extant gametes with sane data.
 */
{
    for (const auto &h : haploids)
    {
        if (!gametes[h].n)
            return false;
        if (!gamete_data_sane(gametes[h], mutations, mutcounts))
            return false;
    }
    return true;
}


template <typename gamete_type,
            typename gamete_cont_type_allocator,
            template <typename, typename> class gamete_cont_type>
std::vector< std::pair<std::size_t, std::size_t> >
group_haps_into_dips(const gsl_rng *r,
                     gamete_cont_type<gamete_type, gamete_cont_type_allocator> &gametes)
{
    
    using diploid_t = std::pair<std::size_t, std::size_t> ;
    using dipvector_t = std::vector<diploid_t> ;

    int popsize = get_sum(gametes) ;
    int gcounts[ popsize ] ;
    int k = 0 ;
    for (int i = 0; i < gametes.size(); i++){
        if(gametes[i].n){
            for(int j=0; j<gametes[i].n; j++){
                gcounts[k] = i ;
                k++ ;
            }
        }
    }
    /*
    for (int i=0; i < k; i++){
        std::cout << gcounts[i] << " " ;
    }
    std::cout << "\n" ;
    */
    //shuffle gcounts  so that "diploids" represent random pairs
    //double randomize for good measure? probably unecessary
    gsl_ran_shuffle( r, gcounts, k, sizeof(int) ) ;
    gsl_ran_shuffle( r, gcounts, k, sizeof(int) ) ;
    /*
    for (int i=0; i < k; i++){
        std::cout << gcounts[i] << " " ;
    }
    std::cout << "\n" ;
    */
    //go thru and put pairs of gametes into diploids
    dipvector_t pseudodips((k/2), diploid_t(0,0)) ;
    for(int i=0; i<k-1; i+=2){
        pseudodips[(i/2)].first = gcounts[i] ;
        pseudodips[(i/2)].second = gcounts[(i+1)] ;
   }
    /*
    for (int i=0; i<pseudodips.size(); i++) {
        std::cout << pseudodips[i].first << "\t" ;
        std::cout << pseudodips[i].second << "\n" ;
    }
    */
    
    return pseudodips ;
}



#endif


