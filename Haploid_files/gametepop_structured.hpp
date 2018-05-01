#ifndef __GAMETEPOP_STRUCTURED_HPP__
#define __GAMETEPOP_STRUCTURED_HPP__

/*
 This class is virtually identical to gametepop
 with the exception of having a vector "haploids"
 that refers to gametes, which are no longer exchangeable
 in a structured population since a gamete with high
 frequency could be distributed amongst multiple populations.
*/
#include <fwdpp/sugar/poptypes/tags.hpp>
#include <fwdpp/sugar/poptypes/popbase.hpp>
#include <fwdpp/sugar/poptypes/popbase.hpp>
#include "misc_functions.hpp"
#include <numeric>

using namespace fwdpp::sugar ;

template <typename mutation_type, typename mcont, typename gcont,
            typename hapvector, typename mvector, typename ftvector,
            typename lookup_table_type>
class gametepopstruct : public popbase<mutation_type, mcont, gcont, mvector,
                                 ftvector, lookup_table_type>
{
  private:
    void
    process_individual_input()
    {
        std::vector<uint_t> gcounts(this->gametes.size(), 0);
        for (auto &&hap : haploids)
        {
            this->validate_individual_keys(hap);
            gcounts[hap]++;
        }
        this->validate_gamete_counts(gcounts);
    }
    
  public:
    virtual ~gametepopstruct() = default;
    gametepopstruct(gametepopstruct &&) = default;
    gametepopstruct(const gametepopstruct &) = default;
    gametepopstruct &operator=(gametepopstruct &&) = default;
    gametepopstruct &operator=(const gametepopstruct &) = default;
    //! Population size
    uint_t N;

    using hapvector_t = hapvector;
    using haploid_t = typename hapvector::value_type;
    //! Typedef for base class
    using popbase_t = popbase<mutation_type, mcont, gcont, mvector,
                              ftvector, lookup_table_type>;
    //! Dispatch tag for other parts of sugar layer
    using popmodel_t = sugar::SINGLELOC_TAG;
    //! Fitness function signature compatible with this type
    //using fitness_t
    //    = fwdpp::traits::fitness_fxn_t<dipvector_t,
    //                                   typename popbase_t::gcont_t,
    //                                   typename popbase_t::mcont_t>;

    //! Container of haploids
    hapvector_t haploids;
    
    //! Constructor
    explicit gametepopstruct(
        const uint_t &popsize,
        typename popbase_t::gamete_t::mutation_container::size_type
            reserve_size
            = 100)
        : popbase_t(popsize, reserve_size), N(popsize),
            // All N haploids contain the only gamete in the pop
            haploids(hapvector_t(popsize, haploid_t(0)))
    {
    }

    //! Constructor for pre-determined population status
    // need function to determine pop size, cycle thru gametes
    template <typename haploids_input, typename gametes_input,
                typename mutations_input>
    explicit gametepopstruct(haploids_input &&h, gametes_input &&g,
                             mutations_input &&m)
        : popbase_t(std::forward<gametes_input>(g),
                    std::forward<mutations_input>(m), 100),
                    N{ static_cast<decltype(N)>( h.size()) },
                    haploids(std::forward<haploids_input>(h))
     //! Constructor for pre-determined population status
    {
        this->process_individual_input();
    }

    bool
    operator==(const gametepopstruct &rhs) const
    {
        return this->haploids == rhs.haploids
        && popbase_t::is_equal(rhs);
    }

    //! Empty all the containers
    void
    clear()
    {
        haploids.clear();
        popbase_t::clear_containers();
    }

};


#endif
