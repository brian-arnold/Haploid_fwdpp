#ifndef __GAMETEPOP_HPP__
#define __GAMETEPOP_HPP__

/*
  A structure representing a single-locus population.
  The user initizializes it with a population size, N
*/
#include <fwdpp/sugar/poptypes/tags.hpp>
#include <fwdpp/sugar/poptypes/popbase.hpp>
#include <fwdpp/sugar/poptypes/popbase.hpp>
#include <misc_functions.hpp>
#include <numeric>

using namespace fwdpp::sugar ;

/*!
  \brief Abstraction of what is needed to simulate a
  single-locus population.

  All that is missing is the mutation_type and the container types.

  See @ref md_md_sugar for rationale, etc.

  \ingroup sugar
*/
template <typename mutation_type, typename mcont, typename gcont,
            typename mvector, typename ftvector,
            typename lookup_table_type>
class gametepop : public popbase<mutation_type, mcont, gcont, mvector,
                                 ftvector, lookup_table_type>
{
  private:

    void
    process_individual_input()
    {
    }
   
    
    /*
    void
    process_individual_input()
    {
        for (auto &&g : this->gametes)
            {
                this->validate_gamete_keys(g);
            }
    }
    
    void
    validate_gamete_keys(const std::size_t key)
    {
        if (key >= this->gametes.size())
        {
            throw std::out_of_range("individual contains out of range keys");
        }
    }
     */
  public:
    virtual ~gametepop() = default;
    gametepop(gametepop &&) = default;
    gametepop(const gametepop &) = default;
    gametepop &operator=(gametepop &&) = default;
    gametepop &operator=(const gametepop &) = default;
    //! Population size
    uint_t N;

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

    //! Constructor
    explicit gametepop(
        const uint_t &popsize,
        typename popbase_t::gamete_t::mutation_container::size_type
            reserve_size
        = 100)
        : popbase_t(popsize, reserve_size), N(popsize)
    {
    }

    //! Constructor for pre-determined population status
    // need function to determine pop size, cycle thru gametes
    template <typename gametes_input,typename mutations_input>
    explicit gametepop(gametes_input &&g, mutations_input &&m)
        : popbase_t(std::forward<gametes_input>(g),
                    std::forward<mutations_input>(m), 100),
                    N{ static_cast<decltype(N)>( get_sum(g)) }
    {
        //this->process_individual_input();
    }

    bool
    operator==(const gametepop &rhs) const
    {
        // is_equal compares gametes, mutations, etc...
        // previous function also included diploids
        return popbase_t::is_equal(rhs);
    }

    //! Empty all the containers
    void
    clear()
    {
        popbase_t::clear_containers();
    }
    
    /*
    //template<typename gcont>
    uint_t get_sum(const gcont_t &gametes)
    {
        uint_t sum ;
        uint_t init = 0 ;
        sum = std::accumulate(gametes.cbegin(), gametes.cend(), init,
                              [](uint_t init,
                                 const typename gcont_t::value_type &__g) {
                                  return init + __g.n;
                              });
        return sum ;
    }
     */
};


#endif
