/*  \include juvenile_migration.cc
 *
 *  Juvenile migration.
 *  Selection.
 *  Constant pop size.
 *  Unequal migration rates.
 *
 * The goal is to show how to use a population
 * object to handle details of population structure.
 *
 * The outline of the scheme is as follows:
 *
 * 1. Offspring are generated in sorted order by deme
 * label.  In other words, deme 1 before deme 2.
 *
 * 2. Each generation, we migrate first.  Thus,
 * we may start a generation with parents like this:
 *
 * Parent Deme
 * 0    0
 * 1    0
 * 2    0
 * 3    1
 * 4    1
 * 5    1
 *
 * After migration, we may end up with:
 *
 * Parent Deme
 * 0    1*
 * 1    0
 * 2    0
 * 3    1
 * 4    0*
 * 5    1
 *
 * An asterisk represent migrant individuals.
 *
 * We will generate efficient lookup tables to sample individuals
 * proportional to their within-deme fitness, post-migration.
 *
 * These lookup tables will map to the indexes of the parents of each
 * deme:
 *
 * Lookup 1:
 *
 * Index Parent
 * 0    1
 * 1    2
 * 2    4
 *
 * Lookup 2:
 * Index Parent
 * 0    0
 * 1    3
 * 1    5
 *
 * The code is not safe for real-world use.  It doesn't handle the corner
 * case of all of deme 1 migrating into deme 2, leaving deme 1 "extinct"
 * in the next generation.  Such edge cases are clearly important, but
 * the goal here it to show the book-keeping of parental fitnesses
 * and the mapping back to deme labels.
 */
#include <fwdpp/recbinder.hpp>
#include <numeric>
#include <functional>
#include <cassert>
#include <unordered_set>
#include <gsl/gsl_randist.h>



using namespace fwdpp;

struct parent_lookup_tables
// Our object for choosing parents each generation
{
    // These return indexes of parents from demes 1 and 2,
    // resp, chosen in O(1) time proportional to
    // relative fitness within each deme
    fwdpp_internal::gsl_ran_discrete_t_ptr seln1, seln2, rand1, rand2;
    // These vectors map indexes returned from sampling
    // lookup1 and lookup2 to diploids in the population
    // object.
    std::vector<std::size_t> parents1, parents2;
};

template <typename structpoptype, typename haploid_fitness_function>
parent_lookup_tables
migrate_and_calc_fitness(const gsl_rng *r, structpoptype &pop,
                         const haploid_fitness_function &ff, const uint_t N1,
                         const uint_t N2, const double m12, const double m21)
// This function will be called at the start of each generation.
// The main goal is to return the lookup tables described above.
// But, "while we're at it", it does some other stuff that
// needs to be done at the start of each generation.
// Neither the most rigorous nor the most efficient:
// 1. Ignores probability of back-migration.
// 2. Allocates 4 vectors each generation.
{
    parent_lookup_tables rv;

    // Temp containers for fitnesses in each deme,
    // post-migration
    std::vector<double> w1, w2;

    // Pick no. migrants 1 -> 2 and 2 -> 1.
    unsigned nmig12 = gsl_ran_poisson(r, static_cast<double>(N1) * m12);
    unsigned nmig21 = gsl_ran_poisson(r, static_cast<double>(N2) * m21);
    // make sure migrants not >= pop size, NO EXTINCTIONS!
    while(nmig12 >= N1){
        nmig12 = gsl_ran_poisson(r, static_cast<double>(N1) * m12) ;
    }
    while(nmig21 >= N2){
        nmig21 = gsl_ran_poisson(r, static_cast<double>(N2) * m21) ;
    }
    // Fill a vector of N1 zeros and N2 ones:
    // The fist N1 indices represent pop1
    std::vector<uint_t> deme_labels(N1, 0);
    deme_labels.resize(N1 + N2, 1);
    assert(deme_labels.size() == pop.haploids.size());
    
    // Set up source and destination containers
    // for sampling w/o replacement
    std::vector<std::size_t> individuals(N1 + N2);
    // fills vector with sequentially increasing values
    std::iota(std::begin(individuals), std::end(individuals), 0);
    std::vector<std::size_t> migrants(std::max(nmig12, nmig21));

    
    // Who is migrating 1 -> 2?
    // gsl function fills migrants array with nmig12 objects from the first
    // N1 elements of individuals array
    gsl_ran_choose(r, migrants.data(), nmig12, individuals.data(), N1,
                   sizeof(std::size_t));
    for (std::size_t i = 0; i < nmig12; ++i)
        {
            deme_labels[migrants[i]] = !deme_labels[migrants[i]];
        }

    // Exact same logic for migrants 2 -> 1
    // migrants array gets recycled, adds nmig21 objects starting at beginning?
    gsl_ran_choose(r, migrants.data(), nmig21, individuals.data() + N1, N2,
                   sizeof(std::size_t));
    for (std::size_t i = 0; i < nmig21; ++i)
        {
            deme_labels[migrants[i]] = !deme_labels[migrants[i]];
        }
    
    // Go over all parents, set gametes counts to zero,
    // and put individual IDs and fitnesses into
    // the right vectors:
    for (std::size_t i = 0; i < deme_labels.size(); ++i)
        {
            // fwdpp requires that we zero out gamete
            // counts each generation.  Since we're looping
            // over diploids here, now is a good time to
            // handle this task, which saves us from having to
            // do another O(N1+N2) loop:
            pop.gametes[pop.haploids[i]].n = 0;
            if (deme_labels[i] == 0)
                {
                    rv.parents1.push_back(i);
                    w1.push_back( ff(pop.gametes[pop.haploids[i]], pop.mutations, deme_labels[i]) );
                }
            else
                {
                    rv.parents2.push_back(i);
                    w2.push_back( ff(pop.gametes[pop.haploids[i]], pop.mutations, deme_labels[i]) );
                }
        }

    // For lookup tables that RANDOMLY sample parents, i.e. for recombination
    std::vector<double> neut1(rv.parents1.size(), 1.0) ;
    std::vector<double> neut2(rv.parents2.size(), 1.0) ;
    // Set up our lookup tables:
    rv.seln1.reset(gsl_ran_discrete_preproc(rv.parents1.size(), w1.data()));
    rv.seln2.reset(gsl_ran_discrete_preproc(rv.parents2.size(), w2.data()));
    rv.rand1.reset(gsl_ran_discrete_preproc(rv.parents1.size(), neut1.data() ));
    rv.rand2.reset(gsl_ran_discrete_preproc(rv.parents2.size(), neut2.data() ));
    return rv;
};

