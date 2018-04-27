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
//#ifdef HAVE_LIBSEQUENCE
#include <Sequence/SimData.hpp>
//#endif
#include <numeric>
#include <functional>
#include <cassert>
#include <fwdpp/sugar/popgenmut.hpp>
#define SINGLEPOP_SIM
// the type of mutation
using mtype = fwdpp::popgenmut;
#include <gsl/gsl_randist.h>
#include "common_ind.hpp"
#include "haploid.hh"
using namespace fwdpp;

struct parent_lookup_tables
// Our object for choosing parents each generation
{
    // These return indexes of parents from demes 1 and 2,
    // resp, chosen in O(1) time proportional to
    // relative fitness within each deme
    fwdpp_internal::gsl_ran_discrete_t_ptr lookup1, lookup2;
    // These vectors map indexes returned from sampling
    // lookup1 and lookup2 to diploids in the population
    // object.
    std::vector<std::size_t> parents1, parents2;
};

template <typename fitness_fxn>
parent_lookup_tables
migrate_and_calc_fitness(const gsl_rng *r, singlepop_t &pop,
                         const fitness_fxn &wfxn, const uint_t N1,
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

    // Fill a vector of N1 zeros and N2 ones:
    std::vector<uint_t> deme_labels(N1, 0);
    deme_labels.resize(N1 + N2, 1);
    assert(deme_labels.size() == get_sum(pop.gametes));

    // A lookup table to facilitate
    // sampling migrants w/o replacement
    std::unordered_set<uint_t> migrants;

    // Who is migrating 1 -> 2?
    for (unsigned i = 0; i < nmig12; ++i)
        {
            // Sample an individual w/o replacement from deme 1:
            auto mig = static_cast<uint_t>(gsl_ran_flat(r, 0, N1));
            while (migrants.find(mig) != migrants.end())
                {
                    mig = static_cast<uint_t>(gsl_ran_flat(r, 0, N1));
                }
            // Flip the deme label:
            deme_labels[i] = !deme_labels[i];

            // Prevent sampling this individual again:
            migrants.insert(mig);
        }

    // Exact same logic for migrants 2 -> 1
    for (unsigned i = 0; i < nmig21; ++i)
        {
            // Mind your ranges for generating indexes:
            auto mig = static_cast<uint_t>(gsl_ran_flat(r, N1, N1 + N2));
            while (migrants.find(mig) != migrants.end())
                {
                    mig = static_cast<uint_t>(gsl_ran_flat(r, N1, N1 + N2));
                }
            deme_labels[i] = !deme_labels[i];
            migrants.insert(mig);
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
            pop.gametes[i].n = 0;
            if (deme_labels[i] == 0)
                {
                    rv.parents1.push_back(i);
                    w1.push_back(
                        wfxn(pop.gametes[i], pop.mutations));
                }
            else
                {
                    rv.parents2.push_back(i);
                    w2.push_back(
                        wfxn(pop.gametes[i], pop.mutations));
                }
        }

    // Set up our lookup tables:
    rv.lookup1.reset(gsl_ran_discrete_preproc(rv.parents1.size(), w1.data()));
    rv.lookup2.reset(gsl_ran_discrete_preproc(rv.parents2.size(), w2.data()));
    return rv;
};

template <typename fitness_fxn, typename rec_fxn, typename mut_fxn>
void
evolve_two_demes(const gsl_rng *r, singlepop_t &pop, const uint_t N1,
                 const uint_t N2, const double m12, const double m21,
                 const double mu, const fitness_fxn &wfxn,
                 const rec_fxn &recfxn, const mut_fxn &mutfxn)
{
    // Handle mutation/gamete "recycling":
    auto mut_recycling_bin = fwdpp_internal::make_mut_queue(pop.mcounts);
    auto gam_recycling_bin = fwdpp_internal::make_gamete_queue(pop.gametes);

    // Migration and build lookup tables:
    auto sampledemes = migrate_and_calc_fitness(r, pop, wfxn, N1, N2, m12, m21);

#ifndef NDEBUG
    for (const auto &g : pop.gametes)
        assert(!g.n);
#endif

    // Copy parents
    //const auto parents(pop.diploids);

    // Fill in the next generation!
    // We generate the offspring for deme 1 first,
    // and then for deme 2
    for (uint_t i = 0; i < N1 + N2; ++i)
        {
            std::size_t p1 = std::numeric_limits<std::size_t>::max();
            std::size_t p2 = std::numeric_limits<std::size_t>::max();
            if (i < N1) // pick parents from pop 1
                {
                    p1 = sampledemes.parents1[gsl_ran_discrete(r, sampledemes.lookup1.get())];
                    // sample random individual as donor
                    p2 = sampledemes.parents1[gsl_rng_uniform_int(r, sampledemes.parents1.size())];
                }
            else // pick parents from pop 2
                {
                    p1 = sampledemes.parents2[gsl_ran_discrete(r, sampledemes.lookup2.get())];
                    p2 = sampledemes.parents2[gsl_rng_uniform_int(r, sampledemes.parents2.size())];
                }
            assert(p1 < pop.gametes.size());
            assert(p2 < pop.gametes.size());

            /*
              These are the gametes from each parent.
            */
 
            mutate_recombine_update(r, pop.gametes, pop.mutations,
                                    std::make_tuple(p1, p2),
                                    recfxn, mutfxn, mu, gam_recycling_bin,
                                    mut_recycling_bin,
                                    pop.neutral, pop.selected);
        }
    assert(check_sum(pop.gametes, (N1 + N2)));
#ifndef NDEBUG
    /*
    for (const auto &dip : pop.diploids)
        {
            assert(pop.gametes[dip.first].n > 0);
            assert(pop.gametes[dip.first].n <= 2 * (N1 + N2));
            assert(pop.gametes[dip.second].n > 0);
            assert(pop.gametes[dip.second].n <= 2 * (N1 + N2));
        }
     */
#endif

    // Update mutation counts
    fwdpp_internal::process_gametes(pop.gametes, pop.mutations, pop.mcounts);

    assert(pop.mcounts.size() == pop.mutations.size());
#ifndef NDEBUG
    for (const auto &mc : pop.mcounts)
        {
            assert(mc <= 2 * (N1 + N2));
        }
#endif
    for (const auto &g : gametes)
    {
        // check if mutation count data sane for extant gametes
        if(g.n){
            assert(gamete_data_sane(g, mutations, mcounts)) ;
        }
    }
    
    // Prune fixations from gametes
    fwdpp_internal::gamete_cleaner(pop.gametes, pop.mutations, pop.mcounts,
                                   (N1 + N2), std::true_type());
}












/*
 
    EVERYTHING BELOW HERE HAS NOT BEEN HAPLOID-IZED
 */

















int
main(int argc, char **argv)
{
    if (argc != 12)
        {
            std::cerr
                << "Too few arguments.\n"
                << "Usage: juvenile_migration N1 N1 m12 m21 theta_neutral "
                   "theta_deleterious rho s h ngens "
                   "seed\n";
            exit(0);
        }
    int argument = 1;
    const unsigned N1 = atoi(argv[argument++]);
    const unsigned N2 = atoi(argv[argument++]);
    const double m12 = atof(argv[argument++]);
    const double m21 = atof(argv[argument++]);
    const double theta_neutral = atof(argv[argument++]);
    const double theta_del = atof(argv[argument++]);
    const double rho = atof(argv[argument++]);
    const double s = atof(argv[argument++]);
    const double h = atof(argv[argument++]);
    const unsigned ngens = atoi(argv[argument++]);
    const unsigned nreps = atoi(argv[argument++]); // Number of replicates to simulate
    const unsigned seed = atoi(argv[argument++]);

    const unsigned N = N1 + N2; // Total metapop size, for convenience
    const double mu_neutral = theta_neutral / double(4 * N);
    const double mu_del = theta_del / double(4 * N);
    const double littler = rho / double(4 * N);

    std::copy(argv, argv + argc,
              std::ostream_iterator<char *>(std::cerr, " "));
    std::cerr << '\n';

    GSLrng r(seed);

    // recombination map is uniform[0,1)
    const auto rec
        = fwdpp::recbinder(fwdpp::poisson_xover(littler, 0., 1.), r.get());

    const double pselected = mu_del / (mu_del + mu_neutral);

    while(nreps--){
        auto wfxn = fwdpp::multiplicative_diploid(1.);
        
        singlepop_t pop(N);
        pop.mutations.reserve(
            size_t(std::ceil(std::log(2 * N) * (theta_neutral + theta_del)
                             + 0.667 * (theta_neutral + theta_del))));
        unsigned generation = 0;
        const auto mmodel = [&pop, &r, &generation, s, h,
                             pselected](std::queue<std::size_t> &recbin,
                                        singlepop_t::mcont_t &mutations) {
            return fwdpp::infsites_popgenmut(
                recbin, mutations, r.get(), pop.mut_lookup, generation, pselected,
                [&r]() { return gsl_rng_uniform(r.get()); }, [s]() { return s; },
                [h]() { return h; });
        };

        double wbar = 1;
        for (generation = 0; generation < ngens; ++generation)
            {
                assert(fwdpp::check_sum(pop.gametes, 2 * (N1 + N2)));

                // Call our fancy new evolve function
                evolve_two_demes(r.get(), pop, N1, N2, m12, m21,
                                 mu_neutral + mu_del, wfxn, rec, mmodel);
                wbar = fwdpp::sample_diploid(
                    r.get(), pop.gametes, pop.diploids, pop.mutations, pop.mcounts,
                    N, mu_neutral + mu_del, mmodel,
                    // The function to generation recombination positions:
                    rec, wfxn, pop.neutral, pop.selected);
                fwdpp::update_mutations(pop.mutations, pop.fixations,
                                        pop.fixation_times, pop.mut_lookup,
                                        pop.mcounts, generation, 2 * N);
                assert(fwdpp::check_sum(pop.gametes, 2 * N));
            }
    }
}
