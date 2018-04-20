/*! \include diploid_ind.cc
  Simulate a single, finite Wright-Fisher diploid population with mutation,
  recombination, and no selection.

  This program illustrates many features of fwdpp:
  1.  Custom mutation classes
  2.  Implementing a mutation model (infinitely-many sites)
  3.  Iterating a population through its life cycle
  4.  Outputting a sample in "ms" format
*/
#include <iostream>
#include <type_traits>
#include <vector>
#ifdef HAVE_LIBSEQUENCE
#include <Sequence/SimData.hpp>
#endif
#include "haploid.hh"
#include <fwdpp/recbinder.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
// typedef mutation_with_age mtype;
using mtype = fwdpp::popgenmut;
#define SINGLEPOP_SIM
#include "common_ind.hpp"

int
main(int argc, char **argv)
{
    if (argc != 8)
        {
            std::cerr << "Too few arguments\n"
                      << "Usage: diploid_ind N theta rho ngens samplesize "
                         "nreps seed\n";
            exit(0);
        }
    int argument = 1;
    const unsigned N = unsigned(atoi(argv[argument++])); // Number of diploids
    const double theta = atof(argv[argument++]); // 4*n*mutation rate.  Note:
                                                 // mutation rate is per
                                                 // REGION, not SITE!!
    const double rho = atof(argv[argument++]);   // 4*n*recombination rate.
                                                 // Note: recombination rate is
                                                 // per REGION, not SITE!!
    const unsigned ngens = unsigned(atoi(argv[argument++])); // Number of generations to simulate
    const unsigned samplesize1 = unsigned(atoi(argv[argument++])); // Sample size to draw from the population
    int nreps = atoi(argv[argument++]); // Number of replicates to simulate
    const unsigned seed = unsigned(atoi(argv[argument++])); // Random number seed
    const double mu = theta / double(4 * N); // per-gamete mutation rate

    /*
      littler r is the recombination rate per region per generation.
    */
    const double littler = rho / double(4 * N);

    // Write the command line to stderr
    std::copy(argv, argv + argc, std::ostream_iterator<char *>(std::cerr, " "));
    std::cerr << '\n';

    // Initiate random number generation system (via fwdpp/sugar/GSLrng_t.hpp)
    GSLrng r(seed);

    unsigned twoN = 2 * N;
    unsigned halfN = N / 2 ;

    // recombination map is uniform[0,1)
    const auto rec = fwdpp::recbinder(fwdpp::poisson_xover(littler, 0., 1.), r.get());
    std::cout << "starting" << "\n" ;
    while (nreps--)
        {
            /* A few hacks to deal with popbase base class
             1.) initialize with N/2 since popbase sets first
             gamete to 2N copies, then manually change back
             to N
             2.) initialize rest of gametes with .n=0
            */
            singlepop_t pop(halfN);
            pop.N = N ;
            std::uint32_t zero = 0 ;
            for(int i=1; i< N; i++){
                pop.gametes.push_back(fwdpp::gamete(zero)) ;
            }



            pop.mutations.reserve(size_t(std::ceil(std::log(2 * N) * theta + 0.667 * theta)));
            unsigned generation = 0;
            double wbar;

            const auto mmodel =
                [&pop, &r, &generation](std::queue<std::size_t> &recbin,
                                        singlepop_t::mcont_t &mutations) {
                    return fwdpp::infsites_popgenmut(
                        recbin, mutations, r.get(), pop.mut_lookup, generation,
                        0.0, [&r]() { return gsl_rng_uniform(r.get()); },
                        []() { return 0.0; }, []() { return 0.0; });
                };
            for (generation = 0; generation < ngens; ++generation)
                {
                    // Iterate the population through 1 generation
                    wbar = sample_haploid(
                        r.get(),
                        pop.gametes,   // non-const reference to gametes
                        pop.mutations, // non-const reference to mutations
                        pop.mcounts,
                        N,  // current pop size, remains constant
                        mu, // mutation rate per gamete
                        /*
                          The mutation model will be applied
                          by
                          sample_diploid in order to add mutations to gametes
                          each generation.
                        */
                        mmodel,
                        // The function to generation recombination positions:
                        rec,
                        /*
                        Fitness function, can only pass pointers to functions
                         or function objects
                        */
                        multiplicative_haploid(),
                        pop.neutral,
                        pop.selected);
                        // 2 more args in template defn but they have defaults
                    fwdpp::update_mutations(pop.mutations, pop.fixations,
                                            pop.fixation_times, pop.mut_lookup,
                                            pop.mcounts, generation, N);
                    assert(fwdpp::check_sum(pop.gametes, N));
                }
		std::cout << "made it!" << "\n" ;
            // Take a sample of size samplesize1 from the population
            // sampleHap of length N, elements = # times gamete sampled
            std::vector<unsigned> sampleHap
                = fwdpp::sample(r.get(), pop.gametes,
                            samplesize1, N);
            
            for(int i=0; i< sampleHap.size() ; i++){
                std::cout << sampleHap[i] << "\t" ; ;
            }
            std::cout << "\n" ;
            /*
            std::vector<std::pair<double, std::string>> mslike
                = fwdpp::ms_sample(r.get(), pop.mutations, pop.gametes,
                                   pop.diploids, samplesize1, true);
            */
// Write the sample date a to libsequence's Sequence::SimData and
// print to screen
#ifdef HAVE_LIBSEQUENCE
            Sequence::SimData sdata;
            if (!mslike.empty())
                {
                    sdata.assign(mslike.begin(), mslike.end());
                    std::cout << sdata << '\n';
                }
            else
                {
                    std::cout << "//\nsegsites: 0\n";
                }
#endif
        }
    return 0;
}
