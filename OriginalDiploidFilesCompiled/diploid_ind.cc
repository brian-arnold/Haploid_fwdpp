/*! \include diploid_ind.cc
  Simulate a single, finite Wright-Fisher diploid population with mutation,
  recombination, and no selection.

  This program illustrates many features of fwdpp:
  1.  Custom mutation classes
  2.  Implementing a mutation model (infinitely-many sites)
  3.  Iterating a population through its life cycle
  4.  Outputting a sample in "ms" format
*/
#include <fstream>
#include <iostream>
#include <type_traits>
#include <vector>
#ifdef HAVE_LIBSEQUENCE
#include <Sequence/SimData.hpp>
#endif
#include <fwdpp/diploid.hh>
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
    const unsigned ngens = unsigned(
        atoi(argv[argument++])); // Number of generations to simulate
    const unsigned samplesize1 = unsigned(
        atoi(argv[argument++])); // Sample size to draw from the population
    int nreps = atoi(argv[argument++]); // Number of replicates to simulate
    const unsigned seed
        = unsigned(atoi(argv[argument++])); // Random number seed

    const double mu = theta / double(4 * N); // per-gamete mutation rate

    /*
      littler r is the recombination rate per region per generation.
    */
    const double littler = rho / double(4 * N);

    // Write the command line to stderr
    std::copy(argv, argv + argc,
              std::ostream_iterator<char *>(std::cerr, " "));
    std::cerr << '\n';

    // Initiate random number generation system (via fwdpp/sugar/GSLrng_t.hpp)
    GSLrng r(seed);

    unsigned twoN = 2 * N;

    // recombination map is uniform[0,1)
    const auto rec
        = fwdpp::recbinder(fwdpp::poisson_xover(littler, 0., 1.), r.get());
    std::cout << "starting" << "\n" ;
	
	std::string filename = "time.txt" ;
    std::ofstream out ;
    out.open(filename.c_str()) ;
    // START CLOCK
    time_t begin, end;
    begin = time(0);

    while (nreps--)
        {
            singlepop_t pop(N);
            pop.mutations.reserve(
                size_t(std::ceil(std::log(2 * N) * theta + 0.667 * theta)));
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
                    wbar = fwdpp::sample_diploid(
                        r.get(),
                        pop.gametes,   // non-const reference to gametes
                        pop.diploids,  // non-const reference to diploids
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
                          Fitness is multiplicative over sites.

                          The fitness model takes two gametes and a
                          vector of mutations as arguments.
                          The gametes are passed to this function by
                          fwdpp::sample_diploid, and the _1 and _2 are
                          placeholders for
                          those gametes (see documentation for boost/bind.hpp
                          for details).
                          The mutation container is passed in as _3.
                          The 2. means that fitnesses are 1, 1+sh, and 1+2s for
                          genotypes
                          AA, Aa, and aa, respectively, where a is a mutation
                          with
                          selection coefficient s and dominance h, and the
                          fitness of
                          the diploid is the product of fitness over sites

                          There is no selection in this simulation, but this
                          function is called anyways to illustrate it as
                          multiplicative
                          models are very common in population genetics
                        */
                        fwdpp::multiplicative_diploid(), pop.neutral,
                        pop.selected);
                    fwdpp::update_mutations(pop.mutations, pop.fixations,
                                            pop.fixation_times, pop.mut_lookup,
                                            pop.mcounts, generation, twoN);
                    assert(fwdpp::check_sum(pop.gametes, twoN));
                }
		std::cout << "made it!" << "\n" ;
            // Take a sample of size samplesize1 from the population
            std::vector<std::pair<double, std::string>> mslike
                = fwdpp::ms_sample(r.get(), pop.mutations, pop.gametes,
                                   pop.diploids, samplesize1, true);

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
        end = time(0);
	int TimeTaken = int(difftime(end,begin)) ;
	out << "hours:min:sec " << TimeTaken/3600 << ":" << (TimeTaken%3600)/60 << ":" << TimeTaken%60 << "\n" ;
    return 0;
}
