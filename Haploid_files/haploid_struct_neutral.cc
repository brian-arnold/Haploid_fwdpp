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
#include <fstream>
#include <type_traits>
#include <vector>
//#ifdef HAVE_LIBSEQUENCE
#include <Sequence/SimData.hpp>
//#endif
#include "haploid.hh"
#include <fwdpp/recbinder.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
// typedef mutation_with_age mtype;
using mtype = fwdpp::popgenmut;
#define MULTIPOP_SIM
#include "common_ind.hpp"

int
main(int argc, char **argv)
{
    if (argc != 13)
        {
            std::cerr << "Too few arguments\n"
                      << "Usage: haploid_ind N1 N2 m12 m21 theta rho fragsize meantrlen "
                         "ngens samplesize nreps seed\n";
            exit(0);
        }
    int argument = 1;
    const unsigned N1 = unsigned(atoi(argv[argument++])); // Number of haploids in N1
    const unsigned N2 = unsigned(atoi(argv[argument++])); // Number of haploids in N2
    const double m12 = atof(argv[argument++]); // migration deme 1->2
    const double m21 = atof(argv[argument++]); // migration deme 2->1
    const double theta = atof(argv[argument++]); // 2*n*u PER REGION
    const double rho = atof(argv[argument++]);   // 2*n*r PER REGION
    const double fragsize = (atof(argv[argument++])); // Size of DNA frag, for homologous rec
    const double meantrlen = (atof(argv[argument++])); // mean tract length of DNA transferred
    const unsigned ngens = unsigned(atoi(argv[argument++])); // Number of generations to simulate
    const unsigned samplesize1 = unsigned(atoi(argv[argument++])); // Sample size to draw from the population
    int nreps = atoi(argv[argument++]); // Number of replicates to simulate
    const unsigned seed = unsigned(atoi(argv[argument++])); // Random number seed
    
    const unsigned N = N1 + N2; // Total metapop size, for convenience
    const double mu = theta / double(2 * N); // per-gamete mutation rate
    const double littler = rho / double(2 * N);// per-gamete reconbination rate

    // Write the command line to stderr
    std::copy(argv, argv + argc, std::ostream_iterator<char *>(std::cerr, " "));
    std::cerr << '\n';
    // Write the command line to output file
    std::copy(argv, argv + argc, std::ostream_iterator<char *>(std::cout, " "));
    std::cout << '\n';
    // Initiate random number generation system (via fwdpp/sugar/GSLrng_t.hpp)
    GSLrng r(seed);

    unsigned twoN = 2 * N;
    unsigned halfN = N / 2 ;

    
    
    
    //double geop = (static_cast<double>(meantrlen))/(static_cast<double>(fragsize)) ;
    // recombination map is uniform[0,1)
    //homologous_rec hrec(littler, 0., 1., meantrlen, fragsize) ;
    //hrec(r.get()) ;
    
    const auto rec = fwdpp::recbinder(homologous_rec(littler, 0., 1., meantrlen, fragsize), r.get());
    /*
    for(int i=0; i< 10; i++){
        std::vector<double> tmp = rec() ;
        for(int j=0; j< tmp.size(); j++){
            std::cout << tmp[j] << "\t" ;
        }
        std::cout << "\n" ;
    }
     */
    std::string filename = "time.txt" ;
    std::ofstream out ;
    out.open(filename.c_str()) ;
    // START CLOCK
    time_t begin, end;
    begin = time(0);
    
    while (nreps--)
        {
            /* A few hacks to deal with popbase base class
             1.) initialize with N/2 since popbase sets first
             gamete to 2N copies, then manually change back
             to N
             2.) initialize rest of haploids, since constructed with N/2
            */
            structpop_t pop(halfN);
            pop.N = N ;
            for(int i=0; i< halfN; i++){
                pop.haploids.emplace_back(0) ;
            }

            pop.mutations.reserve(size_t(std::ceil(std::log(2 * N) * theta + 0.667 * theta)));
            unsigned generation = 0;
            double wbar;

            // MUTATION MODEL, NEUTRAL
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
                    sample_haploid_structure(
                        r.get(),
                        pop,
                        N1,     // Popsize deme1
                        N2,     // Popsize deme2
                        m12,    // mig rate deme1->deme2
                        m21,    // mig rate deme2->deme1
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
                        pop_sign_seln_multiplicative());
                        // 2 more args in template defn but they have defaults
                    fwdpp::update_mutations(pop.mutations, pop.fixations,
                                            pop.fixation_times, pop.mut_lookup,
                                            pop.mcounts, generation, N);
                    assert(fwdpp::check_sum(pop.gametes, N));
                }
            // Take a sample of size samplesize1 from the population
            // sampleHap of length N, elements = # times gamete sampled
            /*
            std::vector<unsigned> sampleHap
                = fwdpp::sample(r.get(), pop.gametes,
                            samplesize1, N);
            
            for(int i=0; i< sampleHap.size() ; i++){
                std::cout << sampleHap[i] << "\t" ; ;
            }
             */
            
            std::vector< std::pair<std::size_t, std::size_t>> pseudodips ;
                pseudodips = group_haps_into_dips(r.get(), pop.gametes) ;
            
            /*
            std::vector<std::pair<double, std::string>> mslike
                = fwdpp::ms_sample(r.get(), pop.mutations, pop.gametes,
                                   pop.diploids, samplesize1, true);
             */
            std::vector<std::pair<double, std::string>> mslike
            = fwdpp::ms_sample(r.get(), pop.mutations, pop.gametes,
                               pseudodips, samplesize1, true);
            
// Write the sample date a to libsequence's Sequence::SimData and
// print to screen
//#ifdef HAVE_LIBSEQUENCE
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
//#endif
            //for(unsigned i = 0; i < pop.mcounts.size(); i++){
            //    std::cout << i << " " << pop.mcounts[i] << "\n" ;
            //}
        }
    end = time(0);
    int TimeTaken = int(difftime(end,begin)) ;
    out << "hours:min:sec " << TimeTaken/3600 << ":" << (TimeTaken%3600)/60 << ":" << TimeTaken%60 << "\n" ;
    return 0;
}
