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
#define SINGLEPOP_SIM
#include "common_ind.hpp"

int
main(int argc, char **argv)
{
    if (argc != 12)
        {
            std::cerr << "Too few arguments\n"
                      << "Usage: diploid_ind N theta_neutral rho fragsize meantrlen ngens samplesize "
                         "nreps seed\n";
            exit(0);
        }
    int argument = 1;
    const unsigned N = unsigned(atoi(argv[argument++])); // Number of haploids
    const double theta_neutral = atof(argv[argument++]); // 2*n*u PER REGION
    const double theta_sel = atof(argv[argument++]); // 2*n*u PER REGION
    const double rho = atof(argv[argument++]);   // 2*n*r PER REGION
    const double fragsize = (atof(argv[argument++])); // Size of DNA frag, for homologous rec
    const double meantrlen = (atof(argv[argument++])); // mean tract length of DNA transferred
    const double s = atof(argv[argument++]);
    const unsigned ngens = unsigned(atoi(argv[argument++])); // Number of generations to simulate
    const unsigned samplesize1 = unsigned(atoi(argv[argument++])); // Sample size to draw from the population
    int nreps = atoi(argv[argument++]); // Number of replicates to simulate
    const unsigned seed = unsigned(atoi(argv[argument++])); // Random number seed
    
    const double mu_neutral = theta_neutral / double(2 * N); // per-gamete mutation rate
    const double mu_sel = theta_sel / double(2 * N); // per-gamete mutation rate
    const double littler = rho / double(2 * N);// per-gamete reconbination rate

    // Write the command line to stderr
    std::copy(argv, argv + argc, std::ostream_iterator<char *>(std::cerr, " "));
    std::cerr << '\n';

    // Initiate random number generation system (via fwdpp/sugar/GSLrng_t.hpp)
    GSLrng r(seed);

    unsigned twoN = 2 * N;
    unsigned halfN = N / 2 ;

    
    
    
    //double geop = (static_cast<double>(meantrlen))/(static_cast<double>(fragsize)) ;
    // recombination map is uniform[0,1)
    //homologous_rec hrec(littler, 0., 1., meantrlen, fragsize) ;
    //hrec(r.get()) ;
    
    const auto rec = fwdpp::recbinder(homologous_rec(littler, 0., 1., meantrlen, fragsize), r.get());
    const double pselected = mu_sel / (mu_sel + mu_neutral);
    /*
    for(int i=0; i< 10; i++){
        std::vector<double> tmp = rec() ;
        for(int j=0; j< tmp.size(); j++){
            std::cout << tmp[j] << "\t" ;
        }
        std::cout << "\n" ;
    }
     */
    // File to print time taken
    std::string filename1 = "time.txt" ;
    std::ofstream out ;
    out.open(filename1.c_str()) ;
    // File to print main output
    std::string filename2 = "results.txt" ;
    std::ofstream results ;
    results.open(filename2.c_str()) ;
    // File to print selected positions
    std::string filename3 = "selected_pos.txt" ;
    std::ofstream selpos ;
    selpos.open(filename3.c_str()) ;
 
    // Write the command line to output file
    std::copy(argv, argv + argc, std::ostream_iterator<char *>(results, " "));
    results << '\n';
    
    // START CLOCK
    time_t begin, end;
    begin = time(0);
    
    // initialize nfds fitness function with pop size and equilibrium freq
    nfds nfds_func(N, 0.5) ;
    
    
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



            pop.mutations.reserve(size_t(std::ceil(std::log(2 * N) * (theta_neutral+theta_sel) + 0.667 * (theta_neutral+theta_sel))));
            unsigned generation = 0;
            double wbar;
            /*
            // MUTATION MODEL, NEUTRAL
            const auto mmodel =
                [&pop, &r, &generation](std::queue<std::size_t> &recbin,
                                        singlepop_t::mcont_t &mutations) {
                    return fwdpp::infsites_popgenmut(
                        recbin, mutations, r.get(), pop.mut_lookup, generation,
                        0.0, [&r]() { return gsl_rng_uniform(r.get()); },
                        []() { return 0.0; }, []() { return 0.0; });
                };
             */
            // MUTATION MODEL, SEL
            const auto mmodel = [&pop, &r, &generation, s, pselected](std::queue<std::size_t> &recbin,
                                singlepop_t::mcont_t &mutations) {
                                    return fwdpp::infsites_popgenmut(
                                         recbin, mutations, r.get(), pop.mut_lookup, generation,
                                         pselected, [&r]() { return gsl_rng_uniform(r.get()); },
                                         [s]() { return s; }, []() { return 0.0; });
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
                        (mu_neutral + mu_sel), // mutation rate per gamete
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
                        multiplicative_negseln_haploid(),
                        pop.neutral,
                        pop.selected);
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
            std::cout << "\n" ;
            */
            
            std::vector< std::pair<std::size_t, std::size_t>> pseudodips
                = group_haps_into_dips(r.get(), pop.gametes) ;
            
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
            Sequence::SimData sdata;
            
            if (!mslike.empty())
                {
                    sdata.assign(mslike.begin(), mslike.end());
                    results << sdata << '\n';
                }
            else
                {
                    results << "//\nsegsites: 0\n";
                }
  
            selpos << "//\n" ;
            selpos << "SELECTED POSITIONS (pos:gen:freq)\n" ;
            for(int i=0; i<pop.mutations.size(); i++ ){
                if(pop.mutations[i].s){
                    selpos << pop.mutations[i].pos << ":" << pop.mutations[i].g << ":" << pop.mcounts[i] <<  "\n" ;
                }
            }

            //for(unsigned i = 0; i < pop.mcounts.size(); i++){
            //    std::cout << i << " " << pop.mcounts[i] << "\n" ;
            //}
        }
    end = time(0);
    int TimeTaken = int(difftime(end,begin)) ;
    out << "hours:min:sec " << TimeTaken/3600 << ":" << (TimeTaken%3600)/60 << ":" << TimeTaken%60 << "\n" ;
    return 0;
}
