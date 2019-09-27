/*!
  Simulate a structured, finite Wright-Fisher population with mutation,
  recombination, and selection. Migration only occurs after ngenMig
  generations

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
    if (argc != 16)
        {
            std::cerr << "Too few arguments\n"
                      << "Usage: haploid_ind N1 N2 m12 m21 theta_neutral theta_sel rho fragsize meantrlen s "
                         "ngens ngenMig samplesize nreps seed\n";
            exit(0);
        }
    int argument = 1;
    const unsigned N1 = unsigned(atoi(argv[argument++])); // Number of haploids in N1
    const unsigned N2 = unsigned(atoi(argv[argument++])); // Number of haploids in N2
    const double m12 = atof(argv[argument++]); // migration deme 1->2
    const double m21 = atof(argv[argument++]); // migration deme 2->1
    const double theta_neutral = atof(argv[argument++]); // 2*n*u PER REGION
    const double theta_sel = atof(argv[argument++]); // 2*n*u PER REGION
    const double rho = atof(argv[argument++]);   // 2*n*r PER REGION
    const double fragsize = (atof(argv[argument++])); // Size of DNA frag, for homologous rec
    const double meantrlen = (atof(argv[argument++])); // mean tract length of DNA transferred
    const double s = atof(argv[argument++]);
    const unsigned ngens = unsigned(atoi(argv[argument++])); // Number of generations to simulate
    const unsigned samplesize1 = unsigned(atoi(argv[argument++])); // Sample size to draw from the population
    int nreps = atoi(argv[argument++]); // Number of replicates to simulate
    const unsigned ngenMig = unsigned(atoi(argv[argument++])); // Number of generations to simulate
    const unsigned seed = unsigned(atoi(argv[argument++])); // Random number seed
    
    const unsigned N = N1 + N2; // Total metapop size, for convenience
    const double mu_neutral = theta_neutral / double(2 * N); // per-gamete mutation rate
    const double mu_sel = theta_sel / double(2 * N); // per-gamete mutation rate
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

    
    // homologous recombination map on interval 0,1
    const auto rec = fwdpp::recbinder(homologous_rec(littler, 0., 1., meantrlen, fragsize), r.get());
    const double pselected = mu_sel / (mu_sel + mu_neutral);
    
    // File to print time taken
    std::string filename = "time.txt" ;
    std::ofstream out ;
    out.open(filename.c_str()) ;
    out.precision(10) ;
   // File to print main output, N1 pop
    std::string filename2 = "results.txt" ;
    std::ofstream results ;
    results.open(filename2.c_str()) ;
    results.precision(10) ;
    // File to print main output, N2 pop
    std::string filename3 = "results2.txt" ;
    std::ofstream results2 ;
    results2.open(filename3.c_str()) ;
    results2.precision(10) ;
    // File to print selected positions
    std::string filename4 = "selected_pos.txt" ;
    std::ofstream selpos ;
    selpos.open(filename4.c_str()) ;
    selpos.precision(10) ;
    // File to print neutral positions
    std::string filename5 = "neutral_pos.txt" ;
    std::ofstream neutpos ;
    neutpos.open(filename5.c_str()) ;
    neutpos.precision(10) ;
    
    // Write the command line to output file
    std::copy(argv, argv + argc, std::ostream_iterator<char *>(results, " "));
    results << '\n';
    std::copy(argv, argv + argc, std::ostream_iterator<char *>(results2, " "));
    results2 << '\n';
    
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
            const auto mmodel
                    = [&pop, &r, &generation, s, pselected](std::queue<std::size_t> &recbin,
                                                            singlepop_t::mcont_t &mutations) {
                    return fwdpp::infsites_popgenmut(
                         recbin, mutations, r.get(), pop.mut_lookup, generation,
                         pselected, [&r]() { return gsl_rng_uniform(r.get()); },
                         [s]() { return s; }, []() { return 0.0; });
            };
            
            for (generation = 0; generation < ngens; ++generation)
                {
                    // Iterate the population through 1 generation
                    sample_haploid_struct_delayed_recmig(
                        r.get(),
                        pop,
                        N1,     // Popsize deme1
                        N2,     // Popsize deme2
                        m12,    // mig rate deme1->deme2
                        m21,    // mig rate deme2->deme1
                        (mu_sel + mu_neutral), // mutation rate per gamete
                        /*
                          The mutation model will be applied
                          by
                          sample_diploid in order to add mutations to gametes
                          each generation.
                        */
                        mmodel,
                        // The function to generation recombination positions:
                        rec,
                        // recombination between populations only occers after ngenMig generations
                        generation,
                        ngenMig,
                        /*
                        Fitness function, can only pass pointers to functions
                         or function objects
                        */
                        pop_neg_seln_multiplicative());
                        // 2 more args in template defn but they have defaults
                    fwdpp::update_mutations(pop.mutations, pop.fixations,
                                            pop.fixation_times, pop.mut_lookup,
                                            pop.mcounts, generation, N);
                    assert(fwdpp::check_sum(pop.gametes, N));
                }
            // group haploids into diploids, in order to use diploid printing function
            std::vector< std::pair<std::size_t, std::size_t>> pseudodips1 ;
            std::vector< std::pair<std::size_t, std::size_t>> pseudodips2 ;
            // collect sample from entire metapopulation
            //pseudodips1 = group_haps_into_dips(r.get(), pop.gametes) ;
            // collect sample from N1
            pseudodips1 = group_N1haps_into_dips(r.get(), pop.gametes, pop.haploids, N1) ;
            pseudodips2 = group_N2haps_into_dips(r.get(), pop.gametes, pop.haploids, N1, N2) ;
                        
            std::vector<std::pair<double, std::string>> mslike1
                    = fwdpp::ms_sample(r.get(), pop.mutations, pop.gametes,
                                       pseudodips1, samplesize1, true);
            std::vector<std::pair<double, std::string>> mslike2
                    = fwdpp::ms_sample(r.get(), pop.mutations, pop.gametes,
                                       pseudodips2, samplesize1, true);
            // Write the sample date a to libsequence's Sequence::SimData and
            // print to file
            Sequence::SimData sdata1;
            if (!mslike1.empty()){
                sdata1.assign(mslike1.begin(), mslike1.end());
                results << sdata1 << '\n';
            }else{
                results << "//\nsegsites: 0\n";
            }
  
            Sequence::SimData sdata2;
            if (!mslike2.empty()){
                sdata2.assign(mslike2.begin(), mslike2.end());
                results2 << sdata2 << '\n';
            }else{
                results2 << "//\nsegsites: 0\n";
            }
            //print file of positions under selection
            selpos << "//\n" ;
            selpos << "SELECTED POSITIONS (pos:gen:freq)\n" ;
            neutpos << "//\n" ;
            neutpos << "NEUTRAL POSITIONS (pos:gen:freq)\n" ;
            for(int i=0; i<pop.mutations.size(); i++ ){
                if(pop.mutations[i].s){
                    selpos << pop.mutations[i].pos << ":" << pop.mutations[i].g << ":" << pop.mcounts[i] <<  "\n" ;
                }else{
                    neutpos << pop.mutations[i].pos << ":" << pop.mutations[i].g << ":" << pop.mcounts[i] <<  "\n" ;
                }
            }
        }
    end = time(0);
    int TimeTaken = int(difftime(end,begin)) ;
    out << "hours:min:sec " << TimeTaken/3600 << ":" << (TimeTaken%3600)/60 << ":" << TimeTaken%60 << "\n" ;
    return 0;
}
