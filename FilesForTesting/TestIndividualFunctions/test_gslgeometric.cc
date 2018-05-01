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
//#ifdef HAVE_LIBSEQUENCE
#include <Sequence/SimData.hpp>
//#endif
#include <fwdpp/diploid.hh>
#include <fwdpp/recbinder.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
// typedef mutation_with_age mtype;
using mtype = fwdpp::popgenmut;
#define SINGLEPOP_SIM
#include <common_ind.hpp>

int
main(int argc, char **argv)
{
    
	int argument = 1;
	const unsigned seed = unsigned(atoi(argv[argument++])); // Random number seed

	std::cout << gsl_ran_geometric(r, 0.002) << "\n" ;
	gsl_rng_free(r) ;

    return 0;
}
