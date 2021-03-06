/*! \file diploid.hh
  \brief Main header for programming using this library
  \warning Do not try to include individual headers a la carte. Life will get
  confusing for you. Just include this file.

  \code
  #include <fwdpp/diploid.hh>
  \endcode
 */
#ifndef __HAPLOID_HH__
#define __HAPLOID_HH__

// namespace std
#include <vector>
#include <list>
#include <cassert>
#include <iterator>
#include <ctime>
#include <cmath>

// gsl
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//
// headers from fwdpp
//
#include <fwdpp/type_traits.hpp>
#include <fwdpp/debug.hpp>
#include <fwdpp/forward_types.hpp>
    // ^ defines mutation_base, mutation, gamete_base, and gamete
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/insertion_policies.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <fwdpp/util.hpp>
#include <fwdpp/general_rec_variation.hpp>
#include <fwdpp/interlocus_recombination.hpp>

// personal haploid headers
#include "fitness_models_haploid.hpp"
//#include "sample_haploid_singlepop.hpp"
#include "sample_haploid_singlepop.cpp"
#include "sample_haploid_structure.cpp"
#include "homologous_recpol.hpp"
#include "migration_functions.cpp"

#endif

/*! \namespace fwdpp
  \brief The primary namespace defined by this library.
 */

/*! \namespace fwdpp::fwdpp_internal
  \brief Nested namespace for nuts and bolts of certain library functions
*/

/*! \namespace fwdpp::traits
  \brief Nested namespace type traits
*/

/*! \namespace fwdpp::traits::internal
  \brief Nested namespace implementation details of type traits
*/

/*! \namespace fwdpp::tags
  \brief Nested namespace for dispatch tags for template functions.
*/

/*! @defgroup sugar Syntactic sugar layer
  \brief Syntactic sugar for easier development of simulations

  See @ref md_md_sugar for a full description of the features that fwdpp's
  sugar layer provides.
 */

/*! @defgroup mlocus Multi-locus/region simulations
 * \brief Functions related to modeling multi-locus/region simulations
 */

/*! \namespace fwdpp::sugar
  \brief Nested namespace for sugar layer.

  This namespace provides the implementation details for @ref sugar.

  See @ref md_md_sugar for a full description of the features that fwdpp's
  sugar layer provides.
 */
