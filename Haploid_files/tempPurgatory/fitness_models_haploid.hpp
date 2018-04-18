#ifndef _FITNESS_MODELS_HAPLOID_HPP_
#define _FITNESS_MODELS_HAPLOID_HPP_

#include <vector>
#include <algorithm>
#include <fwdpp/forward_types.hpp>
#include <numeric>
#include <functional>
#include <fitness_models_haploid.cpp>

using namespace fwdpp ;

template <typename mtype>
double multiplicative_haploid(const gamete & g, const std::vector<mtype> & mutations) ;


#endif


