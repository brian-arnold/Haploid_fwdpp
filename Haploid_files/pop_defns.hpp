#ifndef __POP_DEFNS_HPP__
#define __POP_DEFNS_HPP__

#include <utility>
#include <vector>
#include <unordered_set>
#include <fwdpp/fwd_functional.hpp>
#include "gametepop.hpp"
#include "gametepop_structured.hpp"

using namespace fwdpp ;

/*!
 Many important object defined here
*/
template <typename mtype>
using gametepop_obj
    = gametepop<mtype,
                std::vector<mtype>,
                std::vector<gamete>,    // see fwdpp/forward_types.hpp
                std::vector<mtype>,
                std::vector<uint_t>,
                std::unordered_set<double, std::hash<double>, fwdpp::equal_eps>>;

template <typename mtype, typename haploid_t = std::size_t>
using gametepopstruct_obj
    = gametepopstruct<mtype,
                std::vector<mtype>,
                std::vector<gamete>,    // see fwdpp/forward_types.hpp
                std::vector<haploid_t>,
                std::vector<mtype>,
                std::vector<uint_t>,
                std::unordered_set<double, std::hash<double>, fwdpp::equal_eps>>;

#endif
