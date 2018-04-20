#ifndef __GAMETEPOP_DEFN_HPP__
#define __GAMETEPOP_DEFN_HPP__

#include <utility>
#include <vector>
#include <unordered_set>
#include <fwdpp/fwd_functional.hpp>
#include "gametepop.hpp"

using namespace fwdpp ;

/*!
 Many important object defined here
*/
template <typename mtype>
using gametepop_obj
    = gametepop<mtype,
                std::vector<mtype>,
                std::vector<gamete>,
                std::vector<mtype>,
                std::vector<uint_t>,
                std::unordered_set<double, std::hash<double>, fwdpp::equal_eps>>;

#endif
