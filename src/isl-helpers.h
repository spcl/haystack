/*
* Copyright (c) 2019, ETH Zurich
*/

#ifndef _ISL_HELPERS_H_
#define _ISL_HELPERS_H_

#include <isl/isl-noexceptions.h>

namespace isl {

// val
long get_value(isl::val Value);

// set
long computeMinimum(isl::set Domain, int Dimension);
long computeMaximum(isl::set Domain, int Dimension);

std::pair<isl::aff, isl::val> get_stride_info(isl::set Set, int Dimension);

long cardinality(isl::set Set);

std::string printExpression(isl::aff Aff);

} // namespace isl

#endif