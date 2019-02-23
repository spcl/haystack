#include <cassert>

#include "barvinok/isl.h"
#include "isl-helpers.h"
#include <isl/set.h>
#include <limits>

// safe conversion to integer value
long isl::get_value(isl::val Value) {
  isl::val Denominator = isl::manage(isl_val_get_den_val(Value.get()));
  assert(Denominator.is_one());
  return Value.get_num_si();
}

long isl::computeMinimum(isl::set Domain, int Dimension) {
  // try to find constant bounds
  long Minimum = std::numeric_limits<long>::max();
  auto findMinimum = [&](isl::set Domain, isl::aff Value) {
    if (Value.is_cst()) {
      Minimum = std::min(Minimum, isl::get_value(Value.get_constant_val()));
    } else {
      Minimum = std::numeric_limits<long>::min();
    }
    return isl::stat::ok();
  };
  Domain.dim_min(Dimension).foreach_piece(findMinimum);
  return Minimum;
}

long isl::computeMaximum(isl::set Domain, int Dimension) {
  // try to find constant bounds
  long Maximum = std::numeric_limits<long>::min();
  auto findMaximum = [&](isl::set Domain, isl::aff Value) {
    if (Value.is_cst()) {
      Maximum = std::max(Maximum, isl::get_value(Value.get_constant_val()));
    } else {
      Maximum = std::numeric_limits<long>::max();
    }
    return isl::stat::ok();
  };
  Domain.dim_max(Dimension).foreach_piece(findMaximum);
  return Maximum;
}

std::pair<isl::aff, isl::val> isl::get_stride_info(isl::set Set, int Dimension) {
  isl_stride_info *Info = isl_set_get_stride_info(Set.get(), Dimension);
  isl::val Stride = isl::manage(isl_stride_info_get_stride(Info));
  isl::aff Offset = isl::manage(isl_stride_info_get_offset(Info));
  isl_stride_info_free(Info);
  return std::make_pair(Offset, Stride);
}

long isl::cardinality(isl::set Set) { return get_value(isl::manage(isl_set_card(isl::set(Set).release())).max()); }
