/*
* Copyright (c) 2019, ETH Zurich
*/

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

long isl::cardinality(isl::set Set) {
  return get_value(isl::manage(isl_set_card(isl::set(Set).release())).max());
}

std::string isl::printExpression(isl::aff Expression) {
  bool isFirst = true;
  std::string Result = "";
  auto printValue = [&](isl::val Value, bool isCoefficient) {
    std::string Result = "";
    // print the sign if needed
    int Sign = Value.sgn();
    if (Sign < 0) {
      Value = Value.neg();
      Result = "-";
    } else if (!isFirst) {
      Result = "+";
    }
    if (!(Value.is_one() && isCoefficient)) {
      Result += Value.to_str();
    }
    return Result;
  };
  auto printAff = [&](isl::aff Aff) {
    isl::val Denominator = Aff.get_denominator_val();
    if (!Denominator.is_one()) {
      Result += "(";
    }
    // convert the parameters
    for (int i = 0; i < Aff.dim(isl::dim::param); ++i) {
      isl::val Coefficient = Aff.get_coefficient_val(isl::dim::param, i);
      if (!Coefficient.is_zero()) {
        Result += printValue(Coefficient.mul(Denominator), true);
        Result += Aff.get_dim_name(isl::dim::param, i);
        isFirst = false;
      }
    }
    // convert the variables
    for (int i = 0; i < Aff.dim(isl::dim::in); ++i) {
      isl::val Coefficient = Aff.get_coefficient_val(isl::dim::in, i);
      if (!Coefficient.is_zero()) {
        Result += printValue(Coefficient.mul(Denominator), true);
        Result += Aff.get_dim_name(isl::dim::in, i);
        isFirst = false;
      }
    }
    // add the constant
    isl::val Constant = Aff.get_constant_val();
    if (!Constant.is_zero() || isFirst) {
      Result += printValue(Constant.mul(Denominator), false);
      isFirst = false;
    }
    if (!Denominator.is_one()) {
      isFirst = true;
      Result += ")/" + printValue(Denominator, false);
      isFirst = false;
    }
  };
  // pint the main expression
  printAff(Expression);
  // add possible divisors
  for (int i = 0; i < Expression.dim(isl::dim::div); ++i) {
    isl::aff Divisor = Expression.get_div(i);
    isl::val Coefficient = Expression.get_coefficient_val(isl::dim::div, i);
    if (!Coefficient.is_zero()) {
      Result += printValue(Coefficient, true);
      Result += "\u230A";
      isFirst = true;
      printAff(Divisor);
      Result += "\u230B";
    }
  }
  return Result;
}
