/*
 * Copyright (c) 2019, ETH Zurich
 */

#include <boost/math/common_factor.hpp>

#include <algorithm>
#include <cassert>
#include <limits>
#include <list>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include "Access.h"
#include "Timer.h"

#include "barvinok/isl.h"
#include "isl-helpers.h"
#include "op.h"

void Access::initAccess(std::vector<NamedLong> ParameterValues, isl::set Parameters) {
  // set the parametes
  ParameterValues_ = ParameterValues;
  Parameters_ = Parameters;
  // compute the number of cache levels
  long CacheLevels = MachineModel_.CacheSizes.size();
  // number of cache misses
  Result_.CompulsoryMisses = 0;
  Result_.CapacityMisses = std::vector<long>(CacheLevels, 0);
  Result_.Counted = 0;
  // compute the domain size
  Result_.Total = isl::cardinality(Domain_);
  // clear the stack distances
  Misses_ = 0;
  Affine_.clear();
  NonAffine_.clear();
  Constant_.clear();
}

void Access::countCompulsoryMisses(isl::union_map First) {
  Timer::startTimer("CountCompulsoryMisses");
  // extract access first map
  First = First.intersect_domain(Domain_);
  if (!First.is_empty()) {
    // intersect with the parameter values
    if (!Parameters_.is_null())
      First = First.intersect_params(Parameters_);
    // simplify map
    First = First.coalesce();
    First = First.detect_equalities();
    // count the individual sets
    auto countCompulsory = [&](isl::set Set) {
      Result_.CompulsoryMisses += isl::cardinality(Set);
      return isl::stat::ok();
    };
    First.range().foreach_set(countCompulsory);
  }
  Timer::stopTimer("CountCompulsoryMisses");
}

void Access::computeStackDistances(isl::union_map BetweenMap) {
  Timer::startTimer("ComputeStackDistances");
  // extract access between map
  BetweenMap = BetweenMap.intersect_domain(Domain_);
  if (!BetweenMap.is_empty()) {
    if (!Parameters_.is_null())
      BetweenMap = BetweenMap.intersect_params(Parameters_);
    // simplify between maps
    BetweenMap = BetweenMap.detect_equalities();
    // compute stack distance bounds
    if (ModelOptions_.ComputeBounds) {
      // compute the set of obvious cache misses
      Timer::startTimer("ComputeBounds");
      long Limit = *std::max_element(MachineModel_.CacheSizes.begin(), MachineModel_.CacheSizes.end());
      Limit /= MachineModel_.CacheLineSize;
      isl::union_set Misses = BetweenMap.domain().empty(BetweenMap.domain().get_space());
      isl::union_set Domain = BetweenMap.domain();
      // compute number of maps
      auto countMap = [&](isl::map Map) {
        auto countBasicMap = [&](isl::basic_map BasicMap) {
          // check if the map is already part of the miss set
          if (!Domain.intersect(BasicMap.domain()).is_empty()) {
            // compute the lower bound for the map
            auto Map = isl::map(BasicMap);
            auto Norm = Map.sum(Map.lexmin().neg());
            long UpperBound = isl::cardinality(Norm.range());
            if (UpperBound > Limit) {
              auto Complement = Norm.complement().intersect_domain(Map.domain());
              long LowerBound = isl::cardinality(Complement.range().complement().intersect(Norm.range()));
              if (LowerBound > Limit) {
                Misses = Misses.unite(BasicMap.domain());
                Domain = Domain.subtract(BasicMap.domain());
              }
            }
          }
          return isl::stat::ok();
        };
        Map.foreach_basic_map(countBasicMap);
        return isl::stat::ok();
      };
      BetweenMap.foreach_map(countMap);
      // count the domains with stack distance bound larger than the maximal cache size
      if (!Misses.is_empty()) {
        // count the misses
        auto countMisses = [&](isl::set Set) {
          long Count = isl::cardinality(Set);
          Result_.Counted += Count;
          Misses_ += Count;
          return isl::stat::ok();
        };
        Misses.foreach_set(countMisses);
        // update the between map
        BetweenMap = BetweenMap.subtract_domain(Misses);
      }
      Timer::stopTimer("ComputeBounds");
    }
    // compute the stack distance for the remaining points of the iteration domain
    Timer::startTimer("CountBetweenMap");
    Expression_.clear();
    auto Count = isl::manage(isl_union_map_card(isl::union_map(BetweenMap).release()));
    Count = Count.intersect_domain(BetweenMap.domain());
    auto countMapDomain = [&](isl::set Set) {
      Result_.Counted += isl::cardinality(Set);
      return isl::stat::ok();
    };
    isl::union_map(BetweenMap).domain().foreach_set(countMapDomain);
    Timer::stopTimer("CountBetweenMap");
    Timer::stopTimer("ComputeStackDistances");
    Timer::startTimer("CountCapacityMisses");
    // extract the cache miss expression
    extractStackDistanceExpression(Count);
    // count affine pieces immediately
    storeAffinePieces();
#ifdef EQUALIZATION
    // perform floor term equalization
    applyEqualization();
    // count affine pieces again
    storeAffinePieces();
#endif
#ifdef RASTERIZATION
    // perform floor term rasterization
    applyRasterization();
    // count affine pieces again
    storeAffinePieces();
#endif
    // process all non affine pieces
    for (auto &Piece : Expression_) {
      enumerateNonAffineDimensions(Piece);
    }
    Timer::stopTimer("CountCapacityMisses");
  }
}

void Access::countCapacityMisses() {
  // compute the capacity misses for the machine
  Result_.CapacityMisses = countCapacityMisses(MachineModel_.CacheSizes);
}

std::vector<long> Access::countCapacityMisses(std::vector<long> CacheSizes) {
  Timer::startTimer("CountCapacityMisses");
  if (ModelOptions_.ComputeBounds) {
    // verify the cache sizes do not exceed the maximum cache size of the machine
    // (the bounds were computed for the maximal cache size)
    if (*std::max_element(MachineModel_.CacheSizes.begin(), MachineModel_.CacheSizes.end()) <
        *std::max_element(CacheSizes.begin(), CacheSizes.end())) {
      printf("-> exit(-1) cache size exceeds maximum cache size of the machine\n");
      exit(-1);
    }
  }
  // initializes the limit and result vectors
  std::vector<long> Results;
  std::vector<long> Limits;
  for (auto Size : CacheSizes) {
    Results.push_back(Misses_);
    Limits.push_back(Size / MachineModel_.CacheLineSize);
  }
  // compute the cache misses for the constant domains
  for (auto &Domain : Constant_) {
    for (int i = 0; i < Limits.size(); ++i) {
      Results[i] += Domain.second.gt(isl::val(Domain.second.get_ctx(), Limits[i])) ? Domain.first : 0;
    }
  }
  // compute the cache misses for the affine pieces
  for (auto &Piece : Affine_) {
    auto Misses = countAffineDimensions(Piece, Limits);
    std::transform(Results.begin(), Results.end(), Misses.begin(), Results.begin(), std::plus<long>());
  }
  // compute the cache misses for the non-affine pieces
  for (auto &Piece : NonAffine_) {
    auto Misses = enumerateNonAffinePoints(Piece, Limits);
    std::transform(Results.begin(), Results.end(), Misses.begin(), Results.begin(), std::plus<long>());
  }
  Timer::stopTimer("CountCapacityMisses");
  return Results;
}

void Access::storeAffinePieces() {
  // count and remove all affine pieces
  std::vector<piece> NonAffine;
  for (auto &Piece : Expression_) {
    if (!Piece.Expression.is_null()) {
      long All = getPieceSize(Piece);
      Affine_.push_back(Piece);
    } else {
      NonAffine.push_back(Piece);
    }
  }
  std::swap(Expression_, NonAffine);
}

std::vector<long> Access::enumerateNonAffinePoints(piece Piece, std::vector<long> Limits) const {
  Timer::startTimer("enumerateNonAffinePoints");
  std::vector<long> Results(Limits.size(), 0);
  auto countMisses = [&](isl::point Point) {
    isl::val StackDistance = Piece.Polynomial.eval(Point);
    for (int i = 0; i < Limits.size(); ++i) {
      if (isl::get_value(StackDistance) > Limits[i])
        Results[i]++;
    }
    return isl::stat::ok();
  };
  Piece.Domain.foreach_point(countMisses);
  Timer::stopTimer("enumerateNonAffinePoints");
  return Results;
}

std::vector<long> Access::countAffineDimensions(piece Piece, std::vector<long> Limits) const {
  Timer::startTimer("countAffineDimensions");
  std::vector<long> Results(Limits.size());
  isl::pw_aff LHS = Piece.Expression;
  for (int i = 0; i < Limits.size(); ++i) {
    isl::pw_aff RHS(Piece.Domain, isl::val(Piece.Domain.get_ctx(), Limits[i]));
    isl::set Misses = LHS.gt_set(RHS);
    auto Variable = isl::manage(isl_set_card(Misses.release()));
    assert(isl::get_value(Variable.min()) == isl::get_value(Variable.max()));
    Results[i] = isl::get_value(Variable.max());
  }
  Timer::stopTimer("countAffineDimensions");
  return Results;
}

std::vector<int> Access::findNonAffineDimensions(piece Piece) const {
  // compute the non-affine dimensions and the conflicting dimensions
  std::set<int> NonAffine;
  std::map<int, std::set<int>> Conflicts;
  // first analyze all terms
  auto analyzeTerms = [&](isl::term Term) {
    // compute the parameter and variable exponent vectors
    std::vector<int> Parameters(Term.dim(isl::dim::param), 0);
    std::vector<int> Variables(Term.dim(isl::dim::set), 0);
    // count the parameter exponents
    for (int i = 0; i < Term.dim(isl::dim::param); i++) {
      Parameters[i] += Term.get_exp(isl::dim::param, i);
    }
    // count the variable exponents
    for (int i = 0; i < Term.dim(isl::dim::set); i++) {
      Variables[i] += Term.get_exp(isl::dim::set, i);
    }
    // get the divisor variable
    for (int i = 0; i < Term.dim(isl::dim::div); i++) {
      if (Term.get_exp(isl::dim::div, i) > 0) {
        isl::aff Divisor = Term.get_div(i);
        // extract the parameter and variable information
        for (int j = 0; j < Divisor.dim(isl::dim::param); ++j) {
          isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::param, j);
          if (!Coefficient.is_zero())
            Parameters[j] += Term.get_exp(isl::dim::div, i);
        }
        for (int j = 0; j < Divisor.dim(isl::dim::in); ++j) {
          isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::in, j);
          if (!Coefficient.is_zero())
            Variables[j] += Term.get_exp(isl::dim::div, i);
        }
      }
    }
    // we do not support parameters
    assert(std::accumulate(Parameters.begin(), Parameters.end(), 0) == 0);
    // update the non-affine variables
    std::vector<int> All;
    for (int i = 0; i < Variables.size(); ++i) {
      if (Variables[i] > 1)
        NonAffine.insert(i);
      if (Variables[i] > 0)
        All.push_back(i);
    }
    // update the conflicts
    for (auto Variable : All) {
      for (auto Conflict : All) {
        if (Conflict != Variable)
          Conflicts[Variable].insert(Conflict);
      }
    }
    return isl::stat::ok();
  };
  Piece.Polynomial.foreach_term(analyzeTerms);
  // compute the minimal set of non-affine Dimensions
  std::vector<int> Dimensions;
  auto storeAndRemoveFromConflicts = [&](int Dimension) {
    // remove from conflicts
    Conflicts.erase(Dimension);
    for (auto Iter = Conflicts.begin(); Iter != Conflicts.end();) {
      Iter->second.erase(Dimension);
      if (Iter->second.empty()) {
        Conflicts.erase(Iter++);
      } else {
        ++Iter;
      }
    }
    // store the dimension
    Dimensions.push_back(Dimension);
  };
  // always add all non-affine dimensions
  for (auto Dimension : NonAffine) {
    storeAndRemoveFromConflicts(Dimension);
  }
  // eliminate all conflicts
  while (!Conflicts.empty()) {
    int Maximum = 0;
    int Dimension = 0;
    for (auto Conflict : Conflicts) {
      if (Conflict.second.size() > Maximum) {
        Maximum = Conflict.second.size();
        Dimension = Conflict.first;
      }
    }
    assert(Maximum != 0);
    // store the dimension and
    storeAndRemoveFromConflicts(Dimension);
  }
  // return the sorted dimensions
  std::sort(Dimensions.begin(), Dimensions.end());
  return Dimensions;
}

void Access::enumerateNonAffineDimensions(piece Piece) {
  Timer::startTimer("enumerateNonAffineDimensions");
  // count total number of accesses
  long All = getPieceSize(Piece);
#ifdef ENUMERATE_POINTS
  // count manually if only the memory access dimension is non-affine
  NonAffine_.push_back(Piece);
#else
  // compute the dimensions
  std::vector<int> NonAffine = findNonAffineDimensions(Piece);
  std::vector<int> Affine;
  for (int i = 0; i < Piece.Domain.dim(isl::dim::set); ++i) {
    auto Iter = std::find(NonAffine.begin(), NonAffine.end(), i);
    if (Iter == NonAffine.end())
      Affine.push_back(i);
  }
  // enumerate and count the affine subdomains
  if (Affine.size() <= 1 || All == 0) {
    // count manually if only the memory access dimension is non-affine
    NonAffine_.push_back(Piece);
  } else {
    // compute the enumeration domain
    isl::set Enumeration = Piece.Domain;
    for (auto Iter = Affine.rbegin(); Iter != Affine.rend(); Iter++) {
      Enumeration = Enumeration.project_out(isl::dim::set, *Iter, 1);
    }
    // check if there is enough work for barvinok
    if (All < 16 || All / isl::cardinality(Enumeration) < 16) {
      NonAffine_.push_back(Piece);
    } else {
      // compute the piecewise polynomial
      auto countCacheMisses = [&](isl::basic_set All) {
        isl::set Enumeration = All;
        for (auto Iter = Affine.rbegin(); Iter != Affine.rend(); Iter++) {
          Enumeration = Enumeration.project_out(isl::dim::set, *Iter, 1);
        }
        // enumerate the points of the enumeration domain
        auto countSubdomain = [&](isl::point Point) {
          // compute the updated domain of the piece
          std::map<int, long> Values;
          auto Domain = All;
          for (int i = 0; i < Point.get_space().dim(isl::dim::set); ++i) {
            auto Value = Point.get_coordinate_val(isl::dim::set, i);
            Domain = Domain.fix_val(isl::dim::set, NonAffine[i], Value);
            Values[NonAffine[i]] = isl::get_value(Value);
          }
          auto Expression = extractAffineExpression(Piece.Polynomial, Domain, Values);
          // handle constant expressions efficiently
          if (Expression.is_cst()) {
            long Size = isl::cardinality(Domain);
            auto Constant = Expression.get_constant_val();
            Constant_.push_back(std::make_pair(Size, Constant));
          } else {
            piece Piece;
            Piece.Domain = Domain;
            Piece.Expression = Expression;
            Affine_.push_back(Piece);
          }
          return isl::stat::ok();
        };
        Enumeration.foreach_point(countSubdomain);
        return isl::stat::ok();
      };
      // we evaluate the basic sets separately since otherwise foreach counts some points multiple times
      Piece.Domain = isl::manage(isl_set_make_disjoint(Piece.Domain.release()));
      Piece.Domain.foreach_basic_set(countCacheMisses);
    }
  }
#endif
  Timer::stopTimer("enumerateNonAffineDimensions");
}

void Access::applyEqualization() {
  Timer::startTimer("applyEqualization");
  // convert the terms to to piecewise aff expressions
  std::vector<piece> Expression;
  for (auto Piece : Expression_) {
    // try to eliminate non affine terms
    std::vector<piece> Pieces = {Piece};
    while (!Pieces.empty()) {
      // compute possible candidates and select the one with the highest exponent
      // (doe not split polynoms with a single divisor term)
      auto Candidates = findEqualizationCandidates(Pieces.back());
      std::vector<std::vector<long>> Splits;
      std::vector<int> Indexes;
      for (auto Candidate : Candidates) {
        Splits.push_back(computeSplits(Candidate, Pieces.back()));
      }
      for (int i = 0; i < Candidates.size(); ++i) {
        Indexes.push_back(i);
      }
      std::sort(Indexes.begin(), Indexes.end(), [&](int I1, int I2) { return Splits[I1].size() < Splits[I2].size(); });
      // classify the updates
      std::vector<piece> Done;
      std::vector<piece> Empty;
      std::vector<piece> Improved;
      // try to eliminate all the candidates and stop on success
      for (auto Index : Indexes) {
        std::vector<piece> Updates = equalizeCandidate(Candidates[Index], Splits[Index], Pieces.back());
        // store updates if the exponent improves
        if (!Updates.empty()) {
          for (auto Update : Updates) {
            if (Update.Polynomial.is_null()) {
              Empty.push_back(Update);
            } else if (isPieceAffine(Update)) {
              Update.Expression = extractAffineExpression(Update);
              Done.push_back(Update);
            } else {
              Improved.push_back(Update);
            }
          }
          break;
        }
      }
      // keep the updates or discard everything
      if (Improved.empty() && Done.empty() && Empty.empty()) {
        Expression.push_back(Pieces.back());
        Pieces.pop_back();
      } else {
        Pieces.pop_back();
        Pieces.insert(Pieces.end(), Improved.begin(), Improved.end());
        Expression.insert(Expression.end(), Done.begin(), Done.end());
      }
    }
  }
  // update the expression
  std::swap(Expression_, Expression);
  Timer::stopTimer("applyEqualization");
}

void Access::applyRasterization() {
  Timer::startTimer("applyRasterization");
  // compute the extent of the domain
  std::vector<long> Sizes;
  for (int i = 0; i < Domain_.dim(isl::dim::set); ++i) {
    long LowerBound = isl::computeMinimum(Domain_, i);
    long UpperBound = isl::computeMaximum(Domain_, i);
    Sizes.push_back(1 + UpperBound - LowerBound);
  }
  // raster the dimensions to eliminate even more terms
  std::vector<piece> Expression;
  for (auto &Piece : Expression_) {
    // compute the raster information
    std::vector<int> Dimensions = findRasterDimensions(Piece);
    std::vector<isl::val> Multipliers = computeMultipliers(Dimensions, Piece);
    // filter dimensions with multiplier larger than size
    for (int i = Dimensions.size() - 1; i >= 0; --i) {
      if (Sizes[Dimensions[i]] <= isl::get_value(Multipliers[i])) {
        Dimensions.erase(Dimensions.begin() + i);
        Multipliers.erase(Multipliers.begin() + i);
      }
    }
    // oder the dimensions by the multiplier values
    std::vector<int> Indexes;
    for (int i = 0; i < Dimensions.size(); ++i) {
      Indexes.push_back(i);
    }
    std::sort(Indexes.begin(), Indexes.end(), [&](int I1, int I2) { return Multipliers[I1].lt(Multipliers[I2]); });
    // compute the domain size
    long All = getPieceSize(Piece);
    // raster the pieces dimension by dimension
    std::vector<piece> Current = {Piece};
    std::vector<piece> Next;
    // try to make all dimensions affine
    bool Updated = false;
    for (auto Index : Indexes) {
      // skip the dimension if the average subdomain would be small
      if (All / (Current.size() * isl::get_value(Multipliers[Index])) < 16)
        continue;
      // try to raster the dimension
      Next.clear();
      bool IsSuccess = true;
      for (auto Piece : Current) {
        std::vector<piece> Updates = rasterDimension(Dimensions[Index], Multipliers[Index], Piece);
        if (Updates.empty()) {
          IsSuccess = false;
          break;
        }
        Next.insert(Next.end(), Updates.begin(), Updates.end());
      }
      // swap the containers if successfull
      if (IsSuccess) {
        std::swap(Current, Next);
        Updated = true;
      }
    }
    // if the domain was split preprocess the affine pieces
    if (Updated) {
      for (auto &Piece : Current) {
        if (Piece.Expression.is_null() && isPieceAffine(Piece))
          Piece.Expression = extractAffineExpression(Piece);
      }
    }
    // store the updated pieces
    Expression.insert(Expression.end(), Current.begin(), Current.end());
  }
  // update the expression
  std::swap(Expression_, Expression);
  Timer::stopTimer("applyRasterization");
}

int Access::computeExponent(piece Piece) const {
  int Result = 0;
  // check all terms
  auto analyzeTerms = [&](isl::term Term) {
    int Exponent = 0;
    // count the parameter exponents
    for (int i = 0; i < Term.dim(isl::dim::param); i++) {
      Exponent += Term.get_exp(isl::dim::param, i);
    }
    // count the variable exponents
    for (int i = 0; i < Term.dim(isl::dim::set); i++) {
      Exponent += Term.get_exp(isl::dim::set, i);
    }
    // get the divisor variable
    for (int i = 0; i < Term.dim(isl::dim::div); i++) {
      Exponent += Term.get_exp(isl::dim::div, i);
    }
    Result = std::max(Result, Exponent);
    return isl::stat::ok();
  };
  Piece.Polynomial.foreach_term(analyzeTerms);
  return Result;
}

int Access::computeDimensionExponent(int Dimension, piece Piece) const {
  int Result = 0;
  // check all terms
  auto analyzeTerms = [&](isl::term Term) {
    int Exponent = 0;
    if (Dimension < Term.dim(isl::dim::set))
      Exponent += Term.get_exp(isl::dim::set, Dimension);
    // get the divisor variable
    for (int i = 0; i < Term.dim(isl::dim::div); i++) {
      if (Term.get_exp(isl::dim::div, i) > 0) {
        isl::aff Divisor = Term.get_div(i);
        if (Dimension < Divisor.dim(isl::dim::in) && !Divisor.get_coefficient_val(isl::dim::in, Dimension).is_zero())
          Exponent += Term.get_exp(isl::dim::div, i);
      }
    }
    Result = std::max(Result, Exponent);
    return isl::stat::ok();
  };
  Piece.Polynomial.foreach_term(analyzeTerms);
  return Result;
}

isl::qpolynomial Access::computeReplacement(std::map<int, isl::qpolynomial> Replacements, piece Piece) const {
  isl::space Space = Piece.Polynomial.get_domain_space();
  isl::qpolynomial Polynomial = isl::qpolynomial::zero_on_domain(Space);
  auto updatePolynomial = [&](isl::term Term) {
    bool Replace = false;
    for (int i = 0; i < Term.dim(isl::dim::div); i++) {
      if (Term.get_exp(isl::dim::div, i) > 0 && Replacements.count(i) == 1)
        Replace = true;
    }
    // compute replacement term
    if (Replace) {
      // compute the update term
      isl::qpolynomial Update = isl::qpolynomial::val_on_domain(Space, Term.get_coefficient_val());
      // multiply with the parameters
      for (int i = 0; i < Term.dim(isl::dim::param); ++i) {
        int Exponent = Term.get_exp(isl::dim::param, i);
        Update = Update.mul(isl::qpolynomial::var_on_domain(Space, isl::dim::param, i).pow(Exponent));
      }
      // multiply with variables
      for (int i = 0; i < Term.dim(isl::dim::set); ++i) {
        int Exponent = Term.get_exp(isl::dim::set, i);
        Update = Update.mul(isl::qpolynomial::var_on_domain(Space, isl::dim::set, i).pow(Exponent));
      }
      // compute the divisors
      for (int i = 0; i < Term.dim(isl::dim::div); ++i) {
        int Exponent = Term.get_exp(isl::dim::div, i);
        if (Exponent > 0) {
          if (Replacements.count(i) == 1) {
            Update = Update.mul(Replacements[i].pow(Exponent));
          } else
            Update = Update.mul(isl::qpolynomial::from_aff(Term.get_div(i).floor()).pow(Exponent));
        }
      }
      // update the polynomial
      Polynomial = Polynomial.add(Update);
    } else {
      Polynomial = Polynomial.add(isl::qpolynomial::from_term(Term));
    }
    return isl::stat::ok();
  };
  Piece.Polynomial.foreach_term(updatePolynomial);
  return Polynomial;
}

std::vector<int> Access::findRasterDimensions(piece Piece) const {
  // count the dimension exponent pairs for all and divisor terms
  std::map<std::pair<int, int>, int> Counts;
  std::map<std::pair<int, int>, int> DCounts;
  auto analyzeTerms = [&](isl::term Term) {
    // compute the variable exponents
    int Parameters = 0;
    std::map<int, int> Variables;
    std::map<int, int> DVariables;
    // count the parameter exponents
    for (int i = 0; i < Term.dim(isl::dim::param); i++) {
      Parameters += Term.get_exp(isl::dim::param, i);
    }
    // count the variable exponents
    for (int i = 0; i < Term.dim(isl::dim::set); i++) {
      Variables[i] += Term.get_exp(isl::dim::set, i);
    }
    // get the divisor variable
    for (int i = 0; i < Term.dim(isl::dim::div); i++) {
      if (Term.get_exp(isl::dim::div, i) > 0) {
        isl::aff Divisor = Term.get_div(i);
        // extract the parameter and variable information
        for (int j = 0; j < Divisor.dim(isl::dim::param); ++j) {
          isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::param, j);
          if (!Coefficient.is_zero())
            Parameters += Term.get_exp(isl::dim::div, i);
        }
        for (int j = 0; j < Divisor.dim(isl::dim::in); ++j) {
          isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::in, j);
          if (!Coefficient.is_zero()) {
            Variables[j] += Term.get_exp(isl::dim::div, i);
            DVariables[j] += Term.get_exp(isl::dim::div, i);
          }
        }
      }
    }
    // update the counts for parameter free terms
    if (Parameters == 0) {
      for (auto Variable : Variables) {
        Counts[Variable]++;
        if (DVariables[Variable.first] > 0)
          DCounts[Variable]++;
      }
    }
    return isl::stat::ok();
  };
  Piece.Polynomial.foreach_term(analyzeTerms);
  // compute the candidate dimensions
  std::vector<int> Dimensions;
  for (auto Count : Counts) {
    // select the dimensions that appear multiple times with exponent > 1
    // only consider dimensions with divisors
    if (Count.second > 1 && Count.first.second > 1 && DCounts[Count.first] > 0)
      Dimensions.push_back(Count.first.first);
  }
  // make dimensions unique
  std::sort(Dimensions.begin(), Dimensions.end());
  Dimensions.erase(std::unique(Dimensions.begin(), Dimensions.end()), Dimensions.end());
  return Dimensions;
}

std::vector<isl::val> Access::computeMultipliers(std::vector<int> Dimensions, piece Piece) const {
  // compute the least common multipliers for all dimensions
  std::vector<long> LCM(Piece.Polynomial.dim(isl::dim::in), 1);
  // process every divisor exactly once
  std::set<int> Divisors;
  auto findMultipliers = [&](isl::term Term) {
    // get the divisor variable
    for (int i = 0; i < Term.dim(isl::dim::div); i++) {
      if (Divisors.count(i) == 0 && Term.get_exp(isl::dim::div, i) > 0) {
        Divisors.insert(i);
        isl::aff Divisor = Term.get_div(i);
        // count parameters and find variables
        std::map<int, isl::val> Variables;
        // extract the variable information
        for (int j = 0; j < Divisor.dim(isl::dim::in); ++j) {
          isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::in, j);
          if (!Coefficient.is_zero())
            Variables[j] = Coefficient;
        }
        // compute the multipliers for all variables of the term
        for (auto Dimension : Dimensions) {
          if (Variables.count(Dimension) == 1) {
            long Coefficient = isl::get_value(Variables[Dimension].mul(Divisor.get_denominator_val()));
            long Denominator = isl::get_value(Divisor.get_denominator_val());
            Denominator = Denominator / boost::math::gcd(Coefficient, Denominator);
            LCM[Dimension] = boost::math::lcm(LCM[Dimension], Denominator);
          }
        }
      }
    }
    return isl::stat::ok();
  };
  Piece.Polynomial.foreach_term(findMultipliers);
  // compute the multipliers
  std::vector<isl::val> Multipliers;
  for (auto Dimension : Dimensions) {
    Multipliers.push_back(isl::val(Piece.Domain.get_ctx(), LCM[Dimension]));
  }
  return Multipliers;
}

std::vector<piece> Access::rasterDimension(int Dimension, isl::val Multiplier, piece Piece) const {
  isl::local_space DS = isl::local_space(Piece.Domain.get_space());
  isl::local_space PS = isl::local_space(Piece.Polynomial.get_domain_space());
  // result vector
  std::vector<piece> Pieces;
  // raster the domain with all possible modulo value
  for (int Current = 0; Current < isl::get_value(Multiplier); ++Current) {
    // compute the domain
    isl::set Domain = Piece.Domain;
    isl::pw_aff Modulo = isl::pw_aff::var_on_domain(DS, isl::dim::set, Dimension) % Multiplier;
    isl::set Element = Modulo.eq_set(isl::aff(DS, isl::val(Domain.get_ctx(), Current)));
    Domain = Domain.intersect(Element);
    // compute the polynomial
    isl::qpolynomial Polynomial;
    if (!Domain.is_empty()) {
      // compute replacement divisors
      std::map<int, isl::qpolynomial> Replacements;
      std::set<int> Divisors;
      auto findReplacement = [&](isl::term Term) {
        // does the term contain the critical divisors
        for (int i = 0; i < Term.dim(isl::dim::div); i++) {
          if (Divisors.count(i) == 0 && Term.get_exp(isl::dim::div, i) > 0) {
            Divisors.insert(i);
            isl::aff Divisor = Term.get_div(i);
            // compute the variable information
            int Parameters = 0;
            std::vector<int> Variables;
            std::vector<isl::val> Coefficients;
            // extract the parameter and variable information
            for (int j = 0; j < Divisor.dim(isl::dim::param); ++j) {
              isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::param, j);
              if (!Coefficient.is_zero())
                Parameters++;
            }
            for (int j = 0; j < Divisor.dim(isl::dim::in); ++j) {
              isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::in, j);
              if (!Coefficient.is_zero()) {
                Variables.push_back(j);
                Coefficients.push_back(Coefficient);
              }
            }
            // search the dimension
            int Index = -1;
            auto Iter = std::find(Variables.begin(), Variables.end(), Dimension);
            if (Iter != Variables.end())
              Index = std::distance(Variables.begin(), Iter);
            // replace the divisor if the variable matches
            if (Parameters == 0 && Index >= 0) {
              // compute the replacement
              isl::aff Remainder = isl::aff(PS, isl::val(Domain.get_ctx(), Current));
              isl::aff New = isl::aff(PS, Divisor.get_constant_val());
              // add all other variables
              for (int j = 0; j < Variables.size(); ++j) {
                if (j != Index) {
                  isl::aff Var = isl::aff::var_on_domain(PS, isl::dim::set, Variables[j]);
                  New = New.add(isl::aff(PS, Coefficients[j]).mul(Var));
                }
              }
              // add the remainder and compute the floor
              New = New.add(isl::aff(PS, Coefficients[Index]).mul(Remainder)).floor();
              // add the division outside of the floor
              isl::aff Var = isl::aff::var_on_domain(PS, isl::dim::set, Variables[Index]).sub(Remainder);
              New = New.add(isl::aff(PS, Coefficients[Index]).mul(Var));
              Replacements[i] = isl::qpolynomial::from_aff(New);
            }
          }
        }
        return isl::stat::ok();
      };
      Piece.Polynomial.foreach_term(findReplacement);
      // compute the polynomail
      Polynomial = computeReplacement(Replacements, Piece);
      // store the result
      Pieces.push_back(createPiece(Domain, Polynomial));
    }
  }
  // verify that at least one piece is affine along the dimension
  bool IsAffine = false;
  for (auto Piece : Pieces) {
    if (computeDimensionExponent(Dimension, Piece) <= 1) {
      IsAffine = true;
      break;
    }
  }
  if (!IsAffine)
    return {};
#ifdef VERIFY_RESULT
  // verify the correctness of the elimination
  bool Success = verifySplit(Piece, Pieces);
  assert(Success);
#endif
  return Pieces;
}

std::vector<std::vector<std::tuple<int, long, long>>> Access::findEqualizationCandidates(piece Piece) const {
  // count the number of divisors that differ only by constant
  std::map<std::vector<std::tuple<int, long, long>>, int> Counts;
  std::vector<std::vector<std::tuple<int, long, long>>> Candidates;
  // process every divisor exactly once
  std::set<int> Divisors;
  // analyze the terms
  auto analyzeTerms = [&](isl::term Term) {
    // get the divisor variable
    for (int i = 0; i < Term.dim(isl::dim::div); i++) {
      if (Divisors.count(i) == 0 && Term.get_exp(isl::dim::div, i) > 0) {
        Divisors.insert(i);
        isl::aff Divisor = Term.get_div(i);
        int Parameters = 0;
        std::vector<std::tuple<int, long, long>> Variables;
        isl::val Denominator = Divisor.get_denominator_val();
        // extract the parameter and variable information
        for (int j = 0; j < Divisor.dim(isl::dim::param); ++j) {
          isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::param, j);
          if (!Coefficient.is_zero())
            Parameters++;
        }
        for (int j = 0; j < Divisor.dim(isl::dim::in); ++j) {
          isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::in, j);
          if (!Coefficient.is_zero())
            Variables.push_back(
                std::make_tuple(j, isl::get_value(Coefficient.mul(Denominator)), isl::get_value(Denominator)));
        }
        // count the single variate divisors with coefficient one
        if (Parameters == 0) {
          Counts[Variables]++;
        }
      }
    }
    return isl::stat::ok();
  };
  Piece.Polynomial.foreach_term(analyzeTerms);
  // convert the divisors that apply multiple times
  for (auto Count : Counts) {
    if (Count.second > 1) {
      Candidates.push_back(Count.first);
    }
  }
  return Candidates;
}

std::vector<long> Access::computeSplits(std::vector<std::tuple<int, long, long>> Candidate, piece Piece) const {
  // compute the splits for the candidate variable
  std::vector<long> Splits;
  std::set<int> Divisors;
  auto computeSplits = [&](isl::term Term) {
    // get the divisor variable
    for (int i = 0; i < Term.dim(isl::dim::div); i++) {
      if (Divisors.count(i) == 0 && Term.get_exp(isl::dim::div, i) > 0) {
        Divisors.insert(i);
        isl::aff Divisor = Term.get_div(i);
        int Parameters = 0;
        std::vector<int> Variables;
        std::vector<isl::val> Coefficients;
        isl::val Denominator = Divisor.get_denominator_val();
        // extract the parameter and variable information
        for (int j = 0; j < Divisor.dim(isl::dim::param); ++j) {
          isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::param, j);
          if (!Coefficient.is_zero())
            Parameters++;
        }
        for (int j = 0; j < Divisor.dim(isl::dim::in); ++j) {
          isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::in, j);
          if (!Coefficient.is_zero()) {
            Variables.push_back(j);
            Coefficients.push_back(Coefficient);
          }
        }
        // make sure the term matches the candidate
        if (Parameters == 0 && Variables.size() == Candidate.size()) {
          bool IsMatch = true;
          for (int i = 0; i < Variables.size(); ++i) {
            if (Variables[i] != std::get<0>(Candidate[i]) ||
                isl::get_value(Coefficients[i].mul(Denominator)) != std::get<1>(Candidate[i]) ||
                isl::get_value(Denominator) != std::get<2>(Candidate[i])) {
              IsMatch = false;
              break;
            }
          }
          // if all variables match compute the splits
          if (IsMatch) {
            isl::val Constant = Divisor.get_constant_val().mul(Denominator);
            Splits.push_back(isl::get_value(Denominator.sub(Constant).mod(Denominator)));
          }
        }
      }
    }
    return isl::stat::ok();
  };
  Piece.Polynomial.foreach_term(computeSplits);
  // make the splits unique
  std::sort(Splits.begin(), Splits.end());
  Splits.erase(std::unique(Splits.begin(), Splits.end()), Splits.end());
  return Splits;
}

std::vector<piece> Access::equalizeCandidate(std::vector<std::tuple<int, long, long>> Candidate,
                                             std::vector<long> Splits, piece Piece) const {
  isl::local_space DS = isl::local_space(Piece.Domain.get_space());
  isl::local_space PS = isl::local_space(Piece.Polynomial.get_domain_space());
  // result vector
  std::vector<piece> Pieces;
  assert(!Splits.empty());
  // extract the candidate information
  isl::val Denominator = isl::val(Piece.Domain.get_ctx(), std::get<2>(Candidate[0]));
  isl::aff Dimension = isl::aff(DS, isl::val(Piece.Domain.get_ctx(), 0));
  for (auto Variable : Candidate) {
    Dimension = Dimension.add(isl::aff::var_on_domain(DS, isl::dim::set, std::get<0>(Variable))
                                  .mul(isl::aff(DS, isl::val(Piece.Domain.get_ctx(), std::get<1>(Variable)))));
  }
  // compute the equalization offset
  // (allows us to avoid split at zero)
  isl::val Offset = isl::val(Piece.Domain.get_ctx(), -Splits[0]).mod(Denominator);
  // add the correction terms for each interval
  long Start = Splits.back();
  for (int i = 0; i < Splits.size(); ++i) {
    // compute the range
    long Stop = Splits[i];
    bool Single = Stop - Start == 1 || Stop - Start == 1 - isl::get_value(Denominator);
    // compute the domain
    isl::set Domain = Piece.Domain;
    if (Start != Stop) {
      isl::aff Modulo = Dimension.mod(Denominator);
      isl::set LowerBound = Modulo.ge_set(isl::aff(DS, isl::val(Domain.get_ctx(), Start)));
      isl::set UpperBound = Modulo.lt_set(isl::aff(DS, isl::val(Domain.get_ctx(), Stop)));
      // handle intervals that wrap around
      if (Start < Stop) {
        Domain = Domain.intersect(LowerBound).intersect(UpperBound);
      } else {
        Domain = Domain.intersect(LowerBound).unite(Domain.intersect(UpperBound));
      }
    }
    // compute the polynomial if the domain is not empty
    isl::qpolynomial Polynomial;
    if (!Domain.is_empty()) {
      // compute replacement divisors
      std::map<int, isl::qpolynomial> Replacements;
      std::set<int> Divisors;
      auto findReplacement = [&](isl::term Term) {
        // does the term contain the critical divisors
        for (int i = 0; i < Term.dim(isl::dim::div); i++) {
          if (Divisors.count(i) == 0 && Term.get_exp(isl::dim::div, i) > 0) {
            Divisors.insert(i);
            isl::aff Divisor = Term.get_div(i);
            int Parameters = 0;
            std::vector<int> Variables;
            std::vector<isl::val> Coefficients;
            isl::val Denominator = Divisor.get_denominator_val();
            // extract the parameter and variable information
            for (int j = 0; j < Divisor.dim(isl::dim::param); ++j) {
              isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::param, j);
              if (!Coefficient.is_zero())
                Parameters++;
            }
            for (int j = 0; j < Divisor.dim(isl::dim::in); ++j) {
              isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::in, j);
              if (!Coefficient.is_zero()) {
                Variables.push_back(j);
                Coefficients.push_back(Coefficient);
              }
            }
            // check if we match the candidate
            if (Parameters == 0 && Variables.size() == Candidate.size()) {
              bool IsMatch = true;
              for (int i = 0; i < Variables.size(); ++i) {
                if (Variables[i] != std::get<0>(Candidate[i]) ||
                    isl::get_value(Coefficients[i].mul(Denominator)) != std::get<1>(Candidate[i]) ||
                    isl::get_value(Denominator) != std::get<2>(Candidate[i])) {
                  IsMatch = false;
                  break;
                }
              }
              // if all variables match compute the splits
              if (IsMatch) {
                isl::val Constant = Divisor.get_constant_val().mul(Denominator);
                isl::val Variable = isl::val(Term.get_ctx(), Start);
                isl::val Old = Variable.add(Constant).div(Denominator).floor();
                // compute the replacement divisor
                // (no floor necessary if for single elements)
                if (Single) {
                  isl::val New = Variable.add(Offset).div(Denominator);
                  Replacements[i] = isl::qpolynomial::from_aff(Dimension.add(isl::aff(PS, Offset))
                                                                   .div(isl::aff(PS, Denominator))
                                                                   .add(isl::aff(PS, Old.sub(New))));
                } else {
                  isl::val New = Variable.add(Offset).div(Denominator).floor();
                  Replacements[i] = isl::qpolynomial::from_aff(Dimension.add(isl::aff(PS, Offset))
                                                                   .div(isl::aff(PS, Denominator))
                                                                   .floor()
                                                                   .add(isl::aff(PS, Old.sub(New))));
                }
              }
            }
          }
        }
        return isl::stat::ok();
      };
      Piece.Polynomial.foreach_term(findReplacement);
      // update the polynomial
      Polynomial = computeReplacement(Replacements, Piece);
    }
    // create the piece
    Pieces.push_back(createPiece(Domain, Polynomial));
    Start = Stop;
  }
  // check if the number of divisors improves
  bool Improved = false;
  auto Exponent = computeExponent(Piece);
  for (auto Piece : Pieces) {
    if (!Piece.Polynomial.is_null()) {
      if (computeExponent(Piece) < Exponent) {
        Improved = true;
        break;
      }
    } else {
      Improved = true;
      break;
    }
  }
  if (!Improved) {
    return {};
  }
#ifdef VERIFY_RESULT
  // verify the correctness of the elimination
  bool Success = verifySplit(Piece, Pieces);
  assert(Success);
#endif
  return Pieces;
}

bool Access::verifySplit(piece Piece, std::vector<piece> Updates) const {
  // sum the domain sizes of the pieces and compare the the original piece
  long All = isl::cardinality(Piece.Domain);
  long Total = std::accumulate(Updates.begin(), Updates.end(), 0,
                               [](long Total, piece Current) { return Total + isl::cardinality(Current.Domain); });
  if (All != Total) {
    printf("-> domain split failed (%ld!=%ld)!\n", All, Total);
    Piece.Domain.dump();
    printf("transformed to");
    for (auto Update : Updates) {
      Update.Domain.dump();
    }
    return false;
  }
  // verify the polynomials
  for (auto Update : Updates) {
    bool Match = true;
    if (!Update.Polynomial.is_null()) {
      // compare the polynomials on all points of the domain
      auto comparePolynomials = [&](isl::point Point) {
        isl::val Before = Piece.Polynomial.eval(Point);
        isl::val After = Update.Polynomial.eval(Point);
        if (!Before.eq(After)) {
          Match = false;
        }
        return isl::stat::ok();
      };
      Update.Domain.foreach_point(comparePolynomials);
      if (!Match) {
        printf("-> polynomial split failed!\n");
        Piece.Polynomial.dump();
        Piece.Domain.dump();
        printf("transformed to\n");
        Update.Polynomial.dump();
        Update.Domain.dump();
        printf("Updates.size()==%ld\n", Updates.size());
        return false;
      }
    }
  }
  return true;
}

piece Access::createPiece(isl::set Domain, isl::qpolynomial Count) const {
  // analyze the next piece
  piece Piece;
  Piece.Domain = Domain;
  Piece.Size = -1;
  Piece.Polynomial = Count;
  // analyze the terms
  auto analyzeTerms = [&](isl::term Term) {
    // analyze the term
    term Result;
    // get the coefficient
    Result.Coefficient = Term.get_coefficient_val();
    // get the parameter exponents
    for (int i = 0; i < Term.dim(isl::dim::param); ++i) {
      Result.Parameters.push_back(Term.get_exp(isl::dim::param, i));
    }
    // get the input exponents
    for (int i = 0; i < Term.dim(isl::dim::set); i++) {
      Result.Variables.push_back(Term.get_exp(isl::dim::set, i));
    }
    // get the divisors exponents
    for (int i = 0; i < Term.dim(isl::dim::div); i++) {
      Result.Divisors.push_back(Term.get_exp(isl::dim::div, i));
    }
    Result.Polynomial = isl::qpolynomial::from_term(Term);
    Piece.Terms.push_back(Result);
    return isl::stat::ok();
  };
  Count.foreach_term(analyzeTerms);
  return Piece;
}

void Access::extractStackDistanceExpression(isl::union_pw_qpolynomial Count) {
  // analyze the polynomials
  auto analyzePolynomial = [&](isl::pw_qpolynomial Polynomial) {
    // try to minimize the number of pieces
    Polynomial = Polynomial.coalesce();
    // analyze th individual piece
    auto analyzePiece = [&](isl::set Domain, isl::qpolynomial Polynomial) {
      // analyze the next piece
      auto Piece = createPiece(Domain, Polynomial);
      if (isPieceAffine(Piece)) {
        Piece.Expression = extractAffineExpression(Piece);
      }
      Expression_.push_back(Piece);
      return isl::stat::ok();
    };
    Polynomial.foreach_piece(analyzePiece);
    return isl::stat::ok();
  };
  Count.foreach_pw_qpolynomial(analyzePolynomial);
}

bool Access::isPieceAffine(piece Piece) const {
  // sum the exponents the determine if the piece is affine
  // (note that this function does not consider domains with constant dimensions)
  for (auto Term : Piece.Terms) {
    // count the factors
    int Factors = 0;
    // check the parameters
    for (int i = 0; i < Term.Parameters.size(); ++i) {
      if (Term.Parameters[i] > 0) {
        Factors += Term.Parameters[i];
        if (Factors > 1)
          return false;
      }
    }
    // multiply with variable
    for (int i = 0; i < Term.Variables.size(); ++i) {
      if (Term.Variables[i] > 0) {
        Factors += Term.Variables[i];
        if (Factors > 1)
          return false;
      }
    }
    // multiply with divisor
    for (int i = 0; i < Term.Divisors.size(); ++i) {
      if (Term.Divisors[i] > 0) {
        Factors += Term.Divisors[i];
        if (Factors > 1)
          return false;
      }
    }
  }
  return true;
}

long Access::getPieceSize(piece &Piece) const {
  if (Piece.Size == -1) {
    Piece.Size = isl::cardinality(Piece.Domain);
  }
  return Piece.Size;
}

isl::pw_aff Access::extractAffineExpression(piece Piece) const {
  Timer::startTimer("extractAffineExpression");
  isl::local_space DS = isl::local_space(Piece.Domain.get_space());
  // initialize the expression with the neutral element
  isl::val Zero = isl::val::zero(Piece.Domain.get_ctx());
  isl::pw_aff Sum(Piece.Domain, Zero);
  for (auto &Term : Piece.Terms) {
    // compute the product
    isl::pw_aff Product(Piece.Domain, Term.Coefficient);
    // multiply with parameter
    for (int i = 0; i < Term.Parameters.size(); ++i) {
      if (Term.Parameters[i] > 0) {
        assert(Term.Parameters[i] <= 1);
        isl::pw_aff Parameter = isl::pw_aff::var_on_domain(DS, isl::dim::param, i);
        Product = Product.mul(Parameter);
      }
    }
    // multiply with variable
    for (int i = 0; i < Term.Variables.size(); ++i) {
      if (Term.Variables[i] > 0) {
        assert(Term.Variables[i] <= 1);
        isl::pw_aff Variable;
        Variable = isl::pw_aff::var_on_domain(DS, isl::dim::set, i);
        Product = Product.mul(Variable);
      }
    }
    // multiply with divisor
    auto multiplyTerms = [&](isl::term Term) {
      for (int i = 0; i < Term.dim(isl::dim::div); i++) {
        if (Term.get_exp(isl::dim::div, i) > 0) {
          assert(Term.get_exp(isl::dim::div, i) == 1);
          Product = Product.mul(Term.get_div(i).floor());
        }
      }
      return isl::stat::ok();
    };
    Term.Polynomial.foreach_term(multiplyTerms);
    Sum = Sum.add(Product);
  }
  Timer::stopTimer("extractAffineExpression");
  return Sum;
}

isl::aff Access::extractAffineExpression(isl::qpolynomial Polynomial, isl::set Domain,
                                         std::map<int, long> Values) const {
  Timer::startTimer("extractAffineExpression");
  // buffer the fixed divisors
  std::map<int, isl::aff> Divisors;
  // accumulate the values of all terms
  isl::local_space DS = isl::local_space(Domain.get_space());
  isl::val Zero = isl::val::zero(Domain.get_ctx());
  isl::aff Sum(DS, Zero);
  auto analyzeTerms = [&](isl::term Term) {
    auto Coefficient = Term.get_coefficient_val();
    isl::aff Product(DS, Coefficient);
    // get the parameter exponents
    for (int i = 0; i < Term.dim(isl::dim::param); ++i) {
      int Exponent = Term.get_exp(isl::dim::param, i);
      assert(Exponent <= 1);
      if (Exponent == 1) {
        auto Parameter = isl::aff::var_on_domain(DS, isl::dim::param, i);
        Product = Product.mul(Parameter);
      }
    }
    // get the input exponents
    for (int i = 0; i < Term.dim(isl::dim::set); i++) {
      int Exponent = Term.get_exp(isl::dim::set, i);
      if (Exponent > 0) {
        if (Values.count(i)) {
          long Coefficient = compute_power(Values[i], Exponent);
          isl::aff Value(DS, isl::val(Domain.get_ctx(), Coefficient));
          Product = Product.mul(Value);
        } else {
          assert(Exponent == 1);
          auto Variable = isl::aff::var_on_domain(DS, isl::dim::set, i);
          Product = Product.mul(Variable);
        }
      }
    }
    // get the divisors exponents
    for (int i = 0; i < Term.dim(isl::dim::div); i++) {
      int Exponent = Term.get_exp(isl::dim::div, i);
      if (Exponent > 0) {
        // update the divisor if we don't have it yet
        if (Divisors.count(i) == 0) {
          auto Divisor = Term.get_div(i);
          auto Replacement = isl::aff(DS, Divisor.get_constant_val());
          for (int j = 0; j < Divisor.dim(isl::dim::param); ++j) {
            isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::param, j);
            if (!Coefficient.is_zero()) {
              Replacement = Replacement.add_coefficient_val(isl::dim::param, j, Coefficient);
            }
          }
          for (int j = 0; j < Divisor.dim(isl::dim::in); ++j) {
            isl::val Coefficient = Divisor.get_coefficient_val(isl::dim::in, j);
            if (!Coefficient.is_zero()) {
              if (Values.count(j)) {
                Replacement = Replacement.add_constant_val(Coefficient.mul_ui(Values[j]));
              } else {
                Replacement = Replacement.add_coefficient_val(isl::dim::in, j, Coefficient);
              }
            }
          }
          Divisors[i] = Replacement;
        }
        // distinguish cases
        if (Divisors[i].is_cst()) {
          long Base = isl::get_value(Divisors[i].floor().get_constant_val());
          long Coefficient = compute_power(Base, Exponent);
          isl::aff Value(DS, isl::val(Domain.get_ctx(), Coefficient));
          Product = Product.mul(Value);
        } else {
          assert(Exponent == 1);
          Product = Product.mul(Divisors[i].floor());
        }
      }
    }
    Sum = Sum.add(Product);
    return isl::stat::ok();
  };
  Polynomial.foreach_term(analyzeTerms);
  Timer::stopTimer("extractAffineExpression");
  return Sum;
}
