/*
* Copyright (c) 2019, ETH Zurich
*/

#ifndef _DEFINITIONS_H_
#define _DEFINITIONS_H_

#include <cassert>
#include <isl/isl-noexceptions.h>
#include <map>
#include <string>
#include <utility>
#include <vector>

// computation options
#define RASTERIZATION 1
#define EQUALIZATION 1
#define DIMENSION_COUNTING 1
#define COMPUTE_CONFLICTS 1

//#define TIMERS 1
//#define ENUMERATE_POINTS 1

// flags to control the verification steps
// #ifndef NDEBUG
// #define VERIFY_RESULT 1
// #endif

// struct defining the machine properties
struct machine_model {
  long CacheLineSize;
  std::vector<long> CacheSizes;
};

// struct defining the model options
struct model_options {
  bool ComputeBounds;
};

// struct holding the term information
struct term {
  isl::val Coefficient;
  // the exponents
  std::vector<int> Parameters;
  std::vector<int> Variables;
  std::vector<int> Divisors;
  // the term
  isl::qpolynomial Polynomial;
};

// struct holding the piece information
struct piece {
  isl::set Domain;
  long Size;
  // polynomial
  isl::qpolynomial Polynomial;
  // terms of the polynomial
  std::vector<term> Terms;
  // affine form
  isl::pw_aff Expression;
};

// struct holding the cache misses and compute statistics
struct misses {
  long Total;
  long Counted;
  long CompulsoryMisses;
  std::vector<long> CapacityMisses;
};

// struct holding access information
enum AccessType { Read, Write };
struct access_info {
  std::string Name;
  std::string Access;
  AccessType ReadOrWrite;
  unsigned Start;
  unsigned Stop;
  int Line;
};

// define named types
typedef std::pair<std::string, int> NamedInt;
typedef std::pair<std::string, long> NamedLong;
typedef std::pair<std::string, std::vector<long>> NamedVector;
typedef std::pair<std::string, misses> NamedMisses;

// compute the integer power
inline long compute_power(long base, int exponent) {
  assert(exponent >= 0);
  long Result = 1;
  while (exponent) {
    if (exponent & 1)
      Result *= base;
    base *= base;
    exponent /= 2;
  }
  return Result;
};

#endif