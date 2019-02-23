#ifndef _ACCESS_H_
#define _ACCESS_H_

#include <isl/isl-noexceptions.h>
#include <map>
#include <string>
#include <vector>

#include "Definitions.h"
#include "Program.h"

// class holding the access information
class Access {
public:
  Access() = delete;
  Access(const Access &other) = default;
  Access(std::string Name, machine_model MachineModel, isl::set Domain, std::map<std::string, long> ElementSizes,
         prefetch_info Prefetched)
      : Name_(Name), MachineModel_(MachineModel), Domain_(Domain), ElementSizes_(ElementSizes),
        Prefetched_(Prefetched) {}

  // control the cache miss computation
  void initAccess(std::vector<NamedLong> ParameterValues, isl::set Parameters);

  // compute base cache miss information
  void countCompulsoryMisses(isl::union_map First);
  void computeStackDistances(isl::union_map BetweenMap);
  void countCapacityMisses();
  misses getResult() const { return Result_; };

  // compute the capacity misses for additional cache sizes
  std::vector<long> countCapacityMisses(std::vector<long> CacheSizes);

  // get the access properties
  std::string getName() const { return Name_; }
  isl::set getDomain() const { return Domain_; }

private:
  // helper methods
  piece createPiece(isl::set Domain, isl::qpolynomial Piece) const;
  void extractStackDistanceExpression(isl::union_pw_qpolynomial Count);
  void storeAffinePieces();

  void applyEqualization();
  void applyRasterization();

  // enumerate non-affine dimensions to count all points
  std::vector<int> findNonAffineDimensions(piece Piece) const;
  void enumerateNonAffineDimensions(piece Piece);

  // counting methods
  std::vector<long> countAffineDimensions(piece Piece, std::vector<long> Limits) const;
  std::vector<long> enumerateNonAffinePoints(piece Piece, std::vector<long> Limits) const;

  // helper functions to analyze piece and extract the affine expression
  bool isPieceAffine(piece Piece) const;
  isl::pw_aff extractAffineExpression(piece Piece) const;
  isl::aff extractAffineExpression(isl::qpolynomial, isl::set Domain, std::map<int, long> Values) const;
  long getPieceSize(piece &Piece) const;

  // elimination helper methods
  int computeExponent(piece Piece) const;
  int computeDimensionExponent(int Dimension, piece Piece) const;

  isl::qpolynomial computeReplacement(std::map<int, isl::qpolynomial> Replacements, piece Piece) const;

  // methods to eliminate floor terms due to equalization
  std::vector<std::vector<std::tuple<int, long, long>>> findEqualizationCandidates(piece Piece) const;
  std::vector<long> computeSplits(std::vector<std::tuple<int, long, long>> Candidate, piece Piece) const;
  std::vector<piece> equalizeCandidate(std::vector<std::tuple<int, long, long>> Candidate, std::vector<long> Splits,
                                       piece Piece) const;

  // methods to eliminate floor terms due to rasterization
  std::vector<int> findRasterDimensions(piece Piece) const;
  std::vector<isl::val> computeMultipliers(std::vector<int> Dimensions, piece Piece) const;
  std::vector<piece> rasterDimension(int Dimension, isl::val Multiplier, piece Piece) const;

  // method to verify splits
  bool verifySplit(piece Piece, std::vector<piece> Pieces) const;

  std::string Name_;
  isl::set Domain_;
  machine_model MachineModel_;
  std::map<std::string, long> ElementSizes_;
  prefetch_info Prefetched_;

  isl::set Parameters_;
  std::vector<NamedLong> ParameterValues_;
  std::vector<piece> Expression_;

  // store the stack distance
  long Misses_;
  std::vector<piece> Affine_;
  std::vector<piece> NonAffine_;
  std::vector<std::pair<long, isl::val>> Constant_;

  // result of the cache miss computation
  misses Result_;
};

// support sorting
inline bool operator<(Access const &A, Access const &B) { return A.getName() < B.getName(); }

#endif