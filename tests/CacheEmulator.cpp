/*
* Copyright (c) 2019, ETH Zurich
*/

#include <algorithm>

#include "CacheEmulator.h"

#include <cassert>
#include <iostream>
#include <iterator>

void CacheEmulator::accessMemory(std::string Statement, int &TimeStamp, std::vector<int> CacheLines) {
  // initialize the result if necessary
  if (StackDistances_[Statement].size() == 0)
    StackDistances_[Statement].resize(CacheLines.size(), 0);
  if (CapacityMisses_[Statement].size() == 0)
    CapacityMisses_[Statement].resize(CacheLines.size(), 0);
  if (CompulsoryMisses_[Statement].size() == 0)
    CompulsoryMisses_[Statement].resize(CacheLines.size(), 0);
  // process the accesses one after the other
  for (int i = 0; i < CacheLines.size(); i++) {
    // get the last access
    assert(CacheLines[i] < TimeStamps_.size());
    long LastAccess = TimeStamps_[CacheLines[i]];
    // count the number of cache lines touched since the last access
    long StackDistance = 0;
    for (long j = 0; j < TimeStamps_.size(); ++j) {
      if (LastAccess >= 0 && LastAccess <= TimeStamps_[j])
        StackDistance++;
    }
    // update the time stamp
    TimeStamps_[CacheLines[i]] = TimeStamp;
    // store the results
    if (StackDistance > StackDistances_[Statement][i])
      StackDistances_[Statement][i] = StackDistance;
    if (LastAccess == -1)
      CompulsoryMisses_[Statement][i]++;
    if (StackDistance > CacheSize_)
      CapacityMisses_[Statement][i]++;
    // increment the time stamp
    TimeStamp++;
  }
}

std::map<std::string, std::vector<long>> CacheEmulator::getStackDistances() const {
  // remove zero entries
  std::map<std::string, std::vector<long>> Results;
  std::copy_if(StackDistances_.begin(), StackDistances_.end(), std::inserter(Results, Results.end()),
               [](decltype(Results)::value_type const &Statement) {
                 return *std::max_element(Statement.second.begin(), Statement.second.end()) > 0;
               });
  return Results;
}

std::map<std::string, std::vector<long>> CacheEmulator::getCapacityMisses() const { return CapacityMisses_; }

std::map<std::string, std::vector<long>> CacheEmulator::getCompulsoryMisses() const { return CompulsoryMisses_; }