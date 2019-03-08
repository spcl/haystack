/*
* Copyright (c) 2019, ETH Zurich
*/

#ifndef _CHACHE_EMULATOR_H_
#define _CHACHE_EMULATOR_H_

#include <map>
#include <string>
#include <vector>

#include "../src/Definitions.h"

class CacheEmulator {
public:
  CacheEmulator() = delete;
  CacheEmulator(int CacheLines, int CacheSize) : TimeStamps_(CacheLines, -1), CacheSize_(CacheSize) {}

  void accessMemory(std::string Statement, int &TimeStamp, std::vector<int> CacheLines);

  std::map<std::string, std::vector<long>> getStackDistances() const;
  std::map<std::string, std::vector<long>> getCapacityMisses() const;
  std::map<std::string, std::vector<long>> getCompulsoryMisses() const;

private:
  int CacheSize_;
  // time stamp per cache line
  std::vector<long> TimeStamps_;

  // cache information per statement
  std::map<std::string, std::vector<long>> StackDistances_;
  std::map<std::string, std::vector<long>> CapacityMisses_;
  std::map<std::string, std::vector<long>> CompulsoryMisses_;
};

#endif