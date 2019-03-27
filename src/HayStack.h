/*
 * Copyright (c) 2019, ETH Zurich
 */

#ifndef _HAY_STACK_H_
#define _HAY_STACK_H_

#include <map>
#include <string>
#include <vector>

#include "Access.h"
#include "Definitions.h"
#include "Program.h"
#include <isl/isl-noexceptions.h>

// allocate context with include paths
isl_ctx *allocateContextWithIncludePaths(std::vector<std::string> IncludePaths);

// main class of the cache model
class HayStack {
public:
  HayStack() = delete;
  HayStack(const HayStack &other) = default;
  HayStack(isl::ctx Context, machine_model MachineModel, model_options ModelOptions)
      : Context_(Context), MachineModel_(MachineModel), ModelOptions_(ModelOptions), Program_(Context, MachineModel){};
  HayStack(isl::ctx Context, machine_model MachineModel, model_options ModelOptions, Program P)
      : Context_(Context), MachineModel_(MachineModel), ModelOptions_(ModelOptions), Program_(P){};

  void compileProgram(std::string SourceFile);
  void compileProgram(std::string SourceFile, std::string ScopFunction);

  // initialize the cache model
  void initModel(std::vector<NamedLong> Parameters);
  void initModel();

  // prepare the between maps
  void computeBetweenAndFirstMaps();

  // count the actual performance info
  std::vector<NamedMisses> countCacheMisses();
  std::vector<NamedVector> countCacheMisses(std::vector<long> CacheSizes);

  std::vector<Access> getAccesses() const { return Accesses_; }
  std::map<std::string, std::vector<std::string>> getConflicts() const { return Conflicts_; }

  std::pair<unsigned, unsigned> getScopLoc() const { return Program_.getScopLoc(); }
  std::map<std::string, std::vector<access_info>> getAccessInfos() const { return Program_.getAccessInfos(); }

private:
  void computeGlobalMaps();
  void extractAccesses();
  void addConflicts(isl::union_map BetweenMap);

  // parameters
  isl::ctx Context_;
  machine_model MachineModel_;
  model_options ModelOptions_;

  // program information
  Program Program_;

  // analysis results
  isl::set Parameters_;
  std::vector<NamedLong> ParameterValues_;
  isl::union_map Schedule_;
  isl::map LexSuccEq_;
  isl::union_map SameLineSucc_;
  isl::union_map Before_;
  isl::union_map Forward_;
  isl::union_map First_;

  std::vector<Access> Accesses_;
  std::map<std::string, std::vector<std::string>> Conflicts_;
};

#endif