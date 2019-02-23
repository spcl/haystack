#ifndef _SCOP_H_
#define _SCOP_H_

#include "pet.h"
#include <isl/isl-noexceptions.h>
#include <map>
#include <string>

#include "Definitions.h"

// class parsing the source file and extracting schedule and memory accesses
class Program {
public:
  Program() = delete;
  Program(isl::ctx Context, machine_model MachineModel) : Context_(Context), MachineModel_(MachineModel) {}

  void setScop(isl::schedule Schedule, isl::union_map Reads, isl::union_map Writes,
               std::map<std::string, isl::set> ArrayExtents, std::map<std::string, long> ElementSizes);
  void extractScop(std::string SourceFile);

  void computeAccessToLine(isl::set Parameters);

  isl::union_map getSchedule() const { return ScheduleMap_; }
  isl::union_map getReads() const { return Reads_; }
  isl::union_map getWrites() const { return Writes_; }
  isl::union_set getAccessDomain() const { return AccessDomain_; }
  isl::union_map getAccessToLine() const { return AccessToLine_; }
  isl::union_map getAccessToElement() const { return AccessToElement_; }
  std::map<std::string, long> getElementSizes() const { return ElementSizes_; }

  size_t getNumOfReadReferences(std::string Statement) { return ReadReferences_[Statement].size(); }
  size_t getNumOfWriteReferences(std::string Statement) { return WriteReferences_[Statement].size(); }

private:
  void processScop();
  isl::map introduceCacheLines(std::string Name, isl::map AccessToArray, long ElementsPerCacheLine) const;
  isl::map introduceCacheSets(std::string Name, isl::map AccessToArray, long NumberOfCacheSets) const;

  isl::ctx Context_;
  machine_model MachineModel_;

  isl::schedule Schedule_;
  isl::union_map ScheduleMap_;
  isl::union_map Reads_;
  isl::union_map Writes_;
  isl::union_map AccessToLine_;
  isl::union_map AccessToElement_;
  isl::union_set AccessDomain_;

  // array information
  std::map<std::string, isl::set> ArrayExtents_;
  std::map<std::string, long> ElementSizes_;

  // read and write references per statement
  std::map<std::string, std::vector<std::string>> ReadReferences_;
  std::map<std::string, std::vector<std::string>> WriteReferences_;
};

#endif