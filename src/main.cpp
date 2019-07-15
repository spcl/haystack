/*
 * Copyright (c) 2019, ETH Zurich
 */

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <numeric>
#include <string>
#include <vector>

#include <isl/options.h>

#include "HayStack.h"
#include "Timer.h"

namespace po = boost::program_options;

const int CACHE_LINE_SIZE = 64;

// default
const int CACHE_SIZE1 = 32 * 1024;
const int CACHE_SIZE2 = 512 * 1024;

bool check_path(std::string path) {
  std::ifstream f(path.c_str());
  return f.good();
}

std::map<int, std::string> compute_lines(std::string FileName, std::pair<long, long> ScopLoc) {
  std::map<int, std::string> Result;
  std::ifstream SourceFile;
  SourceFile.open(FileName);
  std::string Line;
  int LineNumber = 0;
  while (std::getline(SourceFile, Line)) {
    LineNumber++;
    if (SourceFile.tellg() > ScopLoc.first && SourceFile.tellg() <= ScopLoc.second) {
      Result[LineNumber] = Line;
    }
  }
  SourceFile.close();
  return Result;
}

std::map<int, std::pair<long, long>> compute_offsets(std::string FileName) {
  std::map<int, std::pair<long, long>> Result;
  std::ifstream SourceFile;
  SourceFile.open(FileName);
  std::string Line;
  int LineNumber = 0;
  long Start = SourceFile.tellg();
  long End = 0;
  while (std::getline(SourceFile, Line)) {
    End = SourceFile.tellg();
    Result[++LineNumber] = std::make_pair(Start, End);
    Start = End;
  }
  SourceFile.close();
  return Result;
}

void print_scop(std::map<int, std::string> &Lines, int Start, int Stop) {
  // compute number of necessary digits
  int Width = std::to_string(Lines.end()->first).length();
  for (int i = Start; i < Stop; ++i) {
    std::cout << std::setw(Width) << std::right << std::to_string(i);
    std::cout << " " << Lines[i] << std::endl;
  }
}

// define print operators
namespace std {
std::ostream &operator<<(std::ostream &os, const std::vector<long> &vec) {
  for (int i = 0; i < vec.size(); ++i) {
    os << vec[i];
    if (i < vec.size() - 1)
      os << " ";
  }
  return os;
}
} // namespace std

void run_model(isl::ctx Context, po::variables_map Variables) {
  // allocate the machine model with default values
  machine_model MachineModel = {Variables["line-size"].as<long>(), Variables["cache-sizes"].as<std::vector<long>>()};
  model_options ModelOptions = {Variables["compute-bounds"].as<bool>()};
  printf("-> setting up cache levels\n");
  std::sort(MachineModel.CacheSizes.begin(), MachineModel.CacheSizes.end());
  for (auto CacheSize : MachineModel.CacheSizes) {
    if (CacheSize % 1024 == 0) {
      printf("   - %ldkB with %ldB cache lines\n", CacheSize / 1024, MachineModel.CacheLineSize);
    } else {
      printf("   - %ldB with %ldB cache lines\n", CacheSize, MachineModel.CacheLineSize);
    }
  }
  printf("-> done\n");
  // compute the total time
  auto StartExecution = std::chrono::high_resolution_clock::now();
  // allocate the cache model and compile the program
  HayStack Model(Context, MachineModel, ModelOptions);
  if (Variables.count("scop-function") == 0) {
    Model.compileProgram(Variables["input-file"].as<std::string>());
  } else {
    Model.compileProgram(Variables["input-file"].as<std::string>(), Variables["scop-function"].as<std::string>());
  }
  // parsing the parameters
  std::vector<NamedLong> Parameters;
  if (Variables.count("define-parameters") != 0) {
    printf("-> parsing the parameters...\n");
    std::vector<std::string> ParamStrings = Variables["define-parameters"].as<std::vector<std::string>>();
    for (auto ParamString : ParamStrings) {
      std::string Name;
      long Value;
      try {
        std::string Delimiter = "=";
        auto Split = ParamString.find(Delimiter);
        if (Split == std::string::npos) {
          throw std::runtime_error("did not find delimiter");
        }
        std::string Name = ParamString.substr(0, ParamString.find(Delimiter));
        long Value = std::stol(ParamString.substr(ParamString.find(Delimiter) + 1));
        printf("   - %s = %ld\n", Name.c_str(), Value);
        Parameters.push_back(std::make_pair(Name, Value));
      } catch (std::exception) {
        printf("-> exit(-1) failed to parse %s\n", ParamString.c_str());
        exit(-1);
      }
    }
    printf("-> done\n");
  }
  // run the preprocessing
  printf("-> start processing...\n");
  auto Start = std::chrono::high_resolution_clock::now();
  Model.initModel(Parameters);
  // execute the cache model
  auto CacheMisses = Model.countCacheMisses();
  auto Stop = std::chrono::high_resolution_clock::now();
  double TotalEvaluation = std::chrono::duration<double, std::milli>(Stop - Start).count();
  printf("-> done after (%.2fms)\n", TotalEvaluation);
  // collect and print result
  long TotalAccesses = 0;
  long TotalCompulsory = 0;
  std::vector<long> TotalCapacity(MachineModel.CacheSizes.size(), 0);
  // sum the cache misses for all accesses
  for (auto &CacheMiss : CacheMisses) {
    TotalAccesses += CacheMiss.second.Total;
    TotalCompulsory += CacheMiss.second.CompulsoryMisses;
    std::transform(TotalCapacity.begin(), TotalCapacity.end(), CacheMiss.second.CapacityMisses.begin(),
                   TotalCapacity.begin(), std::plus<long>());
  };
  // open the input file and seek the start of the scop
  std::map<int, std::string> Lines = compute_lines(Variables["input-file"].as<std::string>(), Model.getScopLoc());
  std::map<int, std::pair<long, long>> Offsets = compute_offsets(Variables["input-file"].as<std::string>());
  long Position = Lines.begin()->first;
  std::string LineStart;
  std::string SingleLine;
  std::string DoubleLine;
  LineStart.resize(std::to_string(Lines.rbegin()->first).length() + 1, ' ');
  SingleLine.resize(80 - LineStart.length(), '-');
  DoubleLine.resize(80, '=');
  // print the access infos sorted by position
  size_t RefWidth = 16;
  std::map<long, std::vector<access_info>> AccessInfosByLn;
  std::map<std::string, access_info> AccessInfoByName;
  for (auto AccessInfos : Model.getAccessInfos()) {
    if (AccessInfos.second.empty())
      continue;
    AccessInfosByLn[AccessInfos.second[0].Line] = AccessInfos.second;
    for (auto AccessInfo : AccessInfos.second) {
      AccessInfoByName[AccessInfo.Name] = AccessInfo;
      RefWidth = std::max(RefWidth, AccessInfo.Name.length() + 1);
    }
  }
  // print the cache info access by access
  std::cout << DoubleLine << std::endl;
  std::cout << "                  relative number of cache misses (statement)" << std::endl;
  std::cout << DoubleLine << std::endl;
  for (auto AccessInfos : AccessInfosByLn) {
    // determine the last line of multiline
    int Next = AccessInfos.first;
    while (Offsets[Next].second < AccessInfos.second[0].Stop) {
      Next++;
    }
    // print the sources
    print_scop(Lines, Position, Next + 1);
    Position = Next + 1;
    // print header
    std::cout << LineStart << SingleLine << std::endl;
    std::cout << std::setw(RefWidth) << std::right << "ref";
    std::cout << "  ";
    std::cout << std::setw(6) << std::left << "type";
    std::cout << std::setw(10) << std::left << "comp[%]";
    for (int i = 1; i <= MachineModel.CacheSizes.size(); ++i) {
      std::string Capacity = "L" + std::to_string(i) + "[%]";
      std::cout << std::setw(10) << std::left << Capacity;
    }
    std::cout << std::setw(10) << std::left << "tot[%]";
    std::cout << std::setw(10) << std::left << "reuse[ln]";
    std::cout << std::endl;
    // print the accesses
    for (auto AccessInfo : AccessInfos.second) {
      // find the actual cache miss info
      auto Iter = std::find_if(CacheMisses.begin(), CacheMisses.end(),
                               [&](NamedMisses Misses) { return Misses.first == AccessInfo.Name; });
      assert(Iter != CacheMisses.end());
      auto Compulsory = Iter->second.CompulsoryMisses;
      auto Capacity = Iter->second.CapacityMisses;
      auto Total = Iter->second.Total;
      // print the access info
      std::string Name = AccessInfo.Access;
      if (Name.length() > RefWidth)
        Name = Name.substr(0, RefWidth);
      std::cout << std::setw(RefWidth) << std::right << Name;
      std::cout << "  ";
      std::cout << std::setw(6) << std::left << (AccessInfo.ReadOrWrite == Read ? "rd" : "wr");
      std::cout << std::setw(10) << std::left << std::setprecision(5) << std::fixed
                << 100.0 * (double)Compulsory / (double)TotalAccesses;
      for (int i = 0; i < MachineModel.CacheSizes.size(); ++i) {
        std::cout << std::setw(10) << std::left << std::setprecision(5) << std::fixed
                  << 100.0 * (double)Capacity[i] / (double)TotalAccesses;
      }
      std::cout << std::setw(10) << std::left << std::setprecision(5) << std::fixed
                << 100.0 * (double)Total / (double)TotalAccesses;
      // compute the reuse line numbers
      auto Conflicts = Model.getConflicts()[AccessInfo.Name];
      // compute the reuse line numbers
      std::vector<int> ReuseLines;
      for (auto Conflict : Conflicts) {
        ReuseLines.push_back(AccessInfoByName[Conflict].Line);
      }
      // sort the line numbers and remove duplicates
      std::sort(ReuseLines.begin(), ReuseLines.end());
      auto Last = std::unique(ReuseLines.begin(), ReuseLines.end());
      for (auto Iter = ReuseLines.begin(); Iter != Last;) {
        std::cout << *Iter;
        if (++Iter != Last)
          std::cout << ",";
      }
      std::cout << std::endl;
    }
    std::cout << LineStart << SingleLine << std::endl;
  }
  print_scop(Lines, Position, Lines.rbegin()->first + 1);
  // print the scop info
  std::cout << DoubleLine << std::endl;
  std::cout << "                     absolute number of cache misses (SCOP)" << std::endl;
  std::cout << DoubleLine << std::endl;
  std::cout.imbue(std::locale(""));
  std::cout << std::setw(16) << std::left << "compulsory:";
  std::cout << std::setw(20) << std::right << TotalCompulsory << std::endl;
  for (int i = 1; i <= MachineModel.CacheSizes.size(); ++i) {
    std::string Capacity = "capacity (L" + std::to_string(i) + "):";
    std::cout << std::setw(16) << std::left << Capacity;
    std::cout << std::setw(20) << std::right << TotalCapacity[i - 1] << std::endl;
  }
  std::cout << std::setw(16) << std::left << "total:";
  std::cout << std::setw(20) << std::right << TotalAccesses << std::endl;
  std::cout << DoubleLine << std::endl;
}

int main(int argc, const char **args) {
  try {
    // define the program options
    po::options_description Descriptor("Program options");
    Descriptor.add_options()                    //
        ("help,h", "print the program options") //
        ("cache-sizes,c", po::value<std::vector<long>>()->multitoken()->default_value({CACHE_SIZE1, CACHE_SIZE2}),
         "cache sizes in byte")                                                                                       //
        ("line-size,l", po::value<long>()->default_value(CACHE_LINE_SIZE), "cache-line size in byte")                 //
        ("input-file,f", po::value<std::string>(), "set the source file [file name]")                                 //
        ("include-path,I", po::value<std::vector<std::string>>(), "set the include path [include path]")              //
        ("define-parameters,d", po::value<std::vector<std::string>>()->multitoken(), "parameter values [N=10 M=100]") //
        ("scop-function,s", po::value<std::string>(), "set the scop function scop")                                   //
        ("compute-bounds,b", po::value<bool>()->default_value(false), "compute stack distance bounds");

    // parse the program options
    po::variables_map Variables;
    po::store(po::parse_command_line(argc, args, Descriptor), Variables);
    po::notify(Variables);
    if (Variables.count("help") || Variables.count("input-file") == 0) {
      std::cout << Descriptor << std::endl;
      return 0;
    }
    // check if the include paths are valid
    for (int i = 0; i < Variables.count("include-path"); ++i) {
      std::string IncludePath = Variables["include-path"].as<std::vector<std::string>>()[i];
      if (!check_path(IncludePath)) {
        printf("-> exit(-1) include path %s not valid\n", IncludePath.c_str());
        exit(-1);
      }
    }
    // check if the source file is valid
    if (!check_path(Variables["input-file"].as<std::string>())) {
      printf("-> exit(-1) input file %s not found\n", Variables["input-file"].as<std::string>().c_str());
      exit(-1);
    }
    // allocate the context outside of the cache model
    std::vector<std::string> IncludePaths;
    if (Variables.count("include-path") > 0)
      IncludePaths = Variables["include-path"].as<std::vector<std::string>>();
    isl::ctx Context = allocateContextWithIncludePaths(IncludePaths);
    isl_options_set_on_error(Context.get(), ISL_ON_ERROR_ABORT);

    // run the cache model
    run_model(Context, Variables);

    isl_ctx_free(Context.get());
  } catch (const boost::program_options::error &ex) {
    printf("-> exit(-1) option parsing error: %s\n", ex.what());
  }
  return 0;
}
