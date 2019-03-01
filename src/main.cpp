
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
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

void print_scop(std::ifstream &SourceFile, int Length) {
  char *Buffer = new char[Length + 1];
  Buffer[Length] = 0;
  SourceFile.read(Buffer, Length);
  printf("%s", &Buffer[0]);
  delete Buffer;
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
  printf("-> setting up cache levels\n");
  std::sort(MachineModel.CacheSizes.begin(), MachineModel.CacheSizes.end());
  for (auto CacheSize : MachineModel.CacheSizes) {
    if (CacheSize % 1024 == 0) {
      printf("   - %ld kB with %ld B cache lines\n", CacheSize / 1024, MachineModel.CacheLineSize);
    } else {
      printf("   - %ld B with %ld B cache lines\n", CacheSize, MachineModel.CacheLineSize);
    }
  }
  printf("-> done\n");
  // compute the total time
  auto StartExecution = std::chrono::high_resolution_clock::now();
  printf("-> searching scop\n");
  // allocate the cache model and compile the program
  HayStack Model(Context, MachineModel);
  Model.compileProgram(Variables["input-file"].as<std::string>());
  printf("-> done\n");
  // run the preprocessing
  printf("-> start preprocessing...\n");
  auto StartPreprocessing = std::chrono::high_resolution_clock::now();
  Model.initModel();
  auto StopPreprocessing = std::chrono::high_resolution_clock::now();
  double TotalPreprocessing = std::chrono::duration<double, std::milli>(StopPreprocessing - StartPreprocessing).count();
  printf("-> done (%.2fms)\n", TotalPreprocessing);
  // execute the cache model
  auto StartEvaluation = std::chrono::high_resolution_clock::now();
  auto CacheMisses = Model.countCacheMisses();
  auto StopEvaluation = std::chrono::high_resolution_clock::now();
  double TotalEvaluation = std::chrono::duration<double, std::milli>(StopEvaluation - StartEvaluation).count();
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
  auto ScopLoc = Model.getScopLoc();
  int Pos = ScopLoc.first;
  std::ifstream SourceFile;
  SourceFile.open(Variables["input-file"].as<std::string>());
  SourceFile.seekg(Pos, std::ios::beg);
  // print the access infos sorted by position
  std::map<long, std::vector<access_info>> Ordered;
  for (auto AccessInfos : Model.getAccessInfos()) {
    if(AccessInfos.second.empty())
      continue;
    Ordered[AccessInfos.second[0].Stop] = AccessInfos.second;
  }
  // print the cache info access by access
  for (auto AccessInfos : Ordered) {
    // print the sources
    print_scop(SourceFile, AccessInfos.first - Pos);
    Pos = AccessInfos.first;
    // print header
    const int Width = 16;
    printf("-------------------------------------------------------------------------------\n");
    std::cout << std::setw(Width) << std::left << "acc";
    std::cout << std::setw(Width/2) << std::left << "rd/wr";
    std::cout << std::setw(Width) << std::left << "comp[%]";
    for (int i = 1; i <= MachineModel.CacheSizes.size(); ++i) {
      std::string Capacity = "L" + std::to_string(i) + "[%]";
      std::cout << std::setw(Width) << std::left << Capacity;
    }
    std::cout << std::setw(Width) << std::left << "tot[%]";
    std::cout << std::endl;
    // print the accesses
    for (auto AccessInfo : AccessInfos.second) {
      std::cout << std::setw(Width) << std::left << AccessInfo.Access;
      std::cout << std::setw(Width/2) << std::left << (AccessInfo.ReadOrWrite == Read ? "rd" : "wr");
      // find the actual cache miss info
      auto Iter = std::find_if(CacheMisses.begin(), CacheMisses.end(),
                               [&](NamedMisses Misses) { return Misses.first == AccessInfo.Name; });
      assert(Iter != CacheMisses.end());
      auto Compulsory = Iter->second.CompulsoryMisses;
      auto Capacity = Iter->second.CapacityMisses;
      auto Total = Iter->second.Total;
      std::cout << std::setw(Width) << std::left << std::setprecision(4) << std::fixed
                << 100.0 * (double)Compulsory / (double)TotalAccesses;
      for (int i = 0; i < MachineModel.CacheSizes.size(); ++i) {
        std::cout << std::setw(Width) << std::left << std::setprecision(4) << std::fixed
                  << 100.0 * (double)Capacity[i] / (double)TotalAccesses;
      }
      std::cout << std::setw(Width) << std::left << std::setprecision(4) << std::fixed
                << 100.0 * (double)Total / (double)TotalAccesses;
      std::cout << std::endl;
    }
    printf("-------------------------------------------------------------------------------\n");
  }
  print_scop(SourceFile, ScopLoc.second - Pos);
  SourceFile.close();
  // print the scop info
  printf("-------------------------------------------------------------------------------\n");
  printf(" - compulsory misses:\t\t%ld\n", TotalCompulsory);
  for (int i = 1; i <= MachineModel.CacheSizes.size(); ++i)
    printf(" - capacity misses (L%d):\t%ld\n", i, TotalCapacity[i - 1]);
  printf(" - memory accesses:\t\t%ld\n", TotalAccesses);
  printf("-------------------------------------------------------------------------------\n");

  //     // print intermediate results if verbose is true
  //     if (Variables["verbose"].as<bool>()) {
  //       std::cout << std::fixed;
  //       std::cout << std::setprecision(2);
  //       std::cout << "   -> " << CacheMiss.first << " counted ";
  //       std::cout << CacheMiss.second.CompulsoryMisses << "/";
  //       for (int i = 0; i < MachineModel.CacheSizes.size(); ++i) {
  //         std::cout << CacheMiss.second.CapacityMisses[i] << "/";
  //       }
  //       std::cout << CacheMiss.second.Total << " (CO/";
  //       for (int i = 0; i < MachineModel.CacheSizes.size(); ++i)
  //         std::cout << "CA" << i << "/";
  //       std::cout << "TO) ";
  // #ifdef PREFETCHING
  //       std::cout << "using ";
  //       for (int i = 0; i < MachineModel.CacheSizes.size(); ++i) {
  //         int Streams = 0;
  //         if (CacheMiss.second.PrefetchInfo.Prefetched[i])
  //           Streams = CacheMiss.second.PrefetchInfo.PrefetchStreams[i];
  //         std::cout << Streams;
  //         if (i != MachineModel.CacheSizes.size() - 1)
  //           std::cout << "/";
  //       }
  //       std::cout << " prefetch streams ";
  // #endif
  //       std::cout << "\n";
  //     }

  // // print the summary
  // std::cout << std::fixed;
  // std::cout << std::setprecision(2);
  // std::cout << "-> done (" << TotalEvaluation << "ms) accumulated ";
  // std::cout << TotalCompulsory << "/";
  // for (int i = 0; i < MachineModel.CacheSizes.size(); ++i)
  //   std::cout << TotalCapacity[i] << "/";
  // std::cout << TotalAccesses << " (CO/";
  // for (int i = 0; i < MachineModel.CacheSizes.size(); ++i)
  //   std::cout << "CA" << i << "/";
  // std::cout << "TO) cache misses\n";
  // // print timing information
  // auto StopExecution = std::chrono::high_resolution_clock::now();
  // double TotalExecution = std::chrono::duration<double, std::milli>(StopExecution - StartExecution).count();
  // printf("-> finished after (%.2fms)\n", TotalExecution);
  // // print the conflicts
  // auto Conflicts = Model.getConflicts();
  // if (Conflicts.size() > 0) {
  //   std::cout << "==================================================" << std::endl;
  //   for (auto Conflict : Conflicts) {
  //     std::cout << " - " << Conflict.first << ": ";
  //     for (auto Name : Conflict.second) {
  //       std::cout << Name << " ";
  //     }
  //     std::cout << std::endl;
  //   }
  //   std::cout << "==================================================" << std::endl;
  // }

  // // compute cache miss curve
  // std::cout << "==================================================" <<
  // std::endl; std::vector<long> CacheSizes; for (int i = 1 * 1024; i <=
  // CACHE_SIZE2; i *= 2)
  //   CacheSizes.push_back(i);
  // // compute the misses
  // auto Accesses = Model.countCacheMisses(CacheSizes);
  // for (int i = 0; i < CacheSizes.size(); ++i) {
  //   // compute total number of memory accesses
  //   long MissCount = TotalCompulsory;
  //   for (auto Access : Accesses) {
  //     MissCount += Access.second[i];
  //   }
  //   double MissRate = 100.0 * (double)MissCount / (double)TotalAccesses;
  //   std::cout << std::setprecision(1) << MissRate << " ";
  // }
  // std::cout << std::endl;
  // std::cout << "==================================================" <<
  // std::endl;
}

int main(int argc, const char **args) {
  try {
    // define the program options
    po::options_description Descriptor("Program options");
    Descriptor.add_options()                    //
        ("help,h", "print the program options") //
        ("cache-sizes,c", po::value<std::vector<long>>()->multitoken()->default_value({CACHE_SIZE1, CACHE_SIZE2}),
         "cache sizes in kilo byte")                                                                         //
        ("line-size,l", po::value<long>()->default_value(CACHE_LINE_SIZE), "cache-line size in byte")        //
        ("input-file,f", po::value<std::string>(), "specify the source file [file name]")                    //
        ("include-path,I", po::value<std::vector<std::string>>(), "specify the include path [include path]") //
        ("verbose,v", po::value<bool>()->default_value(false), "print additional information");

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
