
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <isl/options.h>

#include "HayStack.h"
#include "Timer.h"

namespace po = boost::program_options;

const int CACHE_LINE_SIZE = 64;

// Xeon Gold 6150
const int CACHE_SIZE1 = 32 * 1024; // 8-way
const int CACHE_SIZE2 = 1024 * 1024; // 16-way (non-inclusive which is close to inclusive)
const int CACHE_SIZE3 = 1408 * 1024; // 11-way (non-inclusive which is close to inclusive)

// polycache
// const int CACHE_SIZE1 = 32 * 1024;  // 4-way
// const int CACHE_SIZE2 = 256 * 1024; // 4-way

struct range {
  std::string Name;
  int Current;
  int Start;
  int Stop;
  int Increment;
};

int main(int argc, const char **args) {
  // define the program options
  po::options_description Descriptor("Program options");
  Descriptor.add_options()                                                                                 //
      ("help,h", "print the program options")                                                              //
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

  // allocate the context outside of the cache model
  std::vector<std::string> IncludePaths;
  if (Variables.count("include-path") > 0)
    IncludePaths = Variables["include-path"].as<std::vector<std::string>>();
  isl::ctx Context = allocateContextWithIncludePaths(IncludePaths);
  isl_options_set_on_error(Context.get(), ISL_ON_ERROR_ABORT);
  // TODO test if BV_CHAMBERS_ISL or BV_CHAMBERS_POLYLIB is more efficient

  {
    // allocate the machine model
    machine_model MachineModel = {CACHE_LINE_SIZE, {CACHE_SIZE1, CACHE_SIZE2}};
    // compute the total time
    auto StartExecution = std::chrono::high_resolution_clock::now();
    // allocate the cache model and compile the program
    HayStack Model(Context, MachineModel);
    Model.compileProgram(Variables["input-file"].as<std::string>());
    // run the preprocessing
    printf("-> start preprocessing...\n");
    auto StartPreprocessing = std::chrono::high_resolution_clock::now();
    Model.initModel();
    auto StopPreprocessing = std::chrono::high_resolution_clock::now();
    double TotalPreprocessing =
        std::chrono::duration<double, std::milli>(StopPreprocessing - StartPreprocessing).count();
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
    for (auto &CacheMiss : CacheMisses) {
      TotalAccesses += CacheMiss.second.Total;
      TotalCompulsory += CacheMiss.second.CompulsoryMisses;
      std::transform(TotalCapacity.begin(), TotalCapacity.end(), CacheMiss.second.CapacityMisses.begin(),
                     TotalCapacity.begin(), std::plus<long>());
      // print intermediate results if verbose is true
      if (Variables["verbose"].as<bool>()) {
        std::cout << std::fixed;
        std::cout << std::setprecision(2);
        std::cout << "   -> " << CacheMiss.first << " counted ";
        std::cout << CacheMiss.second.CompulsoryMisses << "/";
        for (int i = 0; i < MachineModel.CacheSizes.size(); ++i) {
          std::cout << CacheMiss.second.CapacityMisses[i] << "/";
        }
        std::cout << CacheMiss.second.Total << " (CO/";
        for (int i = 0; i < MachineModel.CacheSizes.size(); ++i)
          std::cout << "CA" << i << "/";
        std::cout << "TO) ";
#ifdef PREFETCHING
        std::cout << "using ";
        for (int i = 0; i < MachineModel.CacheSizes.size(); ++i) {
          int Streams = 0;
          if (CacheMiss.second.PrefetchInfo.Prefetched[i])
            Streams = CacheMiss.second.PrefetchInfo.PrefetchStreams[i];
          std::cout << Streams;
          if (i != MachineModel.CacheSizes.size() - 1)
            std::cout << "/";
        }
        std::cout << " prefetch streams ";
#endif
        std::cout << "with ";
        std::cout << CacheMiss.second.Methods.Barvinok << "/";
        std::cout << CacheMiss.second.Methods.Bernstein << "/";
        std::cout << CacheMiss.second.Methods.Enumerate << " (BA/BE/EN) after ";
        std::cout << CacheMiss.second.Splits.Barvinok << "/";
        std::cout << CacheMiss.second.Splits.Equalize << "/";
        std::cout << CacheMiss.second.Splits.Raster << "/";
        std::cout << CacheMiss.second.Splits.Enumerate << " (BA/EQ/RA/EN) splits\n";
      }
    }
    // print the summary
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    std::cout << "-> done (" << TotalEvaluation << "ms) accumulated ";
    std::cout << TotalCompulsory << "/";
    for (int i = 0; i < MachineModel.CacheSizes.size(); ++i)
      std::cout << TotalCapacity[i] << "/";
    std::cout << TotalAccesses << " (CO/";
    for (int i = 0; i < MachineModel.CacheSizes.size(); ++i)
      std::cout << "CA" << i << "/";
    std::cout << "TO) cache misses\n";
    // print timing information
    auto StopExecution = std::chrono::high_resolution_clock::now();
    double TotalExecution = std::chrono::duration<double, std::milli>(StopExecution - StartExecution).count();
    printf("-> finished after (%.2fms)\n", TotalExecution);
    Timer::printClocks();
    // print the affinity info
    std::map<std::vector<int>, int> Affinity;
    for (auto &CacheMiss : CacheMisses) {
      for (auto Aff : CacheMiss.second.Affinity) {
        Affinity[Aff.first] += Aff.second;
      }
    }
    if (Affinity.size() > 0) {
      std::cout << "==================================================" << std::endl;
      // print the affinity
      for (auto Aff : Affinity) {
        std::cout << " - NonAffine " << Aff.first[0] << " Affine " << Aff.first[1] << " : " << Aff.second << std::endl;
      }
      std::cout << "==================================================" << std::endl;
    }

    // print the conflicts
    auto Conflicts = Model.getConflicts();
    if (Conflicts.size() > 0) {
      std::cout << "==================================================" << std::endl;
      for (auto Conflict : Conflicts) {
        std::cout << " - " << Conflict.first << ": ";
        for (auto Name : Conflict.second) {
          std::cout << Name << " ";
        }
        std::cout << std::endl;
      }
      std::cout << "==================================================" << std::endl;
    }

    // // compute cache miss curve
    // std::cout << "==================================================" << std::endl;
    // std::vector<long> CacheSizes;
    // for (int i = 1 * 1024; i <= CACHE_SIZE2; i *= 2)
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
    // std::cout << "==================================================" << std::endl;
  }
  isl_ctx_free(Context.get());

  return 0;
}
