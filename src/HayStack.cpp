#include <cassert>
#include <algorithm>
#include <limits>

#include "HayStack.h"
#include "Timer.h"
#include "isl-helpers.h"

#include "barvinok/isl.h"

isl_ctx *allocateContextWithIncludePaths(std::vector<std::string> IncludePaths) {
  // pass the include paths to the context
  std::vector<char *> Arguments;
  char Argument1[] = "program";
  char ArgumentI[] = "-I";
  Arguments.push_back(Argument1);
  for (auto &IncludePath : IncludePaths) {
    Arguments.push_back(ArgumentI);
    Arguments.push_back(const_cast<char *>(IncludePath.c_str()));
  }
  int ArgumentCount = Arguments.size();
  struct pet_options *options;
  options = pet_options_new_with_defaults();
  ArgumentCount = pet_options_parse(options, ArgumentCount, &Arguments[0], ISL_ARG_ALL);
  return isl_ctx_alloc_with_options(&pet_options_args, options);
}

void HayStack::compileProgram(std::string SourceFile) { Program_.extractScop(SourceFile); }

void HayStack::initModel(std::vector<NamedLong> Parameters) {
  assert(!Program_.getSchedule().is_null());
  // extract the parameters
  Parameters_ = Program_.getSchedule().params();
  int NumberOfParameters = Parameters_.dim(isl::dim::param);
  if (Parameters.size() < NumberOfParameters) {
    printf("\n\n-> exit(-1) not enough parameters\n");
    exit(-1);
  }
  std::map<int, NamedLong> SortedParameters;
  for (auto &Parameter : Parameters) {
    int Position = Parameters_.find_dim_by_name(isl::dim::param, Parameter.first);
    if (Position < 0 || Position >= NumberOfParameters) {
      printf("\n\n-> exit(-1) cannot find parameter %s\n", Parameter.first.c_str());
      exit(-1);
    }
    Parameters_ = Parameters_.fix_si(isl::dim::param, Position, Parameter.second);
    SortedParameters[Position] = Parameter;
  }
  // copy the sorted parameters to the parameter value vector
  for (auto &Entry : SortedParameters) {
    ParameterValues_.push_back(Entry.second);
  }
  // initialize the model
  initModel();
}

void HayStack::initModel() {
  // compute the access map
  Program_.computeAccessToLine(Parameters_);
  // compute the between map for all statements
  computeGlobalMaps();
  // extract the access information
  extractAccesses();
}

std::vector<NamedMisses> HayStack::countCacheMisses() {
  // define the result vector
  std::vector<NamedMisses> Results;
  // count the compulsory misses
  for (auto &Current : Accesses_) {
    Current.initAccess(ParameterValues_, Parameters_);
    Current.countCompulsoryMisses(First_);
  }
#ifdef DIMENSION_COUNTING
  // compute the capacity misses dimension by dimension
  isl::space Space = LexSuccEq_.domain().get_space();
  isl::map Universe = isl::map::universe(LexSuccEq_.get_space());
  auto Remaining = SameLineSucc_.reverse();
  for (int i = Space.dim(isl::dim::set) - 1; i >= 0; --i) {
    Timer::startTimer("ComputeBetweenMap");
    printf("-> processing dimension %d\n", i);
    auto Start = std::chrono::high_resolution_clock::now();
    // compute the filter for the dimension
    isl::local_space LSI = isl::local_space(Universe.domain().get_space());
    isl::local_space LSO = isl::local_space(Universe.range().get_space());
    isl::map Filter = Universe;
    for (int j = 0; j < i; ++j) {
      isl::pw_aff VarIn = isl::pw_aff::var_on_domain(LSI, isl::dim::set, j);
      isl::pw_aff VarOut = isl::pw_aff::var_on_domain(LSO, isl::dim::set, j);
      isl::map Constraint = VarIn.eq_map(VarOut);
      Filter = Filter.intersect(Constraint);
    }
    // compute the next map for the level
    auto Next = Remaining.intersect(Filter).lexmax();
    if (i > 0) {
      Remaining = Remaining.subtract_domain(Next.domain());
      // Remaining = Remaining.coalesce();
    }
    // compute the between map
    Next = Schedule_.apply_range(Next).apply_range(Schedule_.reverse()).coalesce();
    auto After = Next.apply_range(Forward_);
    auto BetweenMap = After.intersect(Before_);
    addConflicts(BetweenMap);
    BetweenMap = BetweenMap.apply_range(Program_.getAccessToLine());
    BetweenMap = BetweenMap.detect_equalities(); // important
    Timer::stopTimer("ComputeBetweenMap");
    // compute the cache misses
    for (auto &Current : Accesses_) {
      Current.computeStackDistances(BetweenMap);
    }
    auto Stop = std::chrono::high_resolution_clock::now();
    double Total = std::chrono::duration<double, std::milli>(Stop - Start).count();
    printf("-> done (%.2fms)\n", Total);
  }
#else
  Timer::startTimer("ComputeBetweenMap");
  auto Next = SameLineSucc_.reverse().lexmax();
  Next = Schedule_.apply_range(Next).apply_range(Schedule_.reverse()).coalesce();
  auto After = Next.apply_range(Forward_);
  auto BetweenMap = After.intersect(Before_);
  addConflicts(BetweenMap);
  BetweenMap = BetweenMap.apply_range(Program_.getAccessToLine());
  BetweenMap = BetweenMap.detect_equalities(); // important
  Timer::stopTimer("ComputeBetweenMap");
  // compute the cache misses
  for (auto &Current : Accesses_) {
    Current.computeStackDistances(BetweenMap);
  }
#endif
  // count the capacity misses and collect the results
  for (auto &Current : Accesses_) {
    Current.countCapacityMisses();
    auto CacheMisses = Current.getResult();
    Results.push_back(std::make_pair(Current.getName(), CacheMisses));
  }
#ifdef PREFETCHING
  // string based matching
  std::map<std::vector<int>, std::vector<int>> Streams;
  std::vector<int> Zeros(MachineModel_.CacheSizes.size(), 0);
  for (auto &Result : Results) {
    for (int i = 0; i < MachineModel_.CacheSizes.size(); ++i) {
      if (Result.second.PrefetchInfo.UnitStride && Result.second.CapacityMisses[i]) {
        auto Depth = Result.second.PrefetchInfo.PrefetchDepth;
        while (Depth.size() > 0) {
          if (Streams[Depth].size() == 0)
            Streams[Depth] = Zeros;
          Streams[Depth][i]++;
          Depth.pop_back();
        }
        // set the prefetch flag
        Result.second.PrefetchInfo.Prefetched[i] = true;
      }
    }
  }
  // copy stream info to the prefetch info
  for (auto &Result : Results) {
    auto Depth = Result.second.PrefetchInfo.PrefetchDepth;
    if (Streams.count(Depth) > 0)
      Result.second.PrefetchInfo.PrefetchStreams = Streams[Depth];
  }
#endif
  return Results;
}

std::vector<NamedVector> HayStack::countCacheMisses(std::vector<long> CacheSizes) {
  std::vector<NamedVector> Results;
  // compute the misses for all accesses
  for (auto &Current : Accesses_) {
    auto CacheMisses = Current.countCapacityMisses(CacheSizes);
    Results.push_back(std::make_pair(Current.getName(), CacheMisses));
  }
  return Results;
}

void HayStack::computeGlobalMaps() {
  Timer::startTimer("ComputeBetweenMap");
  // get the schedule and limit parameters
  Schedule_ = Program_.getSchedule();
  if (Parameters_.is_null()) {
    // if parameters are unknown introduce lower bound
    isl::set Parameters = Schedule_.params();
    int NumberOfParameters = Parameters.dim(isl::dim::param);
    for (int i = 0; i < NumberOfParameters; ++i) {
      Parameters = Parameters.lower_bound_si(isl::dim::param, 0, 0);
    }
    Schedule_ = Schedule_.intersect_domain(Parameters);
  } else {
    // otherwise fix the parameters
    Schedule_ = Schedule_.intersect_domain(Parameters_);
  }
  // filter statements without array accesses
  Schedule_ = Schedule_.intersect_domain(Program_.getAccessDomain());
  Schedule_ = Schedule_.coalesce();

  // map the iterations to the cache lines and cache sets
  isl::union_map IterToLine = Program_.getAccessToLine().apply_domain(Schedule_);
  IterToLine = IterToLine.coalesce();

  // get the successor maps
  isl::space Space = isl::set(IterToLine.domain()).get_space();
  isl::map LexSucc = isl::map::lex_lt(Space);
  LexSuccEq_ = isl::map::lex_le(Space);

  // compute accesses of the same cache line
  isl::union_map SameLine = IterToLine.apply_range(IterToLine.reverse());
  SameLineSucc_ = SameLine.intersect(LexSucc);
  SameLineSucc_ = SameLineSucc_.coalesce();

  // compute the before and forward maps
  Before_ = Schedule_.apply_range(LexSuccEq_.reverse()).apply_range(Schedule_.reverse()).coalesce();
  Forward_ = Schedule_.apply_range(LexSuccEq_).apply_range(Schedule_.reverse()).coalesce();
  Timer::stopTimer("ComputeBetweenMap");

  // compute the first map that connects the every memory location to the schedule value that loads it first
  First_ = Program_.getAccessToLine().reverse().apply_range(Schedule_).lexmin();
  First_ = Schedule_.apply_range(First_.reverse());
  First_ = First_.coalesce();
}

void HayStack::extractAccesses() {
  // extract the between map accesses
  auto extractAccess = [&](isl::set Set) {
    std::string Statement = Set.get_tuple_name();
    // compute the reads and writes
    int Reads = Program_.getNumOfReadReferences(Statement);
    int Writes = Program_.getNumOfWriteReferences(Statement);
    // iterate the read and write indexes
    for (int i = 0; i < Reads + Writes; ++i) {
      // set the name
      std::string Name = Statement;
      if (i < Reads)
        Name += "(R" + std::to_string(i) + ")";
      else
        Name += "(W" + std::to_string(i - Reads) + ")";
      // compute the domain
      isl::set Domain = Schedule_.domain().extract_set(Set.get_space());
      Domain = Domain.fix_si(isl::dim::set, Domain.dim(isl::dim::set) - 1, i);
#ifdef PREFETCHING
      Timer::startTimer("ComputePrefetchInfo");
      // compute the array elements accessed by the next statement instance
      isl::union_map Schedule = Schedule_.intersect_domain(Domain);
      isl::map LexSucc = isl::map::lex_lt(isl::set(Schedule.range()).get_space());
      isl::union_map Succ = Schedule.apply_range(isl::union_map(LexSucc).intersect_range(Schedule.range()));
      Succ = Succ.lexmin().apply_range(Schedule.reverse());
      Succ = Succ.apply_range(Program_.getAccessToElement());
      // compute the delta for all sub domains
      int MaxDepth = 0;
      isl::set MaxSet;
      bool UnitStride = false;
      auto searchMaps = [&](isl::map Succ) {
        auto searchBasicMaps = [&](isl::basic_map Succ) {
          // get the array name
          std::string Name = Succ.range().get_tuple_name();
          int ElementsPerCacheLine = Program_.getElementSizes()[Name];
          // compute the last variable dimension
          isl::union_set Depth = Schedule.reverse().apply_range(Succ).domain();
          int LastDepth = 0;
          isl::set LastSet;
          auto computeDepth = [&](isl::set Depth) {
            // drop project out dimensions until the set shrinks
            long Cardinality = isl::cardinality(Depth);
            int Dimension = 0;
            isl::set Set;
            do {
              Set = Depth.project_out(isl::dim::set, Dimension, Depth.dim(isl::dim::set) - Dimension);
              if (isl::cardinality(Set) == Cardinality)
                break;
              Dimension++;
            } while (Dimension < Depth.dim(isl::dim::set));
            if (LastDepth < Dimension) {
              LastDepth = Dimension;
              LastSet = Set;
            }
            return isl::stat::ok();
          };
          Depth.foreach_set(computeDepth);
          // if the last dimension is larger than for the previous successor cases
          if (LastDepth > MaxDepth) {
            // compute the delta
            isl::map AccessToElement = Program_.getAccessToElement().extract_map(Succ.get_space());
            isl::set Delta = AccessToElement.reverse().apply_range(Succ).deltas();
            // compute the prefetchable strides
            isl::set Close = isl::set::universe(Delta.get_space());
            isl::set Same = isl::set::universe(Delta.get_space());
            for (int i = 0; i < Delta.dim(isl::dim::set) - 1; ++i) {
              Close = Close.fix_si(isl::dim::set, i, 0);
              Same = Same.fix_si(isl::dim::set, i, 0);
            }
            Close = Close.upper_bound_si(isl::dim::set, Delta.dim(isl::dim::set) - 1, ElementsPerCacheLine);
            Close = Close.lower_bound_si(isl::dim::set, Delta.dim(isl::dim::set) - 1, -ElementsPerCacheLine);
            Same = Same.fix_si(isl::dim::set, Delta.dim(isl::dim::set) - 1, 0);
            // set prefetchable if not the same is accessed
            if (Close.intersect(Delta).is_empty()) {
              MaxDepth = LastDepth;
              UnitStride = false;
            } else if (Same.intersect(Delta).is_empty()) {
              UnitStride = true;
              MaxDepth = LastDepth;
              MaxSet = LastSet;
            }
          }
          return isl::stat::ok();
        };
        Succ.foreach_basic_map(searchBasicMaps);
        return isl::stat::ok();
      };
      Succ.foreach_map(searchMaps);
      // initialize the prefetch info
      std::vector<int> PrefetchDepth;
      if (!MaxSet.is_null()) {
        for (int i = 0; i < MaxDepth; ++i) {
          long Minimum = isl::computeMinimum(MaxSet, i);
          long Maximum = isl::computeMaximum(MaxSet, i);
          if (Minimum == Maximum) {
            PrefetchDepth.push_back(Minimum);
          }
        }
      }
      std::vector<int> Zeros(MachineModel_.CacheSizes.size(), 0);
      std::vector<bool> NotTrue(MachineModel_.CacheSizes.size(), false);
      prefetch_info Prefetched = {UnitStride, PrefetchDepth, Zeros, NotTrue};
      Timer::stopTimer("ComputePrefetchInfo");
#else
      prefetch_info Prefetched = {false, {}, {}, {}};
#endif
      // create the access
      Access Current(Name, MachineModel_, Domain, Program_.getElementSizes(), Prefetched);
      Accesses_.push_back(Current);
    }
    return isl::stat::ok();
  };
  Schedule_.domain().foreach_set(extractAccess);
  std::sort(Accesses_.begin(), Accesses_.end());
}

void HayStack::addConflicts(isl::union_map BetweenMap) {
#ifdef COMPUTE_CONFLICTS
  // count the statements
  for (auto &Source : Accesses_) {
    auto Destinations = BetweenMap.intersect_domain(Source.getDomain());
    if (!Destinations.is_empty()) {
      for (auto &Destination : Accesses_) {
        if (!Destinations.intersect_range(Destination.getDomain()).is_empty()) {
          Conflicts_[Source.getName()].push_back(Destination.getName());
        }
      }
      std::sort(Conflicts_[Source.getName()].begin(), Conflicts_[Source.getName()].end());
      auto Last = std::unique(Conflicts_[Source.getName()].begin(), Conflicts_[Source.getName()].end());
      Conflicts_[Source.getName()].erase(Last, Conflicts_[Source.getName()].end());
    }
  }
#endif
}
