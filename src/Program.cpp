#include <algorithm>

#include "Program.h"
#include "isl-helpers.h"
#include "op.h"

void Program::extractScop(std::string SourceFile) {
  pet_scop *PetScop = pet_scop_extract_from_C_source(Context_.get(), SourceFile.c_str(), NULL);
  if (PetScop == nullptr) {
    printf("-> exit(-1) cannot extract scope\n");
    exit(-1);
  }

  // get the schedule and access information
  // (the tagged accesses allow us to distinguish multiple accesses of the same array)
  Schedule_ = isl::manage(pet_scop_get_schedule(PetScop));
  Reads_ = isl::manage(pet_scop_get_tagged_may_reads(PetScop));
  Writes_ = isl::manage(pet_scop_get_tagged_may_writes(PetScop));

  // extract the array information
  for (int idx = 0; idx < PetScop->n_array; idx++) {
    isl::set Extent = isl::manage_copy(PetScop->arrays[idx]->extent);
    // ignore scalars
    if (Extent.dim(isl::dim::set) == 0)
      continue;
    // store the array information
    std::string Name = Extent.get_tuple_name();
    ArrayExtents_[Name] = Extent;
    ElementSizes_[Name] = PetScop->arrays[idx]->element_size;
  }

  pet_scop_free(PetScop);
  processScop();
}

void Program::setScop(isl::schedule Schedule, isl::union_map Reads, isl::union_map Writes,
                      std::map<std::string, isl::set> ArrayExtents, std::map<std::string, long> ElementSizes) {
  Schedule_ = Schedule;
  Reads_ = Reads;
  Writes_ = Writes;
  ArrayExtents_ = ArrayExtents;
  ElementSizes_ = ElementSizes;

  processScop();
}

void Program::processScop() {
  // extract the read and write references
  auto extractReads = [&](isl::map Map) {
    if (ElementSizes_.count(Map.get_tuple_name(isl::dim::out)) > 0) {
      std::string Statement = Map.domain().unwrap().get_tuple_name(isl::dim::in);
      std::string Reference = Map.domain().unwrap().get_tuple_name(isl::dim::out);
      ReadReferences_[Statement].push_back(Reference);
    }
    return isl::stat::ok();
  };
  Reads_.foreach_map(extractReads);
  auto extractWrites = [&](isl::map Map) {
    if (ElementSizes_.count(Map.get_tuple_name(isl::dim::out)) > 0) {
      std::string Statement = Map.domain().unwrap().get_tuple_name(isl::dim::in);
      std::string Reference = Map.domain().unwrap().get_tuple_name(isl::dim::out);
      WriteReferences_[Statement].push_back(Reference);
    }
    return isl::stat::ok();
  };
  Writes_.foreach_map(extractWrites);
  // sort the references
  for (auto &ReadReference : ReadReferences_)
    std::sort(ReadReference.second.begin(), ReadReference.second.end());
  for (auto &WriteReference : WriteReferences_)
    std::sort(WriteReference.second.begin(), WriteReference.second.end());

  // extend the schedule with the reference dimension
  ScheduleMap_ = Schedule_.get_map().intersect_domain(Schedule_.get_domain());
  isl::union_map ScheduleExt = isl::map::empty(ScheduleMap_.get_space());
  auto extendSchedule = [&](isl::map Schedule) {
    std::string Statement = Schedule.get_tuple_name(isl::dim::in);
    if (!(ReadReferences_[Statement].empty() && WriteReferences_[Statement].empty())) {
      // extend the schedule
      Schedule = Schedule.insert_dims(isl::dim::in, Schedule.dim(isl::dim::in), 1);
      Schedule = Schedule.insert_dims(isl::dim::out, Schedule.dim(isl::dim::out), 1);
      Schedule = Schedule.set_tuple_name(isl::dim::in, Statement);
      // connect the new dimension
      isl::local_space LSIn = isl::local_space(Schedule.domain().get_space());
      isl::local_space LSOut = isl::local_space(Schedule.range().get_space());
      isl::pw_aff VarIn = isl::pw_aff::var_on_domain(LSIn, isl::dim::set, Schedule.dim(isl::dim::in) - 1);
      isl::pw_aff VarOut = isl::pw_aff::var_on_domain(LSOut, isl::dim::set, Schedule.dim(isl::dim::out) - 1);
      isl::map EqualConstraint = VarIn.eq_map(VarOut);
      Schedule = Schedule.intersect(EqualConstraint);
      // compute the reference range
      isl::pw_aff Min = VarIn * 0;
      isl::pw_aff Max = VarIn * 0 + ReadReferences_[Statement].size() + WriteReferences_[Statement].size();
      isl::set LowerBound = VarIn.ge_set(Min);
      isl::set UpperBound = VarIn.lt_set(Max);
      Schedule = Schedule.intersect_domain(LowerBound);
      Schedule = Schedule.intersect_domain(UpperBound);
      ScheduleExt = ScheduleExt.unite(isl::union_map(Schedule));
    }
    return isl::stat::ok();
  };
  ScheduleMap_.foreach_map(extendSchedule);
  ScheduleMap_ = ScheduleExt;

  // extend the read map
  isl::union_map ReadsExt = isl::map::empty(Reads_.get_space());
  auto extendReads = [&](isl::map Reads) {
    std::string Statement = Reads.domain().unwrap().get_tuple_name(isl::dim::in);
    std::string Reference = Reads.domain().unwrap().get_tuple_name(isl::dim::out);
    std::string Array = Reads.get_tuple_name(isl::dim::out);
    if (ElementSizes_.count(Array) > 0) {
      // get the reference offset
      auto Iter = std::find(ReadReferences_[Statement].begin(), ReadReferences_[Statement].end(), Reference);
      assert(Iter != ReadReferences_[Statement].end());
      // extend the map with the offset
      Reads = Reads.insert_dims(isl::dim::in, Reads.dim(isl::dim::in), 1);
      Reads = Reads.set_tuple_name(isl::dim::in, Statement);
      // set the offset number
      isl::local_space LSIn = isl::local_space(Reads.domain().get_space());
      isl::pw_aff VarIn = isl::pw_aff::var_on_domain(LSIn, isl::dim::set, Reads.dim(isl::dim::in) - 1);
      isl::pw_aff Offset = VarIn * 0 + std::distance(ReadReferences_[Statement].begin(), Iter);
      isl::set EqualConstraint = VarIn.eq_set(Offset);
      Reads = Reads.intersect_domain(EqualConstraint);
      ReadsExt = ReadsExt.unite(isl::union_map(Reads));
    }
    return isl::stat::ok();
  };
  Reads_.foreach_map(extendReads);
  Reads_ = ReadsExt;

  // extend the write map
  isl::union_map WritesExt = isl::map::empty(Writes_.get_space());
  auto extendWrites = [&](isl::map Writes) {
    std::string Statement = Writes.domain().unwrap().get_tuple_name(isl::dim::in);
    std::string Reference = Writes.domain().unwrap().get_tuple_name(isl::dim::out);
    std::string Array = Writes.get_tuple_name(isl::dim::out);
    if (ElementSizes_.count(Array) > 0) {
      // get the reference offset
      auto Iter = std::find(WriteReferences_[Statement].begin(), WriteReferences_[Statement].end(), Reference);
      assert(Iter != WriteReferences_[Statement].end());
      // extend the map with the offset
      Writes = Writes.insert_dims(isl::dim::in, Writes.dim(isl::dim::in), 1);
      Writes = Writes.set_tuple_name(isl::dim::in, Statement);
      // set the offset number
      isl::local_space LSIn = isl::local_space(Writes.domain().get_space());
      isl::pw_aff VarIn = isl::pw_aff::var_on_domain(LSIn, isl::dim::set, Writes.dim(isl::dim::in) - 1);
      auto Distance = std::distance(WriteReferences_[Statement].begin(), Iter) + ReadReferences_[Statement].size();
      isl::pw_aff Offset = VarIn * 0 + Distance;
      isl::set EqualConstraint = VarIn.eq_set(Offset);
      Writes = Writes.intersect_domain(EqualConstraint);
      WritesExt = WritesExt.unite(isl::union_map(Writes));
    }
    return isl::stat::ok();
  };
  Writes_.foreach_map(extendWrites);
  Writes_ = WritesExt;

  // compute the access domain
  AccessDomain_ = Reads_.domain().unite(Writes_.domain()).coalesce();
}

void Program::computeAccessToLine(isl::set Parameters) {
  if (AccessToLine_.is_null()) {
    // compute the access map
    AccessToLine_ = isl::map::empty(Writes_.get_space());
    AccessToElement_ = isl::map::empty(Writes_.get_space());
    for (auto Array : ArrayExtents_) {
      // extract the array information
      std::string Name = Array.first;
      isl::set Extent = Array.second;
      // apply the parameters
      if (!Parameters.is_null()) {
        Extent = Extent.intersect_params(Parameters);
      }
      // map the access to the array offsets
      isl::map AccessToArray = isl::map::identity(Extent.get_space().map_from_set());
      AccessToElement_ = AccessToElement_.unite(isl::union_map(AccessToArray));
      // compute elements per cache line
      long ElementsPerCacheLine = MachineModel_.CacheLineSize / ElementSizes_[Name];
      AccessToArray = introduceCacheLines(Name, AccessToArray, ElementsPerCacheLine);
      AccessToLine_ = AccessToLine_.unite(isl::union_map(AccessToArray));
      AccessToLine_ = AccessToLine_.coalesce();
    }

    // compose the memory accesses with access map
    isl::union_map Accesses = Reads_.unite(Writes_);
    AccessToLine_ = Accesses.apply_range(AccessToLine_);
    AccessToLine_ = AccessToLine_.coalesce();
    AccessToElement_ = Accesses.apply_range(AccessToElement_);
    AccessToElement_ = AccessToElement_.coalesce();
  }
}

isl::map Program::introduceCacheLines(std::string Name, isl::map AccessToArray, long ElementsPerCacheLine) const {
  // introduce additional dimension that divides the innermost dimension by the cache line size
  int Dim = AccessToArray.dim(isl::dim::out) - 1;
  AccessToArray = AccessToArray.add_dims(isl::dim::out, 1);
  isl::local_space LS = isl::local_space(AccessToArray.range().get_space());
  isl::pw_aff Var = isl::pw_aff::var_on_domain(LS, isl::dim::set, Dim);
  isl::pw_aff VarPrime = isl::pw_aff::var_on_domain(LS, isl::dim::set, Dim + 1);
  isl::set ConstraintOne = (ElementsPerCacheLine * VarPrime).le_set(Var);
  isl::set ConstraintTwo = Var.le_set(ElementsPerCacheLine * VarPrime + (ElementsPerCacheLine - 1));
  AccessToArray = AccessToArray.intersect_range(ConstraintOne).intersect_range(ConstraintTwo);
  AccessToArray = AccessToArray.project_out(isl::dim::out, Dim, 1);
  AccessToArray = AccessToArray.set_tuple_name(isl::dim::out, Name);
  // return the resulting array
  return AccessToArray;
}

isl::map Program::introduceCacheSets(std::string Name, isl::map AccessToArray, long NumberOfCacheSets) const {
  // introduce additional dimension that divides the innermost dimension by the cache line size
  int Dim = AccessToArray.dim(isl::dim::out) - 1;
  AccessToArray = AccessToArray.add_dims(isl::dim::out, 1);
  isl::local_space LS = isl::local_space(AccessToArray.range().get_space());
  isl::pw_aff Var = isl::pw_aff::var_on_domain(LS, isl::dim::set, Dim);
  isl::pw_aff VarPrime = isl::pw_aff::var_on_domain(LS, isl::dim::set, Dim + 1);
  isl::pw_aff Modulo = Var.mod(isl::val(Var.get_ctx(), NumberOfCacheSets));
  isl::set Constraint = VarPrime.eq_set(Modulo);
  AccessToArray = AccessToArray.intersect_range(Constraint);
  AccessToArray = AccessToArray.set_tuple_name(isl::dim::out, Name);
  // return the resulting array
  return AccessToArray;
}
