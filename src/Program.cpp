/*
* Copyright (c) 2019, ETH Zurich
*/

#include <algorithm>

#include "Program.h"
#include "isl-helpers.h"
#include "op.h"

void Program::extractScop(std::string SourceFile, std::string ScopFunction) {
  const char *Function = ScopFunction.empty() ? NULL : ScopFunction.c_str();
  pet_scop *PetScop = pet_scop_extract_from_C_source(Context_.get(), SourceFile.c_str(), Function);
  if (PetScop == nullptr) {
    printf("-> exit(-1) cannot extract scope\n");
    exit(-1);
  }

  // get the schedule and access information
  // (the tagged accesses allow us to distinguish multiple accesses of the same array)
  Schedule_ = isl::manage(pet_scop_get_schedule(PetScop));
  Reads_ = isl::manage(pet_scop_get_tagged_may_reads(PetScop));
  Writes_ = isl::manage(pet_scop_get_tagged_may_writes(PetScop));

  // check if the schedule is bounded
  auto checkIfBounded = [](isl::set Set) {
    if (!Set.is_bounded()) {
      printf("-> exit(-1) schedule not bounded\n");
      exit(-1);
    }
    return isl::stat::ok();
  };
  Schedule_.get_domain().foreach_set(checkIfBounded);

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
  // extract the reads and writes
  extractReferences();

  // compute detailed access information
  ScopLoc_ = std::make_pair(pet_loc_get_start(PetScop->loc), pet_loc_get_end(PetScop->loc));
  for (int idx = 0; idx < PetScop->n_stmt; idx++) {
    // extract the statement info
    pet_expr *Expression = pet_tree_expr_get_expr(PetScop->stmts[idx]->body);
    isl::space Space = isl::manage(pet_stmt_get_space(PetScop->stmts[idx]));
    std::string Statement = Space.get_tuple_name(isl::dim::set);
    // extract the access info
    auto printExpression = [](__isl_keep pet_expr *Expr, void *User) {
      if (pet_expr_access_is_read(Expr) || pet_expr_access_is_write(Expr)) {
        isl::id RefId = isl::manage(pet_expr_access_get_ref_id(Expr));
        std::string Name = RefId.to_str();
        isl::multi_pw_aff Index = isl::manage(pet_expr_access_get_index(Expr));
        // filter the array accesses
        if (Index.dim(isl::dim::out) > 0 && Index.has_tuple_id(isl::dim::out)) {
          std::string Access = Index.get_tuple_name(isl::dim::out);
          // process the array dimensions
          for (int i = 0; i < Index.dim(isl::dim::out); ++i) {
            std::vector<std::string> Expressions;
            auto IndexExpr = Index.get_pw_aff(i);
            auto extractExpr = [&](isl::set Set, isl::aff Aff) {
              Access += "[" + isl::printExpression(Aff) + "]";
              return isl::stat::ok();
            };
            IndexExpr.foreach_piece(extractExpr);
          }
          // find the access info
          auto AccessInfos = (std::vector<access_info> *)User;
          auto Iter = AccessInfos->begin();
          do {
            // find the pet references and update the access description
            Iter = std::find_if(Iter, AccessInfos->end(), [&](access_info AccessInfo) {
              if (AccessInfo.Access == Name)
                return true;
              return false;
            });
            if (Iter != AccessInfos->end() && Iter->Access == Name) {
              Iter->Access = Access;
              Iter++;
            }
          } while (Iter != AccessInfos->end());
        }
      }
      return 0;
    };
    pet_expr_foreach_access_expr(Expression, printExpression, &AccessInfos_[Statement]);
    // get the line number
    pet_loc *Loc = pet_tree_get_loc(PetScop->stmts[idx]->body);
    int Line = pet_loc_get_line(Loc);
    int Start = pet_loc_get_start(Loc);
    int Stop = pet_loc_get_end(Loc);
    pet_loc_free(Loc);
    // store the access information
    for (auto &Access : AccessInfos_[Statement]) {
      Access.Line = Line;
      Access.Start = Start;
      Access.Stop = Stop;
    }
    pet_expr_free(Expression);
  }

  // extend the schedule and the read and write maps with an access dimension
  extendSchedule();
  Reads_ = extendAccesses(Reads_, false);
  Writes_ = extendAccesses(Writes_, true);

  // compute the access domain
  AccessDomain_ = Reads_.domain().unite(Writes_.domain()).coalesce();

  // free the pet scop
  pet_scop_free(PetScop);
}

void Program::extendSchedule() {
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
}

void Program::extractReferences() {
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
  auto compare = [](std::string &S1, std::string &S2) {
    if (S1.length() == S2.length())
      return S1 < S2;
    else
      return S1.length() < S2.length();
  };
  for (auto &ReadReferences : ReadReferences_)
    std::sort(ReadReferences.second.begin(), ReadReferences.second.end(), compare);
  for (auto &WriteReferences : WriteReferences_)
    std::sort(WriteReferences.second.begin(), WriteReferences.second.end(), compare);

  // compute the access info
  AccessInfos_.clear();
  for (auto &ReadReferences : ReadReferences_) {
    for (int i = 0; i < ReadReferences.second.size(); ++i) {
      AccessInfos_[ReadReferences.first].push_back(
          {ReadReferences.first + "(R" + std::to_string(i) + ")", ReadReferences.second[i], Read, 0, 0, 0});
    }
  }
  for (auto &WriteReferences : WriteReferences_) {
    for (int i = 0; i < WriteReferences.second.size(); ++i) {
      AccessInfos_[WriteReferences.first].push_back(
          {WriteReferences.first + "(W" + std::to_string(i) + ")", WriteReferences.second[i], Write, 0, 0, 0});
    }
  }
}

isl::union_map Program::extendAccesses(isl::union_map Accesses, bool WriteReferences) {
  // extend the access map
  isl::union_map AccessesExt = isl::map::empty(Accesses.get_space());
  auto extendAccesses = [&](isl::map Access) {
    std::string Statement = Access.domain().unwrap().get_tuple_name(isl::dim::in);
    std::string Reference = Access.domain().unwrap().get_tuple_name(isl::dim::out);
    std::string Array = Access.get_tuple_name(isl::dim::out);
    if (ElementSizes_.count(Array) > 0) {
      // get the reference offset
      long Distance = 0;
      if (WriteReferences) {
        auto Iter = std::find(WriteReferences_[Statement].begin(), WriteReferences_[Statement].end(), Reference);
        assert(Iter != WriteReferences_[Statement].end());
        Distance = std::distance(WriteReferences_[Statement].begin(), Iter) + ReadReferences_[Statement].size();
      } else {
        auto Iter = std::find(ReadReferences_[Statement].begin(), ReadReferences_[Statement].end(), Reference);
        assert(Iter != ReadReferences_[Statement].end());
        Distance = std::distance(ReadReferences_[Statement].begin(), Iter);
      }
      // extend the map with the offset
      Access = Access.insert_dims(isl::dim::in, Access.dim(isl::dim::in), 1);
      Access = Access.set_tuple_name(isl::dim::in, Statement);
      // set the offset number
      isl::local_space LSIn = isl::local_space(Access.domain().get_space());
      isl::pw_aff VarIn = isl::pw_aff::var_on_domain(LSIn, isl::dim::set, Access.dim(isl::dim::in) - 1);
      isl::pw_aff Offset = VarIn * 0 + Distance;
      isl::set EqualConstraint = VarIn.eq_set(Offset);
      Access = Access.intersect_domain(EqualConstraint);
      AccessesExt = AccessesExt.unite(isl::union_map(Access));
    }
    return isl::stat::ok();
  };
  Accesses.foreach_map(extendAccesses);
  return AccessesExt;
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
