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
  // concatenate reads and writes
  AllReferences_ = ReadReferences_;
  for (auto &WriteReference : WriteReferences_)
    std::copy(WriteReference.second.begin(), WriteReference.second.end(),
              std::back_inserter(AllReferences_[WriteReference.first]));

  // extract the access information
  ScopLoc_ = std::make_pair(pet_loc_get_start(PetScop->loc), pet_loc_get_end(PetScop->loc));
  for (int idx = 0; idx < PetScop->n_stmt; idx++) {
    pet_expr *Expr = pet_tree_expr_get_expr(PetScop->stmts[idx]->body);
    std::map<std::string, std::string> Accesses;
    auto printExpr = [](__isl_keep pet_expr *Expr, void *User) {
      if (pet_expr_access_is_read(Expr) || pet_expr_access_is_write(Expr)) {
        isl::id RefId = isl::manage(pet_expr_access_get_ref_id(Expr));
        isl::multi_pw_aff Index = isl::manage(pet_expr_access_get_index(Expr));
        // filter the array accesses
        if (Index.dim(isl::dim::out) > 0 && Index.has_tuple_id(isl::dim::out)) {
          std::string Name = RefId.to_str();
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
          static_cast<std::map<std::string, std::string> *>(User)->operator[](Name) = Access;
        }
      }
      return 0;
    };
    pet_expr_foreach_access_expr(Expr, printExpr, &Accesses);
    // get the line number
    pet_loc *Loc = pet_tree_get_loc(PetScop->stmts[idx]->body);
    int Line = pet_loc_get_line(Loc);
    pet_loc_free(Loc);
    // store the access information
    for (auto Access : Accesses) {
      AccessInfo_[Access.first] = {Access.second, Line};
    }
    pet_expr_free(Expr);
  }

  // TODO
  // map access name to
  // - start and stop indexes instead of line numbers
  // - use line numer in yaml!
  // - array access expression
  // - order for every line the statements by the access name

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
  // extend the write map
  Reads_ = extendAccesses(Reads_);
  Writes_ = extendAccesses(Writes_);
  // compute the access domain
  AccessDomain_ = Reads_.domain().unite(Writes_.domain()).coalesce();

  // free the pet scop
  pet_scop_free(PetScop);
}

isl::union_map Program::extendAccesses(isl::union_map Accesses) {
  // extend the access map
  isl::union_map AccessesExt = isl::map::empty(Accesses.get_space());
  auto extendAccesses = [&](isl::map Access) {
    std::string Statement = Access.domain().unwrap().get_tuple_name(isl::dim::in);
    std::string Reference = Access.domain().unwrap().get_tuple_name(isl::dim::out);
    std::string Array = Access.get_tuple_name(isl::dim::out);
    if (ElementSizes_.count(Array) > 0) {
      // get the reference offset
      auto Iter = std::find(AllReferences_[Statement].begin(), AllReferences_[Statement].end(), Reference);
      assert(Iter != AllReferences_[Statement].end());
      // extend the map with the offset
      Access = Access.insert_dims(isl::dim::in, Access.dim(isl::dim::in), 1);
      Access = Access.set_tuple_name(isl::dim::in, Statement);
      // set the offset number
      isl::local_space LSIn = isl::local_space(Access.domain().get_space());
      isl::pw_aff VarIn = isl::pw_aff::var_on_domain(LSIn, isl::dim::set, Access.dim(isl::dim::in) - 1);
      auto Distance = std::distance(AllReferences_[Statement].begin(), Iter);
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
