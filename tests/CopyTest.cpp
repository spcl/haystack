/*
* Copyright (c) 2019, ETH Zurich
*/

#include "gtest/gtest.h"

#include <isl/options.h>

#include "../src/HayStack.h"
#include "CacheEmulator.h"

// test setup
const int ElementSize = 4;
const int CacheLineSize = 16 / ElementSize;
const int CacheSize = 1024 / ElementSize;

// execute the copy code with the emulator
void emulateCopy(int N1, int N2, int N3, int CacheLineSize, CacheEmulator &Emulator) {
  int TimeStamp = 0;
  int StrideA = (N3 + CacheLineSize - 1) / CacheLineSize;
  int StrideB = (N3 + CacheLineSize - 1) / CacheLineSize;
  int StrideT = (N1 + CacheLineSize - 1) / CacheLineSize;
  int OffsetB = N1 * N2 * StrideA;
  int OffsetT = OffsetB + N1 * N2 * StrideB;
  auto CLA = [&](int i, int j, int k) { return i * StrideA * N2 + j * StrideA + k / CacheLineSize; };
  auto CLB = [&](int i, int j, int k) { return OffsetB + i * StrideB * N2 + j * StrideB + k / CacheLineSize; };
  auto CLT = [&](int i, int j, int k) { return OffsetT + i * StrideT * N2 + j * StrideT + k / CacheLineSize; };
  // run copy
  for (int i = 0; i < N1; i++)
    for (int j = 0; j < N2; j++)
      for (int k = 0; k < N3; k++)
        // T[k][j][i] = A[i][j][k];
        Emulator.accessMemory("S0", TimeStamp, {CLA(i, j, k), CLT(k, j, i)});
  for (int i = 0; i < N1; i++)
    for (int j = 0; j < N2; j++)
      for (int k = 0; k < N3; k++)
        // B[i][j][k] = T[k][j][i];
        Emulator.accessMemory("S1", TimeStamp, {CLT(k, j, i), CLB(i, j, k)});
}

class CopyTest : public ::testing::Test {
protected:
  CopyTest() {
    Context_ = isl_ctx_alloc_with_pet_options();
    isl_options_set_on_error(Context_, ISL_ON_ERROR_ABORT);

    Base_ = new HayStack(Context_, {CacheLineSize * ElementSize, {CacheSize * ElementSize}}, {true});
    Base_->compileProgram("./copy.c");
  }

  virtual ~CopyTest() {
    delete Base_;
    isl_ctx_free(Context_);
  }

  virtual void SetUp() {
    // get fresh copy for every test
    Model_ = new HayStack(*Base_);
  }

  virtual void TearDown() { delete Model_; }

  isl_ctx *Context_;
  HayStack *Base_;
  HayStack *Model_;
};

TEST_F(CopyTest, CapacityMissEven) {
  // define the setup
  int N1 = 32;
  int N2 = 16;
  int N3 = 8;
  int CacheLines =
      2 * N1 * N2 * ((N3 + CacheLineSize - 1) / CacheLineSize) + N2 * N3 * ((N1 + CacheLineSize - 1) / CacheLineSize);
  std::vector<NamedLong> Parameters = {std::make_pair(std::string("N1"), N1), std::make_pair(std::string("N2"), N2),
                                       std::make_pair(std::string("N3"), N3)};
  // emulate the stack distances
  CacheEmulator Emulator(CacheLines, CacheSize / CacheLineSize);
  emulateCopy(N1, N2, N3, CacheLineSize, Emulator);
  auto ExpectedCapacityMisses = Emulator.getCapacityMisses();
  // compute the stack distances
  Model_->initModel(Parameters);
  //
  std::map<std::string, std::vector<long>> ComputedCapacityMisses;
  for (auto ComputedCapacityMiss : Model_->countCacheMisses()) {
    std::string Statement = ComputedCapacityMiss.first;
    auto length = Statement.find_first_of("(");
    Statement = Statement.substr(0, length);
    ComputedCapacityMisses[Statement].push_back(ComputedCapacityMiss.second.CapacityMisses[0]);
  }

  // print computed and expected stack distances
  for (auto ComputedCapacityMiss : ComputedCapacityMisses) {
    printf("Computed %s -> ", ComputedCapacityMiss.first.c_str());
    for (auto Distance : ComputedCapacityMiss.second)
      printf("%ld ", Distance);
    printf("\n");
  }
  for (auto ExpectedCapacityMiss : ExpectedCapacityMisses) {
    printf("Expected %s -> ", ExpectedCapacityMiss.first.c_str());
    for (auto Distance : ExpectedCapacityMiss.second)
      printf("%ld ", Distance);
    printf("\n");
  }

  // make sure the sizes agree
  ASSERT_EQ(ExpectedCapacityMisses.size(), ComputedCapacityMisses.size());

  // compare the stack distances for all statements
  for (auto ComputedCapacityMiss : ComputedCapacityMisses) {
    auto ExpectedCapacityMiss = ExpectedCapacityMisses[ComputedCapacityMiss.first];
    ASSERT_EQ(ExpectedCapacityMiss.size(), ComputedCapacityMiss.second.size());

    for (int i = 0; i < ComputedCapacityMiss.second.size(); ++i)
      EXPECT_EQ(ExpectedCapacityMiss[i], ComputedCapacityMiss.second[i]);
  }
}

TEST_F(CopyTest, CompulsoryMissEven) {
  // define the setup
  int N1 = 32;
  int N2 = 16;
  int N3 = 8;
  int CacheLines =
      2 * N1 * N2 * ((N3 + CacheLineSize - 1) / CacheLineSize) + N3 * N2 * ((N1 + CacheLineSize - 1) / CacheLineSize);
  std::vector<NamedLong> Parameters = {std::make_pair(std::string("N1"), N1), std::make_pair(std::string("N2"), N2),
                                       std::make_pair(std::string("N3"), N3)};
  // emulate the stack distances
  CacheEmulator Emulator(CacheLines, CacheSize / CacheLineSize);
  emulateCopy(N1, N2, N3, CacheLineSize, Emulator);
  auto ExpectedCompulsoryMisses = Emulator.getCompulsoryMisses();
  // compute the stack distances
  Model_->initModel(Parameters);
  //
  std::map<std::string, std::vector<long>> ComputedCompulsoryMisses;
  for (auto ComputedCompulsoryMiss : Model_->countCacheMisses()) {
    std::string Statement = ComputedCompulsoryMiss.first;
    auto length = Statement.find_first_of("(");
    Statement = Statement.substr(0, length);
    ComputedCompulsoryMisses[Statement].push_back(ComputedCompulsoryMiss.second.CompulsoryMisses);
  }

  // print computed and expected stack distances
  for (auto ComputedCompulsoryMiss : ComputedCompulsoryMisses) {
    printf("Computed %s -> ", ComputedCompulsoryMiss.first.c_str());
    for (auto Distance : ComputedCompulsoryMiss.second)
      printf("%ld ", Distance);
    printf("\n");
  }
  for (auto ExpectedCompulsoryMiss : ExpectedCompulsoryMisses) {
    printf("Expected %s -> ", ExpectedCompulsoryMiss.first.c_str());
    for (auto Distance : ExpectedCompulsoryMiss.second)
      printf("%ld ", Distance);
    printf("\n");
  }

  // make sure the sizes agree
  ASSERT_EQ(ExpectedCompulsoryMisses.size(), ComputedCompulsoryMisses.size());

  // compare the stack distances for all statements
  for (auto ComputedCompulsoryMiss : ComputedCompulsoryMisses) {
    auto ExpectedCompulsoryMiss = ExpectedCompulsoryMisses[ComputedCompulsoryMiss.first];
    ASSERT_EQ(ExpectedCompulsoryMiss.size(), ComputedCompulsoryMiss.second.size());

    for (int i = 0; i < ComputedCompulsoryMiss.second.size(); ++i)
      EXPECT_EQ(ExpectedCompulsoryMiss[i], ComputedCompulsoryMiss.second[i]);
  }
}

TEST_F(CopyTest, CapacityMissOdd) {
  // define the setup
  int N1 = 11;
  int N2 = 22;
  int N3 = 33;
  int CacheLines =
      2 * N1 * N2 * ((N3 + CacheLineSize - 1) / CacheLineSize) + N3 * N2 * ((N1 + CacheLineSize - 1) / CacheLineSize);
  std::vector<NamedLong> Parameters = {std::make_pair(std::string("N1"), N1), std::make_pair(std::string("N2"), N2),
                                       std::make_pair(std::string("N3"), N3)};
  // emulate the stack distances
  CacheEmulator Emulator(CacheLines, CacheSize / CacheLineSize);
  emulateCopy(N1, N2, N3, CacheLineSize, Emulator);
  auto ExpectedCapacityMisses = Emulator.getCapacityMisses();
  // compute the stack distances
  Model_->initModel(Parameters);
  //
  std::map<std::string, std::vector<long>> ComputedCapacityMisses;
  for (auto ComputedCapacityMiss : Model_->countCacheMisses()) {
    std::string Statement = ComputedCapacityMiss.first;
    auto length = Statement.find_first_of("(");
    Statement = Statement.substr(0, length);
    ComputedCapacityMisses[Statement].push_back(ComputedCapacityMiss.second.CapacityMisses[0]);
  }

  // print computed and expected stack distances
  for (auto ComputedCapacityMiss : ComputedCapacityMisses) {
    printf("Computed %s -> ", ComputedCapacityMiss.first.c_str());
    for (auto Distance : ComputedCapacityMiss.second)
      printf("%ld ", Distance);
    printf("\n");
  }
  for (auto ExpectedCapacityMiss : ExpectedCapacityMisses) {
    printf("Expected %s -> ", ExpectedCapacityMiss.first.c_str());
    for (auto Distance : ExpectedCapacityMiss.second)
      printf("%ld ", Distance);
    printf("\n");
  }

  // make sure the sizes agree
  ASSERT_EQ(ExpectedCapacityMisses.size(), ComputedCapacityMisses.size());

  // compare the stack distances for all statements
  for (auto ComputedCapacityMiss : ComputedCapacityMisses) {
    auto ExpectedCapacityMiss = ExpectedCapacityMisses[ComputedCapacityMiss.first];
    ASSERT_EQ(ExpectedCapacityMiss.size(), ComputedCapacityMiss.second.size());

    for (int i = 0; i < ComputedCapacityMiss.second.size(); ++i)
      EXPECT_EQ(ExpectedCapacityMiss[i], ComputedCapacityMiss.second[i]);
  }
}

TEST_F(CopyTest, CompulsoryMissOdd) {
  // define the setup
  int N1 = 11;
  int N2 = 22;
  int N3 = 33;
  int CacheLines =
      2 * N1 * N2 * ((N3 + CacheLineSize - 1) / CacheLineSize) + N3 * N2 * ((N1 + CacheLineSize - 1) / CacheLineSize);
  std::vector<NamedLong> Parameters = {std::make_pair(std::string("N1"), N1), std::make_pair(std::string("N2"), N2),
                                       std::make_pair(std::string("N3"), N3)};
  // emulate the stack distances
  CacheEmulator Emulator(CacheLines, CacheSize / CacheLineSize);
  emulateCopy(N1, N2, N3, CacheLineSize, Emulator);
  auto ExpectedCompulsoryMisses = Emulator.getCompulsoryMisses();
  // compute the stack distances
  Model_->initModel(Parameters);
  //
  std::map<std::string, std::vector<long>> ComputedCompulsoryMisses;
  for (auto ComputedCompulsoryMiss : Model_->countCacheMisses()) {
    std::string Statement = ComputedCompulsoryMiss.first;
    auto length = Statement.find_first_of("(");
    Statement = Statement.substr(0, length);
    ComputedCompulsoryMisses[Statement].push_back(ComputedCompulsoryMiss.second.CompulsoryMisses);
  }

  // print computed and expected stack distances
  for (auto ComputedCompulsoryMiss : ComputedCompulsoryMisses) {
    printf("Computed %s -> ", ComputedCompulsoryMiss.first.c_str());
    for (auto Distance : ComputedCompulsoryMiss.second)
      printf("%ld ", Distance);
    printf("\n");
  }
  for (auto ExpectedCompulsoryMiss : ExpectedCompulsoryMisses) {
    printf("Expected %s -> ", ExpectedCompulsoryMiss.first.c_str());
    for (auto Distance : ExpectedCompulsoryMiss.second)
      printf("%ld ", Distance);
    printf("\n");
  }

  // make sure the sizes agree
  ASSERT_EQ(ExpectedCompulsoryMisses.size(), ComputedCompulsoryMisses.size());

  // compare the stack distances for all statements
  for (auto ComputedCompulsoryMiss : ComputedCompulsoryMisses) {
    auto ExpectedCompulsoryMiss = ExpectedCompulsoryMisses[ComputedCompulsoryMiss.first];
    ASSERT_EQ(ExpectedCompulsoryMiss.size(), ComputedCompulsoryMiss.second.size());

    for (int i = 0; i < ComputedCompulsoryMiss.second.size(); ++i)
      EXPECT_EQ(ExpectedCompulsoryMiss[i], ComputedCompulsoryMiss.second[i]);
  }
}
