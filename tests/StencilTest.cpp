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

// execute stencil with the emulator
void emulateStencil(int N1, int N2, int CacheLineSize, CacheEmulator &Emulator) {
  int TimeStamp = 0;
  // define the index to cachline conversion
  int Stride = (N2 + CacheLineSize - 1) / CacheLineSize;
  int Offset = N2 * Stride;
  auto CLA = [&](int i, int j) { return i * Stride + j / CacheLineSize; };
  auto CLB = [&](int i, int j) { return Offset + i * Stride + j / CacheLineSize; };
  // run stencil
  for (int t = 0; t < N1; t++) {
    for (int i = 1; i < N2 - 1; i++)
      for (int j = 1; j < N2 - 1; j++)
        // B[i][j] = 0.2 * (A[i][j] + A[i][j - 1] + A[i][1 + j] + A[1 + i][j] + A[i - 1][j]);
        Emulator.accessMemory("S0", TimeStamp,
                              {CLA(i, j), CLA(i, j - 1), CLA(i, j + 1), CLA(i + 1, j), CLA(i - 1, j), CLB(i, j)});
    for (int i = 1; i < N2 - 1; i++)
      for (int j = 1; j < N2 - 1; j++)
        // A[i][j] = 0.2 * (B[i][j] + B[i][j - 1] + B[i][1 + j] + B[1 + i][j] + B[i - 1][j]);
        Emulator.accessMemory("S1", TimeStamp,
                              {CLB(i, j), CLB(i, j - 1), CLB(i, j + 1), CLB(i + 1, j), CLB(i - 1, j), CLA(i, j)});
  }
}

class StencilTest : public ::testing::Test {
protected:
  StencilTest() {
    Context_ = isl_ctx_alloc_with_pet_options();
    isl_options_set_on_error(Context_, ISL_ON_ERROR_ABORT);

    Base_ = new HayStack(Context_, {CacheLineSize * ElementSize, {CacheSize * ElementSize}}, {true});
    Base_->compileProgram("./stencil.c");
  }

  virtual ~StencilTest() {
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

TEST_F(StencilTest, CapacityMissEven) {
  // define the setup
  int N1 = 2;
  int N2 = 32;
  int CacheLines = 2 * N2 * ((N2 + CacheLineSize - 1) / CacheLineSize);
  std::vector<NamedLong> Parameters = {std::make_pair(std::string("N1"), N1), std::make_pair(std::string("N2"), N2)};
  // emulate the stack distances
  CacheEmulator Emulator(CacheLines, CacheSize / CacheLineSize);
  emulateStencil(N1, N2, CacheLineSize, Emulator);
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

TEST_F(StencilTest, CompulsoryMissEven) {
  // define the setup
  int N1 = 2;
  int N2 = 32;
  int CacheLines = 2 * N2 * ((N2 + CacheLineSize - 1) / CacheLineSize);
  std::vector<NamedLong> Parameters = {std::make_pair(std::string("N1"), N1), std::make_pair(std::string("N2"), N2)};
  // emulate the stack distances
  CacheEmulator Emulator(CacheLines, CacheSize / CacheLineSize);
  emulateStencil(N1, N2, CacheLineSize, Emulator);
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

TEST_F(StencilTest, CapacityMissOdd) {
  // define the setup
  int N1 = 2;
  int N2 = 33;
  int CacheLines = 2 * N2 * ((N2 + CacheLineSize - 1) / CacheLineSize);
  std::vector<NamedLong> Parameters = {std::make_pair(std::string("N1"), N1), std::make_pair(std::string("N2"), N2)};
  // emulate the stack distances
  CacheEmulator Emulator(CacheLines, CacheSize / CacheLineSize);
  emulateStencil(N1, N2, CacheLineSize, Emulator);
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

TEST_F(StencilTest, CompulsoryMissOdd) {
  // define the setup
  int N1 = 2;
  int N2 = 33;
  int CacheLines = 2 * N2 * ((N2 + CacheLineSize - 1) / CacheLineSize);
  std::vector<NamedLong> Parameters = {std::make_pair(std::string("N1"), N1), std::make_pair(std::string("N2"), N2)};
  // emulate the stack distances
  CacheEmulator Emulator(CacheLines, CacheSize / CacheLineSize);
  emulateStencil(N1, N2, CacheLineSize, Emulator);
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
