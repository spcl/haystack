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
const int CacheSize = 256 / ElementSize;

// execute the toy code with the emulator
void emulateToy(int N, int CacheLineSize, CacheEmulator &Emulator) {
  int TimeStamp = 0;
  // define the index to cachline conversion
  auto CL = [&](int i) { return i / CacheLineSize; };
  // run toy
  for (int i = 0; i < N; i++) {
    // A[i] = 0;
    Emulator.accessMemory("S0", TimeStamp, {CL(i)});
    // A[N-1-i] = 1;
    Emulator.accessMemory("S1", TimeStamp, {CL(N - 1 - i)});
    if (i < N / 2)
      // A[2*i] = 2;
      Emulator.accessMemory("S2", TimeStamp, {CL(2 * i)});
  }
}

class ToyTest : public ::testing::Test {
protected:
  ToyTest() {
    Context_ = isl_ctx_alloc_with_pet_options();
    isl_options_set_on_error(Context_, ISL_ON_ERROR_ABORT);

    Base_ = new HayStack(Context_, {CacheLineSize * ElementSize, {CacheSize * ElementSize}}, {true});
    Base_->compileProgram("./toy.c");
  }

  virtual ~ToyTest() {
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

TEST_F(ToyTest, CapacityMissEven) {
  // define the setup
  int N = 128;
  int CacheLines = (N + CacheLineSize - 1) / CacheLineSize;
  std::vector<NamedLong> Parameters = {std::make_pair(std::string("N"), N)};
  // emulate the stack distances
  CacheEmulator Emulator(CacheLines, CacheSize / CacheLineSize);
  emulateToy(N, CacheLineSize, Emulator);
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

TEST_F(ToyTest, CompulsoryMissEven) {
  // define the setup
  int N = 128;
  int CacheLines = (N + CacheLineSize - 1) / CacheLineSize;
  std::vector<NamedLong> Parameters = {std::make_pair(std::string("N"), N)};
  // emulate the stack distances
  CacheEmulator Emulator(CacheLines, CacheSize / CacheLineSize);
  emulateToy(N, CacheLineSize, Emulator);
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

TEST_F(ToyTest, CapacityMissOdd) {
  // define the setup
  int N = 177;
  int CacheLines = (N + CacheLineSize - 1) / CacheLineSize;
  std::vector<NamedLong> Parameters = {std::make_pair(std::string("N"), N)};
  // emulate the stack distances
  CacheEmulator Emulator(CacheLines, CacheSize / CacheLineSize);
  emulateToy(N, CacheLineSize, Emulator);
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

TEST_F(ToyTest, CompulsoryMissOdd) {
  // define the setup
  int N = 177;
  int CacheLines = (N + CacheLineSize - 1) / CacheLineSize;
  std::vector<NamedLong> Parameters = {std::make_pair(std::string("N"), N)};
  // emulate the stack distances
  CacheEmulator Emulator(CacheLines, CacheSize / CacheLineSize);
  emulateToy(N, CacheLineSize, Emulator);
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
