/*
* Copyright (c) 2019, ETH Zurich
*/

#ifndef _TIMER_H_
#define _TIMER_H_

#include <chrono>
#include <iostream>
#include <map>
#include <string>

// class that allows to profile the execution times
class Timer {
public:
  typedef std::pair<std::chrono::time_point<std::chrono::high_resolution_clock>, double> Entry;

  static inline void startTimer(std::string Name) {
#ifdef TIMERS
    Clocks[Name].first = std::chrono::high_resolution_clock::now();
#endif
  }

  static inline void stopTimer(std::string Name) {
#ifdef TIMERS
    auto stop = std::chrono::high_resolution_clock::now();
    Clocks[Name].second += std::chrono::duration<double, std::milli>(stop - Clocks[Name].first).count();
#endif
  }

  static inline void printClocks() {
#ifdef TIMERS
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(3);
    std::cout << "==================================================" << std::endl;
    for (auto Clock : Clocks) {
      std::cout << " - " << Clock.first << ":\t" << Clock.second.second << "ms" << std::endl;
    }
    std::cout << "==================================================" << std::endl;
#endif
  }

private:
  static std::map<std::string, Entry> Clocks;
};

#endif