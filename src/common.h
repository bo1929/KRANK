#ifndef _COMMON_H
#define _COMMON_H

#include "aixlog.hpp"
#include <algorithm>
#include <assert.h>
#include <bitset>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <dirent.h>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <new>
#include <ostream>
#include <random>
#include <set>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <type_traits>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

/* #define UINT32_MAX ((uint32_t)-1) */

#ifdef _OPENMP
  #include <omp.h>
#else
  #warning                                                                                         \
    "OpenMP not found, multi-threading will be DISABLED and --num-threads option will be ignored!"
#endif

extern unsigned int num_threads;
extern std::random_device rd;
extern std::mt19937 gen;

template<typename T>
using vvec = std::vector<std::vector<T>>;

// TODO: Check if types are adequate.
// TODO: Is it possible to set types automatically?
typedef uint16_t tT;
typedef uint64_t encT;
typedef uint16_t scT;

class AddTimeStamp : public std::streambuf
{
public:
  AddTimeStamp(std::basic_ios<char>& out)
    : out_(out)
    , sink_()
    , newline_(false)
  {
    sink_ = out_.rdbuf(this);
    assert(sink_);
  }
  ~AddTimeStamp() { out_.rdbuf(sink_); }

protected:
  int_type overflow(int_type m = traits_type::eof())
  {
    if (traits_type::eq_int_type(m, traits_type::eof()))
      return sink_->pubsync() == -1 ? m : traits_type::not_eof(m);
    if (newline_) {
      std::ostream str(sink_);
      auto curr_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
      if (!(str << std::ctime(&curr_time)))
        return traits_type::eof();
    }
    newline_ = traits_type::to_char_type(m) == '\n';
    return sink_->sputc(m);
  }

private:
  AddTimeStamp(const AddTimeStamp&);
  AddTimeStamp& operator=(const AddTimeStamp&);
  std::basic_ios<char>& out_;
  std::streambuf* sink_;
  bool newline_;
};

#endif
