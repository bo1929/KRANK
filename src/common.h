#ifndef _COMMON_H
#define _COMMON_H

#include <algorithm>
#include <assert.h>
#include <bitset>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <ostream>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <type_traits>
#include <dirent.h>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <cstdlib>
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
#include <zlib.h>
#include "aixlog.hpp"

#ifdef _OPENMP
  #include <omp.h>
#else
  #warning "OpenMP not found, multi-threading will be DISABLED and --num-threads option will be ignored!"
#endif

#define VERSION "v0.0.1"
#define PRINT_VERSION printf("KRANK version: " VERSION "\n");

#define MAX_NUM_TASKS 8

extern unsigned int num_threads;
extern std::random_device rd;
extern std::mt19937 gen;
extern std::bernoulli_distribution ct;
extern const unsigned char seq_nt4_table[128];

template<typename T>
using vvec = std::vector<std::vector<T>>;

#ifdef LARGE_TAXONOMY
typedef uint32_t tT;
typedef uint32_t scT;
#else
typedef uint16_t tT;
typedef uint16_t scT;
#endif
#ifdef SHORT_TABLE
typedef uint64_t encT;
#else
typedef uint32_t encT;
#endif

namespace std {
  template<class T>
  struct _Unique_if
  {
    typedef unique_ptr<T> _Single_object;
  };

  template<class T>
  struct _Unique_if<T[]>
  {
    typedef unique_ptr<T[]> _Unknown_bound;
  };

  template<class T, size_t N>
  struct _Unique_if<T[N]>
  {
    typedef void _Known_bound;
  };

  template<class T, class... Args>
  typename _Unique_if<T>::_Single_object make_unique(Args &&...args)
  {
    return unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

  template<class T>
  typename _Unique_if<T>::_Unknown_bound make_unique(size_t n)
  {
    typedef typename remove_extent<T>::type U;
    return unique_ptr<T>(new U[n]());
  }

  template<class T, class... Args>
  typename _Unique_if<T>::_Known_bound make_unique(Args &&...) = delete;
}

template<typename Map>
bool map_compare(Map const &lhs, Map const &rhs)
{
  return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

#endif
