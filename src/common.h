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
#include <zlib.h>

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
extern std::bernoulli_distribution ct;
extern const unsigned char seq_nt4_table[128];

template<typename T>
using vvec = std::vector<std::vector<T>>;

// TODO: Check if types are adequate.
// TODO: Is it possible to set types automatically?
typedef uint16_t tT;
typedef uint32_t encT;
typedef uint16_t scT;

template<class T, class... Args>
typename _Unique_if<T>::_Known_bound
make_unique(Args&&...) = delete;
} // namespace std

template<typename Map>
bool
map_compare(Map const& lhs, Map const& rhs)
{
  return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

#endif
