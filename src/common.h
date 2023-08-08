#ifndef _COMMON_H
#define _COMMON_H

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
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

/* #define UINT32_MAX ((uint32_t)-1) */

#ifdef _OPENMP
  #include <omp.h>
#else
  #warning "OpenMP not found, multi-threading will be DISABLED and --num-threads option will be ignored!"
#endif

extern unsigned int num_threads;

template<typename T>
using vvec = std::vector<std::vector<T>>;

#endif
