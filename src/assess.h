#ifndef _ASSESS_H
#define _ASSESS_H

#include "common.h"

inline unsigned int computeFastModulo(const unsigned int input, const unsigned int ceil)
{
  return input >= ceil ? input % ceil : input;
}

void getIxsRandom(std::vector<size_t> &ixs, size_t b, size_t n);

template<typename T>
void getIxsArgsort(std::vector<size_t> &ixs, const std::vector<T> &v, size_t n, bool reverse = true);

template<typename T>
void vecRemoveIxs(std::vector<T> &v, std::vector<size_t> &ixs);

template<typename T>
void arrRemoveIxs(T *arr, uint8_t last, std::vector<size_t> &ixs);

template<typename T>
void vecArgsort(std::vector<size_t> &ixs, const std::vector<T> &v, bool reverse = false);

template<typename T>
void arrArgsort(std::vector<size_t> &ixs, const T *arr, uint8_t last, bool reverse = false);

template<typename T>
T histArgmax(const std::map<T, uint64_t> val_counts, uint64_t n, bool reverse = false);

template<typename T>
T vvecArgmax(const vvec<T> &vv, uint64_t n, bool reverse = false);

#endif
