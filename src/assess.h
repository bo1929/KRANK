#ifndef _ASSESS_H
#define _ASSESS_H
#include "common.h"

inline unsigned int
computeFastModulo(const unsigned int input, const unsigned int ceil)
{
  return input >= ceil ? input % ceil : input;
}

void
getIxsRandom(std::vector<unsigned int>& ixs, uint8_t b_curr, uint8_t b_max);

template<typename T>
void
vecRemoveIxs(std::vector<T>& v, std::vector<unsigned int>& ixs);

template<typename T>
void
arrRemoveIxs(T* arr, uint8_t last, std::vector<unsigned int>& ixs);

template<typename T>
void
vecArgsort1D(std::vector<unsigned int>& ixs, const std::vector<T>& v, bool reverse = false);

template<typename T>
void
arrArgsort1D(std::vector<unsigned int>& ixs, const T* arr, uint8_t last, bool reverse = false);

template<typename T>
T
vvecArgmax2D(const vvec<T>& vv, uint64_t n, bool reverse = false);

template<typename T>
T
arrArgmax2D(const T* arr, const uint8_t* ind_arr, uint32_t num_rows, uint8_t b, uint64_t n, bool reverse = false);

template<typename T>
void
vecIxsNumber(std::vector<unsigned int>& ixs, const std::vector<T>& s_v, uint8_t number);

template<typename T>
void
vecIxsThreshold(std::vector<unsigned int>& ixs, const std::vector<T>& s_v, T threshold);

template<typename T>
void
arrIxsNumber(std::vector<unsigned int>& ixs, const T* s_arr, uint8_t number);

template<typename T>
void
arrIxsThreshold(std::vector<unsigned int>& ixs, const T* s_arr, uint8_t last, T threshold);

template<typename T>
void
vvecSizeOrder(std::vector<unsigned int>& ixs, const vvec<T>& vv, bool reverse = true);

#endif
