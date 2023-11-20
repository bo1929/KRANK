#include "assess.h"

template<typename T>
void vecRemoveIxs(std::vector<T> &v, std::vector<size_t> &ixs)
{
  if (!v.empty() && (ixs.size() <= v.size())) {
    if (ixs.size() > 0) {
      std::vector<int> its(std::begin(ixs), std::end(ixs));
      its.push_back(-1);
      its.push_back(v.size());
      std::sort(its.begin(), its.end());
      typename std::vector<T>::iterator last = v.begin();
      for (size_t i = 1; i != its.size(); ++i) {
        size_t range_begin = its[i - 1] + 1;
        size_t range_end = its[i];
        if ((range_begin >= 0) && (range_end <= v.size())) {
          std::copy(v.begin() + range_begin, v.begin() + range_end, last);
          last += range_end - range_begin;
        } else {
          std::puts("The attempt of removing non-existent index has been ignored.\n");
        }
      }
      v.erase(last, v.end());
    }
  } else {
    std::puts("Given vector is empty or too small, no index has been removed.\n");
  }
}

template<typename T>
void arrRemoveIxs(T *r, uint8_t last, std::vector<size_t> &ixs)
{
  if ((last > 0) && (last >= ixs.size())) {
    if (!ixs.empty()) {
      std::vector<T> v(last - ixs.size());
      std::sort(ixs.begin(), ixs.end());
      unsigned int ix = 0;
      unsigned int iv = 0;
      for (unsigned int i = 0; i < last; ++i) {
        if ((ix < ixs.size()) && (i == ixs[ix])) {
          ix++;
        } else {
          v[iv] = r[i];
          iv++;
        }
      }
      std::copy(v.begin(), v.end(), r);
      if (ixs[ix - 1] != ixs.back()) {
        std::puts("The attempt of removing non-existent index has been ignored.\n");
      }
    }
  } else {
    std::puts("Given vector is empty or too small, no index has been removed.\n");
  }
}

template<typename T>
void vecArgsort(std::vector<size_t> &ixs, const std::vector<T> &v, bool reverse)
{
  ixs.resize(v.size());
  std::iota(ixs.begin(), ixs.end(), 0);
  std::random_shuffle(ixs.begin(), ixs.end());
  if (reverse)
    std::sort(ixs.begin(), ixs.end(), [&v](unsigned int i1, unsigned int i2) { return v[i1] > v[i2]; });
  else
    std::sort(ixs.begin(), ixs.end(), [&v](unsigned int i1, unsigned int i2) { return v[i1] < v[i2]; });
}

template<typename T>
void arrArgsort(std::vector<size_t> &ixs, const T *r, uint8_t last, bool reverse)
{
  ixs.resize(last);
  std::iota(ixs.begin(), ixs.end(), 0);
  std::random_shuffle(ixs.begin(), ixs.end());
  if (reverse)
    std::sort(ixs.begin(), ixs.end(), [&r](unsigned int i1, unsigned int i2) { return r[i1] > r[i2]; });
  else
    std::sort(ixs.begin(), ixs.end(), [&r](unsigned int i1, unsigned int i2) { return r[i1] < r[i2]; });
}

template<typename T>
T vvecArgmax(const vvec<T> &vv, uint64_t n, bool reverse)
{
  std::map<T, uint64_t> val_counts;
  for (auto &v : vv) {
    for (auto val : v)
      val_counts[val]++;
  }
  return histArgmax(val_counts, n, reverse);
}

template<typename T>
T histArgmax(const std::map<T, uint64_t> val_counts, uint64_t n, bool reverse)
{
  uint64_t cumsum = 0;
  T val_prev;
  if (reverse) {
    val_prev = 0;
    for (auto it = val_counts.cbegin(); (it != val_counts.cend()) && (cumsum < n); ++it) {
      cumsum += it->second;
      val_prev = it->first;
    }
  } else {
    val_prev = std::numeric_limits<T>::max();
    for (auto it = val_counts.crbegin(); (it != val_counts.crend()) && (cumsum < n); ++it) {
      cumsum += it->second;
      val_prev = it->first;
    }
  }
  assert(val_prev > 0);
  assert(val_prev < std::numeric_limits<T>::max());
  return val_prev;
}

void getIxsRandom(std::vector<size_t> &ixs, size_t b, size_t n)
{
  assert(b > 0);
  assert(n <= b);
  if (ixs.size() != b)
    ixs.resize(b);
  std::iota(ixs.begin(), ixs.end(), 0);
  std::random_shuffle(ixs.begin(), ixs.end());
  ixs.resize(n);
  std::sort(ixs.begin(), ixs.end());
}

template<typename T>
void getIxsArgsort(std::vector<size_t> &ixs, const std::vector<T> &v, size_t n, bool reverse)
{
  // Sort indices:
  ixs.resize(v.size());
  std::iota(ixs.begin(), ixs.end(), 0);
  std::random_shuffle(ixs.begin(), ixs.end());
  if (reverse)
    std::sort(ixs.begin(), ixs.end(), [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });
  else
    std::sort(ixs.begin(), ixs.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
  // Get first n:
  ixs.resize(n);
  std::sort(ixs.begin(), ixs.end());
}

template void getIxsArgsort(std::vector<size_t> &ixs, const std::vector<float> &v, size_t n, bool reverse);

template void vecRemoveIxs(std::vector<float> &vec, std::vector<size_t> &ixs);

template void vecRemoveIxs(std::vector<uint16_t> &vec, std::vector<size_t> &ixs);

template void vecRemoveIxs(std::vector<uint32_t> &vec, std::vector<size_t> &ixs);

template void vecRemoveIxs(std::vector<uint64_t> &vec, std::vector<size_t> &ixs);

template void arrRemoveIxs(uint16_t *arr, uint8_t last, std::vector<size_t> &ixs);

template void arrRemoveIxs(uint32_t *arr, uint8_t last, std::vector<size_t> &ixs);

template void arrRemoveIxs(uint64_t *arr, uint8_t last, std::vector<size_t> &ixs);

template void vecArgsort(std::vector<size_t> &ixs, const std::vector<uint16_t> &vec, bool reverse);

template void vecArgsort(std::vector<size_t> &ixs, const std::vector<uint32_t> &vec, bool reverse);

template void vecArgsort(std::vector<size_t> &ixs, const std::vector<uint64_t> &vec, bool reverse);

template void arrArgsort(std::vector<size_t> &ixs, const uint16_t *arr, uint8_t last, bool reverse);

template void arrArgsort(std::vector<size_t> &ixs, const uint32_t *arr, uint8_t last, bool reverse);

template void arrArgsort(std::vector<size_t> &ixs, const uint64_t *arr, uint8_t last, bool reverse);

template uint16_t histArgmax(const std::map<uint16_t, uint64_t> val_counts, uint64_t n, bool reverse);

template uint32_t histArgmax(const std::map<uint32_t, uint64_t> val_counts, uint64_t n, bool reverse);

template uint64_t histArgmax(const std::map<uint64_t, uint64_t> val_counts, uint64_t n, bool reverse);

template uint16_t vvecArgmax(const vvec<uint16_t> &vv, uint64_t n, bool reverse);

template uint32_t vvecArgmax(const vvec<uint32_t> &vv, uint64_t n, bool reverse);

template uint64_t vvecArgmax(const vvec<uint64_t> &vv, uint64_t n, bool reverse);
