#include "assess.h"

template<typename T>
void
vecRemoveIxs(std::vector<T>& v, std::vector<unsigned int>& ixs)
{
  ixs.push_back(-1);
  ixs.push_back(v.size());
  std::sort(ixs.begin(), ixs.end());
  typename std::vector<T>::iterator last = v.begin();
  for (unsigned int i = 1; i != ixs.size(); ++i) {
    unsigned int range_begin = ixs[i - 1] + 1;
    unsigned int range_end = ixs[i];
    std::copy(v.begin() + range_begin, v.begin() + range_end, last);
    last += range_end - range_begin;
  }
  v.erase(last, v.end());
}

template<typename T>
void
arrRemoveIxs(T* r, uint8_t last, std::vector<unsigned int>& ixs)
{
  std::sort(ixs.begin(), ixs.end());
  unsigned int ix = 0;
  if (!ixs.empty()) {
    std::vector<unsigned int>::iterator it = ixs.begin();
    for (unsigned int i = 0; i < last; ++i) {
      if (i == *it) {
        it++;
      } else {
        r[ix] = r[i];
        ix++;
      }
    }
  }
}

template<typename T>
void
vecArgsort1D(std::vector<unsigned int>& ixs, const std::vector<T>& v, bool reverse)
{
  ixs.resize(v.size());
  std::iota(ixs.begin(), ixs.end(), 0);
  if (reverse)
    std::sort(ixs.begin(), ixs.end(), [&v](unsigned int i1, unsigned int i2) { return v[i1] < v[i2]; });
  else
    std::sort(ixs.begin(), ixs.end(), [&v](unsigned int i1, unsigned int i2) { return v[i1] > v[i2]; });
}

template<typename T>
void
arrArgsort1D(std::vector<unsigned int>& ixs, const T* r, uint8_t last, bool reverse)
{
  ixs.resize(last);
  std::iota(ixs.begin(), ixs.end(), 0);
  if (reverse)
    std::sort(ixs.begin(), ixs.end(), [&r](unsigned int i1, unsigned int i2) { return r[i1] < r[i2]; });
  else
    std::sort(ixs.begin(), ixs.end(), [&r](unsigned int i1, unsigned int i2) { return r[i1] > r[i2]; });
}

template<typename T>
T
vvecArgmax2D(const vvec<T>& vv, uint64_t n, bool reverse)
{
  std::map<T, uint64_t> val_counts;
  for (auto& v : vv) {
    for (auto val : v)
      val_counts[val]++;
  }
  uint64_t cumsum = 0;
  T val_prev;
  if (reverse) {
    for (auto it = val_counts.cbegin(); (it != val_counts.cend()) && (cumsum < n); ++it) {
      val_prev = 0;
      cumsum += it->second;
      val_prev = it->first;
    }
  } else {
    for (auto it = val_counts.crbegin(); (it != val_counts.crend()) && (cumsum < n); ++it) {
      val_prev = std::numeric_limits<T>::max();
      cumsum += it->second;
      val_prev = it->first;
    }
  }
  return val_prev;
}

template<typename T>
T
arrArgmax2D(const T* r, const uint8_t* ind_r, uint32_t num_rows, uint8_t b, uint64_t n, bool reverse)
{
  std::map<T, uint64_t> val_counts;
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    for (uint8_t i = 0; i < ind_r[rix]; ++i)
      val_counts[r[b * rix + i]]++;
  }
  uint64_t cumsum = 0;
  T val_prev;
  if (reverse) {
    val_prev = 0;
    for (auto it = val_counts.cbegin(); (it != val_counts.cend()) && (cumsum < n); ++it) {
      cumsum += it->second;
      val_prev = it->first;
    }
  } else {
    for (auto it = val_counts.crbegin(); (it != val_counts.crend()) && (cumsum < n); ++it) {
      val_prev = std::numeric_limits<T>::max();
      cumsum += it->second;
      val_prev = it->first;
    }
  }
  return val_prev;
}

void
getIxsRandom(std::vector<unsigned int>& ixs, uint8_t b_curr, uint8_t b_max)
{
  if (ixs.size() != b_curr)
    ixs.resize(b_curr);
  std::iota(ixs.begin(), ixs.end(), 0);
  std::random_shuffle(ixs.begin(), ixs.end());
  ixs.resize(b_max);
  std::sort(ixs.begin(), ixs.end());
}

template<typename T>
void
vecIxsNumber(std::vector<unsigned int>& ixs, const std::vector<T>& s_v, uint8_t number)
{
  vecArgsort1D(ixs, s_v, true);
  ixs.resize(number);
  std::sort(ixs.begin(), ixs.end());
}

template<typename T>
void
vecIxsThreshold(std::vector<unsigned int>& ixs, const std::vector<T>& s_v, T threshold)
{
  ixs.clear();
  for (unsigned int ix = 0; ix < s_v.size(); ++ix) {
    if (s_v[ix] > threshold)
      ixs.push_back(ix);
  }
}

template<typename T>
void
arrIxsNumber(std::vector<unsigned int>& ixs, const T* s_r, uint8_t number)
{
  arrArgsort1D(ixs, s_r, true);
  ixs.resize(number);
  std::sort(ixs.begin(), ixs.end());
}

template<typename T>
void
arrIxsThreshold(std::vector<unsigned int>& ixs, const T* s_r, uint8_t last, T threshold)
{
  ixs.clear();
  for (unsigned int ix = 0; ix < (unsigned int)last; ++ix) {
    if (s_r[ix] > threshold)
      ixs.push_back(ix);
  }
}

template<typename T>
void
vvecSizeOrder(std::vector<unsigned int>& ixs, const vvec<T>& vv, bool reverse)
{
  ixs.resize(vv.size());
  std::iota(ixs.begin(), ixs.end(), 0);
  if (reverse)
    std::sort(ixs.begin(), ixs.end(), [&vv](unsigned int i1, unsigned int i2) { return vv[i1].size() < vv[i2].size(); });
  else
    std::sort(ixs.begin(), ixs.end(), [&vv](unsigned int i1, unsigned int i2) { return vv[i1].size() > vv[i2].size(); });
}

#include "assessins.cpp"
