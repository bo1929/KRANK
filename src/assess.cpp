#include "assess.h"
#include <ostream>

template<typename T1, typename T2>
void
vecInformationScores(std::vector<T1>& s,
                     std::vector<T2>& v,
                     std::unordered_map<T2, std::vector<T1>>& values_map)
{
  s.resize(v.size());
  for (unsigned int i = 0; i < v.size(); ++i) {
    std::sort(values_map[v[i]].begin(), values_map[v[i]].end());
    if (values_map[v[i]].size() > 1)
      s[i] = 1 + (3 * values_map[v[i]].back()) - (2 * values_map[v[i]].end()[-2]) -
             (values_map[v[i]][0]);
    else
      s[i] = 1;
  }
}

template<typename T1, typename T2>
void
arrInformationScores(std::vector<T1>& s,
                     T2* r,
                     uint8_t last,
                     std::unordered_map<T2, std::vector<T1>>& values_map)
{
  s.resize(last);
  for (unsigned int i = 0; i < last; ++i) {
    std::sort(values_map[r[i]].begin(), values_map[r[i]].end());
    if (values_map[r[i]].size() > 1)
      s[i] = 1 + (3 * values_map[r[i]].back()) - (2 * values_map[r[i]].end()[-2]) -
             (values_map[r[i]][0]);
    else
      s[i] = 1;
  }
}

template<typename T>
void
vecRemoveIxs(std::vector<T>& v, std::vector<unsigned int>& ixs)
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
void
arrRemoveIxs(T* r, uint8_t last, std::vector<unsigned int>& ixs)
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
void
vecArgsort1D(std::vector<unsigned int>& ixs, const std::vector<T>& v, bool reverse)
{
  ixs.resize(v.size());
  std::iota(ixs.begin(), ixs.end(), 0);
  std::random_shuffle(ixs.begin(), ixs.end());
  if (reverse)
    std::sort(
      ixs.begin(), ixs.end(), [&v](unsigned int i1, unsigned int i2) { return v[i1] > v[i2]; });
  else
    std::sort(
      ixs.begin(), ixs.end(), [&v](unsigned int i1, unsigned int i2) { return v[i1] < v[i2]; });
}

template<typename T>
void
arrArgsort1D(std::vector<unsigned int>& ixs, const T* r, uint8_t last, bool reverse)
{
  ixs.resize(last);
  std::iota(ixs.begin(), ixs.end(), 0);
  std::random_shuffle(ixs.begin(), ixs.end());
  if (reverse)
    std::sort(
      ixs.begin(), ixs.end(), [&r](unsigned int i1, unsigned int i2) { return r[i1] > r[i2]; });
  else
    std::sort(
      ixs.begin(), ixs.end(), [&r](unsigned int i1, unsigned int i2) { return r[i1] < r[i2]; });
}

template<typename T>
T
mapArgmax(std::map<T, uint64_t>& val_counts, uint64_t n, bool reverse)
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

template<typename T>
T
arrArgmax2D(const T* r,
            const uint8_t* ind_r,
            uint32_t num_rows,
            uint8_t b,
            uint64_t n,
            bool reverse)
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

void
getIxsRandom(std::vector<unsigned int>& ixs, uint8_t b, uint8_t number)
{
  if (number > b)
    std::cout << std::endl << (int)number << "," << (int)b << std::endl;
  assert(b > 0);
  assert(number <= b);
  if (ixs.size() != b)
    ixs.resize(b);
  std::iota(ixs.begin(), ixs.end(), 0);
  std::random_shuffle(ixs.begin(), ixs.end());
  ixs.resize(number);
  std::sort(ixs.begin(), ixs.end());
}

template<typename T>
void
vecIxsNumber(std::vector<unsigned int>& ixs,
             const std::vector<T>& s_v,
             uint8_t number,
             bool reverse)
{
  vecArgsort1D(ixs, s_v, reverse);
  ixs.resize(number);
  std::sort(ixs.begin(), ixs.end());
}

template<typename T>
void
vecIxsThreshold(std::vector<unsigned int>& ixs,
                const std::vector<T>& s_v,
                T threshold,
                bool reverse)
{
  ixs.clear();
  for (unsigned int ix = 0; ix < s_v.size(); ++ix) {
    if (reverse) {
      if (s_v[ix] <= threshold)
        ixs.push_back(ix);
    } else {
      if (s_v[ix] >= threshold)
        ixs.push_back(ix);
    }
    /* if ((s_v.size() - ixs.size() <= 2)) */
    /*   break; */
  }
}

template<typename T>
void
arrIxsNumber(std::vector<unsigned int>& ixs, const T* s_r, uint8_t number, bool reverse)
{
  arrArgsort1D(ixs, s_r, reverse);
  ixs.resize(number);
  std::sort(ixs.begin(), ixs.end());
}

template<typename T>
void
arrIxsThreshold(std::vector<unsigned int>& ixs,
                const T* s_r,
                uint8_t last,
                T threshold,
                bool reverse)
{
  ixs.clear();
  for (unsigned int ix = 0; ix < last; ++ix) {
    if (reverse) {
      if (s_r[ix] <= threshold)
        ixs.push_back(ix);
    } else {
      if (s_r[ix] >= threshold)
        ixs.push_back(ix);
    }
    /* if ((last - ixs.size() <= 2)) */
    /*   break; */
  }
}

template<typename T>
void
vvecSizeOrder(std::vector<unsigned int>& ixs, const vvec<T>& vv, bool reverse)
{
  ixs.resize(vv.size());
  std::iota(ixs.begin(), ixs.end(), 0);
  if (reverse)
    std::sort(ixs.begin(), ixs.end(), [&vv](unsigned int i1, unsigned int i2) {
      return vv[i1].size() > vv[i2].size();
    });
  else
    std::sort(ixs.begin(), ixs.end(), [&vv](unsigned int i1, unsigned int i2) {
      return vv[i1].size() < vv[i2].size();
    });
}

void
arrSizeOrder(std::vector<unsigned int>& ixs, const uint8_t* ind_r, uint32_t num_rows, bool reverse)
{
  ixs.resize(num_rows);
  std::iota(ixs.begin(), ixs.end(), 0);
  if (reverse)
    std::sort(ixs.begin(), ixs.end(), [&ind_r](unsigned int i1, unsigned int i2) {
      return ind_r[i1] > ind_r[i2];
    });
  else
    std::sort(ixs.begin(), ixs.end(), [&ind_r](unsigned int i1, unsigned int i2) {
      return ind_r[i1] < ind_r[i2];
    });
}

#include "assessins.cpp"
