#ifndef _SCOUNTS_H
#define _SCOUNTS_H

#pragma GCC optimize("Og,inline")

template<typename T>
inline std::unordered_map<T, scT>
mapTotalCountsHTd(std::vector<T>& r1, std::vector<T>& r2, std::vector<scT>& c1, std::vector<scT>& c2)
{
  std::unordered_map<T, scT> count_map{};
  for (unsigned int i = 0; i < r1.size(); ++i)
    count_map[r1[i]] = c1[i];
  std::unordered_map<T, scT> inc_map{};
  for (unsigned int i = 0; i < r2.size(); ++i)
    inc_map[r2[i]] = c2[i];
  count_map = std::accumulate(
    inc_map.begin(), inc_map.end(), count_map, [](std::unordered_map<T, scT>& m, const std::pair<const T, scT>& p) { return (m[p.first] += p.second, m); });
  return count_map;
}

template<typename T>
inline std::unordered_map<T, scT>
mapTotalCountsHTs(T* r1, scT* c1, uint8_t last1, T* r2, scT* c2, uint8_t last2)
{
  std::unordered_map<T, scT> count_map{};
  for (unsigned int i = 0; i < last1; ++i)
    count_map[r1[i]] = c1[i];
  std::unordered_map<T, scT> inc_map{};
  for (unsigned int i = 0; i < last2; ++i)
    inc_map[r2[i]] = c2[i];
  count_map = std::accumulate(
    inc_map.begin(), inc_map.end(), count_map, [](std::unordered_map<T, scT>& m, const std::pair<const T, scT>& p) { return (m[p.first] += p.second, m); });
  return count_map;
}

template<typename T>
inline std::unordered_map<T, scT>
mapCountsHTd(std::vector<T>& r, std::vector<scT>& c)
{
  std::unordered_map<T, scT> count_map{};
  for (unsigned int i = 0; i < r.size(); ++i) {
    count_map[r[i]] = c[i];
  }
  return count_map;
}

template<typename T>
inline std::unordered_map<T, scT>
mapCountsHTs(T* r, scT* c, uint8_t last)
{
  std::unordered_map<T, scT> count_map{};
  for (unsigned int i = 0; i < last; ++i) {
    count_map[r[i]] = c[i];
  }
  return count_map;
}

template<typename T>
inline void
updateCountsHTd(std::vector<T>& r, std::vector<scT>& c, std::unordered_map<T, scT> count_map)
{
  c.clear();
  c.resize(r.size());
  for (unsigned int i = 0; i < r.size(); ++i) {
    c[i] = count_map[r[i]];
  }
}

template<typename T>
inline void
updateCountsHTs(T* r, scT* c, uint8_t last, std::unordered_map<T, scT> count_map)
{
  for (unsigned int i = 0; i < last; ++i) {
    c[i] = count_map[r[i]];
  }
}

#endif
