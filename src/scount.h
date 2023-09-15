#ifndef _SCOUNTS_H
#define _SCOUNTS_H

#include <cstdint>
#pragma GCC optimize("Og,inline")

template<typename T>
inline std::unordered_map<T, scT>
mapTotalCountsHTd(std::vector<T>& r1,
                  std::vector<T>& r2,
                  std::vector<scT>& c1,
                  std::vector<scT>& c2)
{
  std::unordered_map<T, scT> count_map{};
  for (unsigned int i = 0; i < r1.size(); ++i)
    count_map[r1[i]] = c1[i];
  std::unordered_map<T, scT> inc_map{};
  for (unsigned int i = 0; i < r2.size(); ++i)
    inc_map[r2[i]] = c2[i];
  count_map = std::accumulate(inc_map.begin(),
                              inc_map.end(),
                              count_map,
                              [](std::unordered_map<T, scT>& m, const std::pair<const T, scT>& p) {
                                return (m[p.first] += p.second, m);
                              });
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
  count_map = std::accumulate(inc_map.begin(),
                              inc_map.end(),
                              count_map,
                              [](std::unordered_map<T, scT>& m, const std::pair<const T, scT>& p) {
                                return (m[p.first] += p.second, m);
                              });
  return count_map;
}

template<typename T>
inline std::unordered_map<T, std::vector<float>>
mapValuesProbabilityHTd(std::vector<HTd<T>>& td_vec, std::uint32_t rix)
{
  std::unordered_map<T, std::vector<float>> prob_map{};
  for (auto& td : td_vec) {
    for (unsigned int i = 0; i < td.enc_vvec[rix].size(); ++i)
      prob_map[td.enc_vvec[rix][i]].push_back(static_cast<float>(td.scount_vvec[rix][i]) /
                                              static_cast<float>(td.num_species));
  }
  for (auto& kv : prob_map) {
    for (unsigned int i = 0; i < td_vec.size() - kv.second.size(); ++i)
      kv.second.push_back(0);
  }
  return prob_map;
}

template<typename T>
inline std::unordered_map<T, std::vector<float>>
mapValuesProbabilityHTs(std::vector<HTs<T>>& ts_vec, std::uint32_t rix)
{
  std::unordered_map<T, std::vector<float>> prob_map{};
  for (auto& ts : ts_vec) {
    for (unsigned int i = 0; i < ts.ind_arr[rix]; ++i)
      prob_map[ts.enc_arr[rix * ts.b + i]].push_back(
        static_cast<float>(ts.scount_arr[rix * ts.b + i]) / static_cast<float>(ts.num_species));
  }
  for (auto& kv : prob_map) {
    for (unsigned int i = 0; i < ts_vec.size() - kv.second.size(); ++i)
      kv.second.push_back(0);
  }
  return prob_map;
}

template<typename T>
inline std::unordered_map<T, std::vector<float>>
mapValuesSelectHTd(std::vector<HTd<T>>& td_vec, std::uint32_t rix)
{
  std::unordered_map<T, std::vector<float>> select_map{};
  uint64_t total_num_species = 0;
  for (auto& td : td_vec)
    total_num_species += td.num_species;
  for (auto& td : td_vec) {
    for (unsigned int i = 0; i < td.enc_vvec[rix].size(); ++i)
      select_map[td.enc_vvec[rix][i]].push_back(static_cast<float>(td.num_species) /
                                                static_cast<float>(total_num_species));
  }
  for (auto& kv : select_map) {
    for (unsigned int i = 0; i < td_vec.size() - kv.second.size(); ++i)
      kv.second.push_back(0);
  }
  return select_map;
}

template<typename T>
inline std::unordered_map<T, std::vector<float>>
mapValuesSelectHTs(std::vector<HTs<T>>& ts_vec, std::uint32_t rix)
{
  std::unordered_map<T, std::vector<float>> select_map{};
  uint64_t total_num_species = 0;
  for (auto& ts : ts_vec)
    total_num_species += ts.num_species;
  for (auto& ts : ts_vec) {
    for (unsigned int i = 0; i < ts.ind_arr[rix]; ++i)
      select_map[ts.enc_arr[rix * ts.b + i]].push_back(static_cast<float>(ts.num_species) /
                                                       static_cast<float>(total_num_species));
  }
  for (auto& kv : select_map) {
    for (unsigned int i = 0; i < ts_vec.size() - kv.second.size(); ++i)
      kv.second.push_back(0);
  }
  return select_map;
}

template<typename T>
inline std::unordered_map<T, std::vector<tT>>
mapValuesLCAtHTd(std::vector<HTd<T>>& td_vec, std::uint32_t rix)
{
  std::unordered_map<T, std::vector<tT>> tlca_map{};
  for (auto& td : td_vec) {
    for (unsigned int i = 0; i < td.enc_vvec[rix].size(); ++i) {
      tlca_map[td.enc_vvec[rix][i]].push_back(td.tlca_vvec[rix][i]);
    }
  }
  return tlca_map;
}

template<typename T>
inline std::unordered_map<T, std::vector<tT>>
mapValuesLCAtHTs(std::vector<HTs<T>>& ts_vec, std::uint32_t rix)
{
  std::unordered_map<T, std::vector<tT>> tlca_map{};
  for (auto& ts : ts_vec) {
    for (unsigned int i = 0; i < ts.ind_arr[rix]; ++i)
      tlca_map[ts.enc_arr[rix * ts.b + i]].push_back(ts.tlca_arr[rix * ts.b + i]);
  }
  return tlca_map;
}

template<typename T>
inline std::unordered_map<T, std::vector<scT>>
mapValuesCountsHTd(std::vector<HTd<T>>& td_vec, std::uint32_t rix)
{
  std::unordered_map<T, std::vector<scT>> values_map{};
  for (auto& td : td_vec) {
    for (unsigned int i = 0; i < td.enc_vvec[rix].size(); ++i)
      values_map[td.enc_vvec[rix][i]].push_back(td.scount_vvec[rix][i]);
  }
  for (auto& kv : values_map) {
    for (unsigned int i = 0; i < td_vec.size() - kv.second.size(); ++i)
      kv.second.push_back(0);
  }
  return values_map;
}

template<typename T>
inline std::unordered_map<T, std::vector<bool>>
mapValuesBinaryHTd(std::vector<HTd<T>>& td_vec, std::uint32_t rix)
{
  std::unordered_map<T, std::vector<bool>> values_map{};
  for (auto& td : td_vec) {
    for (unsigned int i = 0; i < td.enc_vvec[rix].size(); ++i)
      values_map[td.enc_vvec[rix][i]].push_back(true);
  }
  for (auto& kv : values_map) {
    for (unsigned int i = 0; i < td_vec.size() - kv.second.size(); ++i)
      kv.second.push_back(false);
  }
  return values_map;
}

template<typename T>
inline std::unordered_map<T, std::vector<scT>>
mapValuesCountsHTs(std::vector<HTs<T>>& ts_vec, std::uint32_t rix)
{
  std::unordered_map<T, std::vector<scT>> values_map{};
  for (auto& ts : ts_vec) {
    for (unsigned int i = 0; i < ts.ind_arr[rix]; ++i)
      values_map[ts.enc_arr[rix * ts.b + i]].push_back(ts.scount_arr[rix * ts.b + i]);
  }
  for (auto& kv : values_map) {
    for (unsigned int i = 0; i < ts_vec.size() - kv.second.size(); ++i)
      kv.second.push_back(0);
  }
  return values_map;
}

template<typename T>
inline std::unordered_map<T, std::vector<bool>>
mapValuesBinaryHTs(std::vector<HTs<T>>& ts_vec, std::uint32_t rix)
{
  std::unordered_map<T, std::vector<bool>> values_map{};
  for (auto& ts : ts_vec) {
    for (unsigned int i = 0; i < ts.ind_arr[rix]; ++i)
      values_map[ts.enc_arr[rix * ts.b + i]].push_back(true);
  }
  for (auto& kv : values_map) {
    for (unsigned int i = 0; i < ts_vec.size() - kv.second.size(); ++i)
      kv.second.push_back(false);
  }
  return values_map;
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
inline std::unordered_map<T, tT>
mapLCAtHTd(std::vector<T>& r, std::vector<tT>& t)
{
  std::unordered_map<T, tT> tlca_map{};
  for (unsigned int i = 0; i < r.size(); ++i) {
    tlca_map[r[i]] = t[i];
  }
  return tlca_map;
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
inline std::unordered_map<T, tT>
mapLCAtHTs(T* r, tT* t, uint8_t last)
{
  std::unordered_map<T, tT> tlca_map{};
  for (unsigned int i = 0; i < last; ++i) {
    tlca_map[r[i]] = t[i];
  }
  return tlca_map;
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
updateLCAtHTd(std::vector<T>& r, std::vector<tT>& t, std::unordered_map<T, tT> tlca_map)
{
  t.clear();
  t.resize(r.size());
  for (unsigned int i = 0; i < r.size(); ++i) {
    t[i] = tlca_map[r[i]];
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

template<typename T>
inline void
updateLCAtHTs(T* r, tT* t, uint8_t last, std::unordered_map<T, tT> tlca_map)
{
  for (unsigned int i = 0; i < last; ++i) {
    t[i] = tlca_map[r[i]];
  }
}

#endif
