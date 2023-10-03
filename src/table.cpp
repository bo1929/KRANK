#include "table.h"

#undef DEBUG
#define FL true

template<typename encT>
std::unordered_map<uint8_t, uint64_t>
HTd<encT>::histRowSizes()
{
  std::unordered_map<uint8_t, uint64_t> hist_map;
  for (auto& row : enc_vvec) {
    hist_map[row.size()]++;
  }
  return hist_map;
}

template<typename encT>
std::unordered_map<uint8_t, uint64_t>
HTs<encT>::histRowSizes()
{
  std::unordered_map<uint8_t, uint64_t> hist_map;
  for (unsigned int rix = 0; rix < num_rows; ++rix) {
    hist_map[ind_arr[rix]]++;
  }
  return hist_map;
}

template<typename encT>
void
HTd<encT>::makeUnique(bool update_size)
{
#ifdef DEBUG
  if (!HTd<encT>::areColumnsSorted()) {
    HTd<encT>::sortColumns();
    std::puts("HTd has unsorted columns before making columns unique.\n");
  }
#endif
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!enc_vvec[rix].empty()) {
      std::unordered_map<encT, scT> count_map = mapCountsHTd(enc_vvec[rix], scount_vvec[rix]);
      std::unordered_map<encT, tT> tlca_map = mapLCAtHTd(enc_vvec[rix], tlca_vvec[rix]);
      enc_vvec[rix].erase(std::unique(enc_vvec[rix].begin(), enc_vvec[rix].end()),
                          enc_vvec[rix].end());
      updateCountsHTd(enc_vvec[rix], scount_vvec[rix], count_map);
      updateLCAtHTd(enc_vvec[rix], tlca_vvec[rix], tlca_map);
    }
  }
  if (update_size)
    HTd<encT>::updateSize();
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void
HTs<encT>::makeUnique(bool update_size)
{
#ifdef DEBUG
  if (!HTs<encT>::areColumnsSorted()) {
    HTs<encT>::sortColumns();
    std::puts("HTs has unsorted columns before making columns unique.\n");
  }
#endif
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (ind_arr[rix] > 0) {
      std::unordered_map<encT, scT> count_map =
        mapCountsHTs(enc_arr + (rix * b), scount_arr + (rix * b), ind_arr[rix]);
      std::unordered_map<encT, tT> tlca_map =
        mapLCAtHTs(enc_arr + (rix * b), tlca_arr + (rix * b), ind_arr[rix]);
      auto last =
        std::unique(enc_arr + (rix * b), enc_arr + (rix * b) + ind_arr[rix]) - enc_arr - (rix * b);
      ind_arr[rix] = last;
      updateCountsHTs(enc_arr + (rix * b), scount_arr + (rix * b), ind_arr[rix], count_map);
      updateLCAtHTs(enc_arr + (rix * b), tlca_arr + (rix * b), ind_arr[rix], tlca_map);
    }
  }
  if (update_size)
    HTs<encT>::updateSize();
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void
HTd<encT>::sortColumns()
{
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!enc_vvec[rix].empty()) {
      std::unordered_map<encT, scT> count_map = mapCountsHTd(enc_vvec[rix], scount_vvec[rix]);
      std::unordered_map<encT, tT> tlca_map = mapLCAtHTd(enc_vvec[rix], tlca_vvec[rix]);
      std::sort(enc_vvec[rix].begin(), enc_vvec[rix].end());
      updateCountsHTd(enc_vvec[rix], scount_vvec[rix], count_map);
      updateLCAtHTd(enc_vvec[rix], tlca_vvec[rix], tlca_map);
    }
  }
}

template<typename encT>
void
HTs<encT>::sortColumns()
{
  uint8_t b = HTs<encT>::b;
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (ind_arr[rix] > 0) {
      std::unordered_map<encT, scT> count_map =
        mapCountsHTs(enc_arr + (rix * b), scount_arr + (rix * b), ind_arr[rix]);
      std::unordered_map<encT, tT> tlca_map =
        mapLCAtHTs(enc_arr + (rix * b), tlca_arr + (rix * b), ind_arr[rix]);
      std::sort(enc_arr + (rix * b), enc_arr + (rix * b) + ind_arr[rix]);
      updateCountsHTs(enc_arr + (rix * b), scount_arr + (rix * b), ind_arr[rix], count_map);
      updateLCAtHTs(enc_arr + (rix * b), tlca_arr + (rix * b), ind_arr[rix], tlca_map);
    }
  }
}

template<typename encT>
bool
HTd<encT>::areColumnsSorted()
{
  bool is_ok = true;
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!enc_vvec[rix].empty() && !std::is_sorted(enc_vvec[rix].begin(), enc_vvec[rix].end())) {
#pragma omp atomic write
      is_ok = false;
    }
  }
  return is_ok;
}

template<typename encT>
bool
HTs<encT>::areColumnsSorted()
{
  bool is_ok = true;
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if ((ind_arr[rix] > 0) &&
        !std::is_sorted(enc_arr + (rix * b), enc_arr + (rix * b) + ind_arr[rix])) {
#pragma omp atomic write
      is_ok = false;
    }
  }
  return is_ok;
}

template<typename encT>
void
HTd<encT>::clearRows()
{
#ifdef DEBUG
  LOG(INFO) << "The current size, before clearing, of the table is " << num_kmers << std::endl;
#endif
  enc_vvec.clear();
  enc_vvec.resize(num_rows);
  scount_vvec.clear();
  scount_vvec.resize(num_rows);
  tlca_vvec.clear();
  tlca_vvec.resize(num_rows);
  HTd<encT>::updateSize();
  HTd<encT>::num_species = 0;
  HTd<encT>::tIDsBasis = {};
#ifdef DEBUG
  LOG(INFO) << "The current size, after clearing, of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void
HTs<encT>::clearRows()
{
#ifdef DEBUG
  LOG(INFO) << "The current size, before clearing, of the table is " << num_kmers << std::endl;
#endif
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    ind_arr[rix] = 0;
  }
  HTs<encT>::updateSize();
  num_species = 0;
  tIDsBasis = {};
#ifdef DEBUG
  LOG(INFO) << "The current size, after clearing, of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void
HTd<encT>::updateSize()
{
  uint64_t num_kmers_sum = 0;
  /* #pragma omp parallel for num_threads(num_threads), schedule(dynamic) */
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    /* if (!enc_vvec[rix].empty()) { */
    /* #pragma omp atomic update */
    num_kmers_sum += enc_vvec[rix].size();
    /* } */
  }
  num_kmers = num_kmers_sum;
  table_size = num_kmers * sizeof(encT) + num_rows * sizeof(std::vector<encT>);
  table_size += num_kmers * sizeof(scT) + num_rows * sizeof(std::vector<scT>);
  table_size += num_kmers * sizeof(tT) + num_rows * sizeof(std::vector<tT>);
}

template<typename encT>
void
HTs<encT>::updateSize()
{
  uint64_t num_kmers_sum = 0;
  /* #pragma omp parallel for num_threads(num_threads), schedule(dynamic) */
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    /* #pragma omp atomic update */
    num_kmers_sum += ind_arr[rix];
  }
  num_kmers = num_kmers_sum;
  table_size = table_size;
}

template<typename encT>
void
HTd<encT>::initBasis(tT tID)
{
  HTd<encT>::updateSize();
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
  uint32_t value;
  if (num_kmers > 0)
    value = 1;
  else
    value = 0;
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!enc_vvec[rix].empty()) {
      scount_vvec[rix].resize(enc_vvec[rix].size());
      std::fill(scount_vvec[rix].begin(), scount_vvec[rix].end(), value);
      tlca_vvec[rix].resize(enc_vvec[rix].size());
      std::fill(tlca_vvec[rix].begin(), tlca_vvec[rix].end(), tID);
    } else {
      scount_vvec[rix].clear();
      tlca_vvec[rix].clear();
    }
  }
  num_species = value;
  tIDsBasis = { tID };
}

template<typename encT>
void
HTd<encT>::trimColumns(uint8_t b_max)
{
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (enc_vvec[rix].size() > b_max) {
      enc_vvec[rix].resize(b_max);
      scount_vvec[rix].resize(b_max);
      tlca_vvec[rix].resize(b_max);
    }
  }
#ifdef DEBUG
  HTd<encT>::updateSize();
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void
HTd<encT>::pruneColumns(uint8_t b_max)
{
#ifdef DEBUG
  LOG(INFO) << "The current size, before pruning, of the table is " << num_kmers << std::endl;
#endif
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (enc_vvec[rix].size() >= std::numeric_limits<uint8_t>::max()) {
      enc_vvec[rix].resize(std::numeric_limits<uint8_t>::max() - 1);
      scount_vvec[rix].resize(std::numeric_limits<uint8_t>::max() - 1);
      tlca_vvec[rix].resize(std::numeric_limits<uint8_t>::max() - 1);
    }
    if (enc_vvec[rix].size() > b_max) {
      std::vector<unsigned int> ixs;
      switch (ranking_method) {
        case random_kmer: {
          getIxsRandom(ixs, enc_vvec[rix].size(), enc_vvec[rix].size() - b_max);
          break;
        }
        case species_count: {
          vecIxsNumber(ixs, scount_vvec[rix], enc_vvec[rix].size() - b_max, FL);
          break;
        }
        case information_score: {
          std::unordered_map<encT, std::vector<scT>> values_map =
            mapValuesCountsHTd(childrenHT, rix);
          std::vector<scT> scores_vec;
          vecInformationScores(scores_vec, enc_vvec[rix], values_map);
          vecIxsNumber(ixs, scores_vec, enc_vvec[rix].size() - b_max, FL);
          break;
        }
        case taxa_count: {
          std::unordered_map<encT, std::vector<bool>> values_map =
            mapValuesBinaryHTd(childrenHT, rix);
          std::vector<scT> tcount_vec;
          vecTaxaCounts(tcount_vec, enc_vvec[rix], values_map);
          vecIxsNumber(ixs, tcount_vec, enc_vvec[rix].size() - b_max, FL);
          break;
        }
      }
      vecRemoveIxs(enc_vvec[rix], ixs);
      vecRemoveIxs(scount_vvec[rix], ixs);
      vecRemoveIxs(tlca_vvec[rix], ixs);
    }
  }
#ifdef DEBUG
  HTd<encT>::updateSize();
  if (!HTd<encT>::areColumnsSorted()) {
    HTd<encT>::sortColumns();
    std::puts("HTd has unsorted columns after pruning.\n");
  }
  LOG(INFO) << "The current size, after pruning, of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void
HTd<encT>::unionRows(HTd<encT>& sibling, bool update_size)
{
#ifdef DEBUG
  if (!HTd<encT>::areColumnsSorted()) {
    HTd<encT>::sortColumns();
    std::puts("HTd (main) has unsorted columns before computing the union of "
              "rows.\n");
  }
  if (!sibling.areColumnsSorted()) {
    sibling.sortColumns();
    std::puts("HTd (sibling) has unsorted columns before computing the union "
              "of rows.\n");
  }
#endif
  assert(num_rows == sibling.num_rows);
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!enc_vvec[rix].empty() && !sibling.enc_vvec[rix].empty()) {
      std::unordered_map<encT, scT> count_map = mapTotalCountsHTd(
        enc_vvec[rix], sibling.enc_vvec[rix], scount_vvec[rix], sibling.scount_vvec[rix]);
      enc_vvec[rix].insert(
        enc_vvec[rix].end(), sibling.enc_vvec[rix].begin(), sibling.enc_vvec[rix].end());
      std::inplace_merge(enc_vvec[rix].begin(),
                         enc_vvec[rix].begin() +
                           (enc_vvec[rix].size() - sibling.enc_vvec[rix].size()),
                         enc_vvec[rix].end());
      enc_vvec[rix].erase(std::unique(enc_vvec[rix].begin(), enc_vvec[rix].end()),
                          enc_vvec[rix].end());
      updateCountsHTd(enc_vvec[rix], scount_vvec[rix], count_map);
    } else if (!sibling.enc_vvec[rix].empty()) {
      enc_vvec[rix] = sibling.enc_vvec[rix];
      scount_vvec[rix] = sibling.scount_vvec[rix];
      tlca_vvec[rix] = sibling.tlca_vvec[rix];
    }
  }
  num_species += sibling.num_species;
  tIDsBasis.insert(sibling.tIDsBasis.begin(), sibling.tIDsBasis.end());
  HTd<encT>::pruneColumns(std::numeric_limits<uint8_t>::max());
  if (update_size)
    HTd<encT>::updateSize();
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void
HTs<encT>::unionRows(HTs<encT>& sibling, bool update_size)
{
#ifdef DEBUG
  if (!HTs<encT>::areColumnsSorted()) {
    HTs<encT>::sortColumns();
    std::puts("HTs (main) has unsorted columns before computing the union of "
              "rows.\n");
  }
  if (!sibling.areColumnsSorted()) {
    sibling.sortColumns();
    std::puts("HTs (sibling) has unsorted columns before computing the union "
              "of rows.\n");
  }
#endif
  assert(b == sibling.b);
  assert(num_rows * b == sibling.num_rows * sibling.b);
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    unsigned int ix = rix * b;
    if ((sibling.ind_arr[rix] > 0) && (ind_arr[rix] > 0)) {
      std::unordered_map<encT, scT> count_map = mapTotalCountsHTs(enc_arr + ix,
                                                                  scount_arr + ix,
                                                                  ind_arr[rix],
                                                                  sibling.enc_arr + ix,
                                                                  sibling.scount_arr + ix,
                                                                  sibling.ind_arr[rix]);
      std::vector<encT> v;
      std::set_union(enc_arr + ix,
                     enc_arr + (ix + ind_arr[rix]),
                     sibling.enc_arr + ix,
                     sibling.enc_arr + (ix + sibling.ind_arr[rix]),
                     std::back_inserter(v));
      if (v.size() > b) {
        std::vector<unsigned int> ixs;
        switch (ranking_method) {
          case random_kmer: {
            getIxsRandom(ixs, v.size(), v.size() - b);
            break;
          }
          case species_count: {
            std::vector<scT> s_v;
            std::transform(v.begin(), v.end(), back_inserter(s_v), [&count_map](encT val) {
              return count_map[val];
            });
            vecIxsNumber(ixs, s_v, v.size() - b, FL);
            break;
          }
          case information_score: {
            std::unordered_map<encT, std::vector<scT>> values_map =
              mapValuesCountsHTs(childrenHT, rix);
            std::vector<scT> scores_vec;
            vecInformationScores(scores_vec, v, values_map);
            vecIxsNumber(ixs, scores_vec, v.size() - b, FL);
            break;
          }
          case taxa_count: {
            std::unordered_map<encT, std::vector<bool>> values_map =
              mapValuesBinaryHTs(childrenHT, rix);
            std::vector<scT> tcount_vec;
            vecTaxaCounts(tcount_vec, v, values_map);
            vecIxsNumber(ixs, tcount_vec, v.size() - b, FL);
            break;
          }
        }
        vecRemoveIxs(v, ixs);
      }
      std::copy(v.begin(), v.end(), enc_arr + ix);
      ind_arr[rix] = v.size();
      updateCountsHTs(enc_arr + ix, scount_arr + ix, ind_arr[rix], count_map);
    } else if (sibling.ind_arr[rix] > 0) {
      std::copy(sibling.enc_arr + ix, sibling.enc_arr + (ix + sibling.ind_arr[rix]), enc_arr + ix);
      std::copy(
        sibling.scount_arr + ix, sibling.scount_arr + (ix + sibling.ind_arr[rix]), scount_arr + ix);
      std::copy(
        sibling.tlca_arr + ix, sibling.tlca_arr + (ix + sibling.ind_arr[rix]), tlca_arr + ix);
      ind_arr[rix] = sibling.ind_arr[rix];
    }
  }
  num_species += sibling.num_species;
  tIDsBasis.insert(sibling.tIDsBasis.begin(), sibling.tIDsBasis.end());
  if (update_size)
    HTs<encT>::updateSize();
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void
HTd<encT>::mergeRows(HTd<encT>& sibling, bool update_size)
{
#ifdef DEBUG
  if (!HTd<encT>::areColumnsSorted()) {
    HTd<encT>::sortColumns();
    std::puts("HTd (main) has unsorted columns before merging the rows.\n");
  }
  if (!sibling.areColumnsSorted()) {
    sibling.sortColumns();
    std::puts("HTd (sibling) has unsorted columns before merging the rows.\n");
  }
#endif
  assert(num_rows == sibling.num_rows);
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!enc_vvec[rix].empty() && !sibling.enc_vvec[rix].empty()) {
      std::unordered_map<encT, scT> count_map = mapTotalCountsHTd(
        enc_vvec[rix], sibling.enc_vvec[rix], scount_vvec[rix], sibling.scount_vvec[rix]);
      enc_vvec[rix].insert(
        enc_vvec[rix].end(), sibling.enc_vvec[rix].begin(), sibling.enc_vvec[rix].end());
      std::inplace_merge(enc_vvec[rix].begin(),
                         enc_vvec[rix].begin() +
                           (enc_vvec[rix].size() - sibling.enc_vvec[rix].size()),
                         enc_vvec[rix].end());
      updateCountsHTd(enc_vvec[rix], scount_vvec[rix], count_map);
    } else if (!sibling.enc_vvec[rix].empty()) {
      enc_vvec[rix] = sibling.enc_vvec[rix];
      scount_vvec[rix] = sibling.scount_vvec[rix];
      tlca_vvec[rix] = sibling.tlca_vvec[rix];
    }
  }
  num_species += sibling.num_species;
  tIDsBasis.insert(sibling.tIDsBasis.begin(), sibling.tIDsBasis.end());
  HTd<encT>::pruneColumns(std::numeric_limits<uint8_t>::max());
  if (update_size)
    HTd<encT>::updateSize();
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void
HTs<encT>::mergeRows(HTs<encT>& sibling, bool update_size)
{
#ifdef DEBUG
  if (!HTs<encT>::areColumnsSorted()) {
    HTs<encT>::sortColumns();
    std::puts("HTs (main) has unsorted columns before merging the rows.\n");
  }
  if (!sibling.areColumnsSorted()) {
    sibling.sortColumns();
    std::puts("HTs (sibling) has unsorted columns before merging the rows.\n");
  }
#endif
  assert(b == sibling.b);
  assert(num_rows * b == sibling.num_rows * sibling.b);
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    unsigned int ix = rix * b;
    if ((sibling.ind_arr[rix] > 0) && (ind_arr[rix] > 0)) {
      std::unordered_map<encT, scT> count_map = mapTotalCountsHTs(enc_arr + ix,
                                                                  scount_arr + ix,
                                                                  ind_arr[rix],
                                                                  sibling.enc_arr + ix,
                                                                  sibling.scount_arr + ix,
                                                                  sibling.ind_arr[rix]);
      std::vector<encT> v;
      std::merge(enc_arr + ix,
                 enc_arr + (ix + ind_arr[rix]),
                 sibling.enc_arr + ix,
                 sibling.enc_arr + (ix + sibling.ind_arr[rix]),
                 std::back_inserter(v));
      if (v.size() > b) {
        v.erase(std::unique(v.begin(), v.end()), v.end());
      }
      if (v.size() > b) {
        std::vector<unsigned int> ixs;
        switch (ranking_method) {
          case random_kmer: {
            getIxsRandom(ixs, v.size(), v.size() - b);
            break;
          }
          case species_count: {
            std::vector<scT> s_v;
            std::transform(v.begin(), v.end(), back_inserter(s_v), [&count_map](encT val) {
              return count_map[val];
            });
            vecIxsNumber(ixs, s_v, v.size() - b, FL);
            break;
          }
          case information_score: {
            std::unordered_map<encT, std::vector<scT>> values_map =
              mapValuesCountsHTs(childrenHT, rix);
            std::vector<scT> scores_vec;
            vecInformationScores(scores_vec, v, values_map);
            vecIxsNumber(ixs, scores_vec, v.size() - b, FL);
            break;
          }
          case taxa_count: {
            std::unordered_map<encT, std::vector<bool>> values_map =
              mapValuesBinaryHTs(childrenHT, rix);
            std::vector<scT> tcount_vec;
            vecTaxaCounts(tcount_vec, v, values_map);
            vecIxsNumber(ixs, tcount_vec, v.size() - b, FL);
            break;
          }
        }
        vecRemoveIxs(v, ixs);
      }
      std::copy(v.begin(), v.end(), enc_arr + ix);
      ind_arr[rix] = v.size();
      updateCountsHTs(enc_arr + ix, scount_arr + ix, ind_arr[rix], count_map);
    } else if (sibling.ind_arr[rix] > 0) {
      std::copy(sibling.enc_arr + ix, sibling.enc_arr + (ix + sibling.ind_arr[rix]), enc_arr + ix);
      std::copy(
        sibling.scount_arr + ix, sibling.scount_arr + (ix + sibling.ind_arr[rix]), scount_arr + ix);
      std::copy(
        sibling.tlca_arr + ix, sibling.tlca_arr + (ix + sibling.ind_arr[rix]), tlca_arr + ix);
      ind_arr[rix] = sibling.ind_arr[rix];
    }
  }
  num_species += sibling.num_species;
  tIDsBasis.insert(sibling.tIDsBasis.begin(), sibling.tIDsBasis.end());
  if (update_size)
    HTs<encT>::updateSize();
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void
HTd<encT>::shrinkHT(uint64_t num_rm, uint8_t b_max)
{
  assert(num_rm < num_kmers);
  uint64_t init_num_kmers = num_kmers;
  HTd<encT>::pruneColumns(b_max);
  HTd<encT>::updateSize();
  volatile int64_t to_rm = num_rm;
  to_rm = to_rm - (init_num_kmers - num_kmers);
  volatile int64_t pto_rm = to_rm;
#ifdef DEBUG
  LOG(INFO) << "Initial number of k-mers to remove is" << to_rm << std::endl;
#endif
  while (to_rm > 0) {
    assert(to_rm < std::numeric_limits<int64_t>::max());
    assert(to_rm < num_kmers);
    std::vector<unsigned int> row_order;
    vvecSizeOrder(row_order, scount_vvec, true);
    uint8_t n = static_cast<uint64_t>(to_rm) / num_rows + 1;
    switch (ranking_method) {
      case random_kmer: {
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          uint32_t& rix = row_order[ix];
          if (scount_vvec[rix].size() >= n) {
            int64_t r_to_rm;
#pragma omp atomic read
            r_to_rm = to_rm;
            if (r_to_rm > 0) {
              std::vector<unsigned int> ixs;
              getIxsRandom(ixs, enc_vvec[rix].size(), n);
              vecRemoveIxs(scount_vvec[rix], ixs);
              vecRemoveIxs(enc_vvec[rix], ixs);
              vecRemoveIxs(tlca_vvec[rix], ixs);
#pragma omp atomic update
              to_rm = to_rm - ixs.size();
            }
          }
        }
        break;
      }
      case species_count: {
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          uint32_t& rix = row_order[ix];
          if (scount_vvec[rix].size() >= n) {
            int64_t r_to_rm;
#pragma omp atomic read
            r_to_rm = to_rm;
            if (r_to_rm > 0) {
              std::vector<unsigned int> ixs;
              vecIxsNumber(ixs, scount_vvec[rix], n, FL);
              vecRemoveIxs(scount_vvec[rix], ixs);
              vecRemoveIxs(enc_vvec[rix], ixs);
              vecRemoveIxs(tlca_vvec[rix], ixs);
#pragma omp atomic update
              to_rm = to_rm - ixs.size();
            }
          }
        }
        break;
      }
      case information_score: {
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          uint32_t& rix = row_order[ix];
          if (enc_vvec[rix].size() >= n) {
            int64_t r_to_rm;
#pragma omp atomic read
            r_to_rm = to_rm;
            if (r_to_rm > 0) {
              std::vector<unsigned int> ixs;
              std::unordered_map<encT, std::vector<scT>> values_map =
                mapValuesCountsHTd(childrenHT, rix);
              std::vector<scT> scores_vec;
              vecInformationScores(scores_vec, enc_vvec[rix], values_map);
              vecIxsNumber(ixs, scores_vec, n, FL);
              vecRemoveIxs(scount_vvec[rix], ixs);
              vecRemoveIxs(enc_vvec[rix], ixs);
              vecRemoveIxs(tlca_vvec[rix], ixs);
#pragma omp atomic update
              to_rm = to_rm - ixs.size();
            }
          }
        }
        break;
      }
      case taxa_count: {
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          uint32_t& rix = row_order[ix];
          if (enc_vvec[rix].size() >= n) {
            int64_t r_to_rm;
#pragma omp atomic read
            r_to_rm = to_rm;
            if (r_to_rm > 0) {
              std::vector<unsigned int> ixs;
              std::unordered_map<encT, std::vector<bool>> values_map =
                mapValuesBinaryHTd(childrenHT, rix);
              std::vector<scT> tcount_vec;
              vecTaxaCounts(tcount_vec, enc_vvec[rix], values_map);
              vecIxsNumber(ixs, tcount_vec, n, FL);
              vecRemoveIxs(scount_vvec[rix], ixs);
              vecRemoveIxs(enc_vvec[rix], ixs);
              vecRemoveIxs(tlca_vvec[rix], ixs);
#pragma omp atomic update
              to_rm = to_rm - ixs.size();
            }
          }
        }
        break;
      }
    }
#ifdef DEBUG
    LOG(INFO) << "Number of k-mers removed in this round is " << pto_rm - to_rm << std::endl;
#endif
    assert(pto_rm != to_rm);
    pto_rm = to_rm;
  }
  HTd<encT>::updateSize();
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void
HTs<encT>::shrinkHT(uint64_t num_rm)
{
  assert(num_rm <= num_kmers);
  assert(num_rm < std::numeric_limits<int64_t>::max());
  volatile int64_t to_rm = static_cast<int64_t>(num_rm);
  volatile int64_t pto_rm = to_rm;
#ifdef DEBUG
  LOG(INFO) << "Initial number of k-mers to remove is" << to_rm << std::endl;
#endif
  while (to_rm > 0) {
    std::vector<unsigned int> row_order;
    arrSizeOrder(row_order, ind_arr, num_rows, true);
    uint8_t n = to_rm / num_rows + 1;
    switch (ranking_method) {
      case random_kmer: {
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          uint32_t& rix = row_order[ix];
          if (ind_arr[rix] >= n) {
            int64_t r_to_rm;
#pragma omp atomic read
            r_to_rm = to_rm;
            if (r_to_rm > 0) {
              std::vector<unsigned int> ixs;
              getIxsRandom(ixs, ind_arr[rix], n);
              arrRemoveIxs(enc_arr + (rix * b), ind_arr[rix], ixs);
              arrRemoveIxs(scount_arr + (rix * b), ind_arr[rix], ixs);
              arrRemoveIxs(tlca_arr + (rix * b), ind_arr[rix], ixs);
              ind_arr[rix] = ind_arr[rix] - ixs.size();
#pragma omp atomic update
              to_rm = to_rm - ixs.size();
            }
          }
        }
        break;
      }
      case species_count: {
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          uint32_t& rix = row_order[ix];
          if (ind_arr[rix] >= n) {
            int64_t r_to_rm;
#pragma omp atomic read
            r_to_rm = to_rm;
            if (r_to_rm > 0) {
              std::vector<unsigned int> ixs;
              arrIxsNumber(ixs, scount_arr + (rix * b), ind_arr[rix], n, true);
              arrRemoveIxs(enc_arr + (rix * b), ind_arr[rix], ixs);
              arrRemoveIxs(scount_arr + (rix * b), ind_arr[rix], ixs);
              arrRemoveIxs(tlca_arr + (rix * b), ind_arr[rix], ixs);
              ind_arr[rix] = ind_arr[rix] - ixs.size();
#pragma omp atomic update
              to_rm = to_rm - ixs.size();
            }
          }
        }
        break;
      }
      case information_score: {
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          uint32_t& rix = row_order[ix];
          if (ind_arr[rix] >= n) {
            int64_t r_to_rm;
#pragma omp atomic read
            r_to_rm = to_rm;
            if (r_to_rm > 0) {
              std::unordered_map<encT, std::vector<scT>> values_map =
                mapValuesCountsHTs(childrenHT, rix);
              std::vector<scT> scores_vec;
              arrInformationScores(scores_vec, enc_arr + rix * b, ind_arr[rix], values_map);
              std::vector<unsigned int> ixs;
              vecIxsNumber(ixs, scores_vec, n, FL);
              arrRemoveIxs(enc_arr + (rix * b), ind_arr[rix], ixs);
              arrRemoveIxs(scount_arr + (rix * b), ind_arr[rix], ixs);
              arrRemoveIxs(tlca_arr + (rix * b), ind_arr[rix], ixs);
              ind_arr[rix] = ind_arr[rix] - ixs.size();
#pragma omp atomic update
              to_rm = to_rm - ixs.size();
            }
          }
        }
        break;
      }
      case taxa_count: {
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          uint32_t& rix = row_order[ix];
          if (ind_arr[rix] > 0) {
            int64_t r_to_rm;
#pragma omp atomic read
            r_to_rm = to_rm;
            if (r_to_rm > 0) {
              std::unordered_map<encT, std::vector<bool>> values_map =
                mapValuesBinaryHTs(childrenHT, rix);
              std::vector<scT> tcount_vec;
              arrTaxaCounts(tcount_vec, enc_arr + rix * b, ind_arr[rix], values_map);
              std::vector<unsigned int> ixs;
              vecIxsNumber(ixs, tcount_vec, n, FL);
              arrRemoveIxs(enc_arr + (rix * b), ind_arr[rix], ixs);
              arrRemoveIxs(scount_arr + (rix * b), ind_arr[rix], ixs);
              arrRemoveIxs(tlca_arr + (rix * b), ind_arr[rix], ixs);
              ind_arr[rix] = ind_arr[rix] - ixs.size();
#pragma omp atomic update
              to_rm = to_rm - ixs.size();
            }
          }
        }
        break;
      }
    }
#ifdef DEBUG
    LOG(INFO) << "Number of k-mers removed in this round is " << pto_rm - to_rm << std::endl;
#endif
    assert(pto_rm != to_rm);
    pto_rm = to_rm;
  }
  HTs<encT>::updateSize();
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void
HTd<encT>::updateLCA()
{
  if (childrenHT.size() > 1) {
#pragma omp parallel for num_threads(num_threads) schedule(dynamic) private(gen)
    for (uint32_t rix = 0; rix < num_rows; ++rix) {
      tlca_vvec[rix].resize(enc_vvec[rix].size());
      std::unordered_map<encT, std::vector<float>> prob_map =
        mapValuesProbabilityHTd(childrenHT, rix);
      std::unordered_map<encT, std::vector<tT>> tlca_map = mapValuesLCAtHTd(childrenHT, rix);
      for (unsigned int i = 0; i < enc_vvec[rix].size(); ++i) {
        if (tlca_map[enc_vvec[rix][i]].size() == 1) {
          tlca_vvec[rix][i] = tlca_map[enc_vvec[rix][i]][0];
        } else {
          std::unordered_map<encT, std::vector<float>> select_map =
            mapValuesSelectHTd(childrenHT, rix);
          std::bernoulli_distribution d(
            updateLCAtProbabilityAPN(prob_map[enc_vvec[rix][i]], select_map[enc_vvec[rix][i]]));
          /* std::bernoulli_distribution d(updateLCAtProbabilityC2N(prob_map[enc_vvec[rix][i]])); */
          if (d(gen)) {
            tlca_vvec[rix][i] = tID;
          } else {
            std::discrete_distribution<> d(prob_map[enc_vvec[rix][i]].begin(),
                                           prob_map[enc_vvec[rix][i]].begin() +
                                             tlca_map[enc_vvec[rix][i]].size());
            tlca_vvec[rix][i] = tlca_map[enc_vvec[rix][i]][d(gen)];
          }
        }
      }
    }
  }
}

template<typename encT>
void
HTs<encT>::updateLCA()
{
  if (childrenHT.size() > 1) {
#pragma omp parallel for num_threads(num_threads) schedule(dynamic) private(gen)
    for (uint32_t rix = 0; rix < num_rows; ++rix) {
      std::unordered_map<encT, std::vector<float>> prob_map =
        mapValuesProbabilityHTs(childrenHT, rix);
      std::unordered_map<encT, std::vector<tT>> tlca_map = mapValuesLCAtHTs(childrenHT, rix);
      for (unsigned int i = 0; i < ind_arr[rix]; ++i) {
        if (tlca_map[enc_arr[b * rix + i]].size() == 1) {
          tlca_arr[rix * b + i] = tlca_map[enc_arr[rix * b + i]][0];
        } else {
          std::unordered_map<encT, std::vector<float>> select_map =
            mapValuesSelectHTs(childrenHT, rix);
          std::bernoulli_distribution d(updateLCAtProbabilityAPN(prob_map[enc_arr[rix * b + i]],
                                                                 select_map[enc_arr[rix * b + i]]));
          /* std::bernoulli_distribution d(updateLCAtProbabilityC2N(prob_map[enc_arr[rix * b +
           * i]])); */
          if (d(gen)) {
            tlca_arr[rix * b + i] = tID;
          } else {
            std::discrete_distribution<> d(prob_map[enc_arr[rix * b + i]].begin(),
                                           prob_map[enc_arr[rix * b + i]].begin() +
                                             tlca_map[enc_arr[rix * b + i]].size());
            tlca_arr[rix * b + i] = tlca_map[enc_arr[rix * b + i]][d(gen)];
          }
        }
      }
    }
  }
}

template<typename encT>
void
HTd<encT>::convertHTs(HTs<encT>* new_table)
{
  assert(new_table->k == k);
  assert(new_table->h == h);
  assert(new_table->num_rows == num_rows);
  assert(new_table->ptr_lsh_vg == ptr_lsh_vg);
  uint8_t b = new_table->b;
  uint32_t num_rows = new_table->num_rows;
  new_table->clearRows();
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (enc_vvec[rix].size() < b) {
      new_table->ind_arr[rix] = enc_vvec[rix].size();
      if (new_table->ind_arr[rix]) {
        std::copy(enc_vvec[rix].begin(), enc_vvec[rix].end(), new_table->enc_arr + (rix * b));
        std::copy(
          scount_vvec[rix].begin(), scount_vvec[rix].end(), new_table->scount_arr + (rix * b));
        std::copy(tlca_vvec[rix].begin(), tlca_vvec[rix].end(), new_table->tlca_arr + (rix * b));
      }
    } else {
      new_table->ind_arr[rix] = b;
      std::copy(enc_vvec[rix].begin(), enc_vvec[rix].begin() + b, new_table->enc_arr + (rix * b));
      std::copy(
        scount_vvec[rix].begin(), scount_vvec[rix].begin() + b, new_table->scount_arr + (rix * b));
      std::copy(
        tlca_vvec[rix].begin(), tlca_vvec[rix].begin() + b, new_table->tlca_arr + (rix * b));
    }
  }
  new_table->num_species = num_species;
  new_table->tIDsBasis = tIDsBasis;
  new_table->updateSize();
}

#include "tableins.cpp"
