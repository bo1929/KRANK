#include "table.h"

/* #undef DEBUG */
#define FL false

template<typename encT>
std::map<size_t, uint64_t> HTd<encT>::histRowSizes()
{
  std::map<size_t, uint64_t> hist_map;
  for (auto &row : enc_vvec) {
    hist_map[row.size()]++;
  }
  return hist_map;
}

template<typename encT>
std::map<uint8_t, uint64_t> HTs<encT>::histRowSizes()
{
  std::map<uint8_t, uint64_t> hist_map;
  for (unsigned int rix = 0; rix < num_rows; ++rix) {
    hist_map[ind_arr[rix]]++;
  }
  return hist_map;
}

template<typename encT>
inline void
HTd<encT>::mapBefore(std::unordered_map<encT, scT> &scount_map, std::unordered_map<encT, tT> &tlca_map, uint32_t rix)
{
  for (unsigned int i = 0; i < enc_vvec[rix].size(); ++i) {
    scount_map[enc_vvec[rix][i]] = scount_vvec[rix][i];
    /* tlca_map[enc_vvec[rix][i]] = tlca_vvec[rix][i]; */
  }
}

template<typename encT>
inline void
HTd<encT>::updateAfter(std::unordered_map<encT, scT> &scount_map, std::unordered_map<encT, tT> &tlca_map, uint32_t rix)
{
  scount_vvec[rix].resize(enc_vvec[rix].size());
  /* tlca_vvec[rix].resize(enc_vvec[rix].size()); */
  for (unsigned int i = 0; i < enc_vvec[rix].size(); ++i) {
    scount_vvec[rix][i] = scount_map[enc_vvec[rix][i]];
    /* tlca_vvec[rix][i] = tlca_map[enc_vvec[rix][i]]; */
  }
}

template<typename encT>
void HTd<encT>::makeUnique(bool update_size)
{
#ifdef DEBUG
  if (!HTd<encT>::areColumnsSorted()) {
    HTd<encT>::sortColumns();
    std::puts("HTd has unsorted columns before making rows unique.\n");
  }
#endif
#pragma omp parallel for schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!enc_vvec[rix].empty()) {
      std::unordered_map<encT, scT> scount_map{};
      std::unordered_map<encT, tT> tlca_map{};
      mapBefore(scount_map, tlca_map, rix);
      enc_vvec[rix].erase(std::unique(enc_vvec[rix].begin(), enc_vvec[rix].end()), enc_vvec[rix].end());
      updateAfter(scount_map, tlca_map, rix);
    }
  }
  if (update_size)
    HTd<encT>::updateSize();
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
inline void
HTs<encT>::mapBefore(std::unordered_map<encT, scT> &scount_map, std::unordered_map<encT, tT> &tlca_map, uint32_t rix)
{
  uint8_t cix = ind_arr[rix];
  for (unsigned int i = 0; i < cix; ++i) {
    scount_map[enc_arr[rix * b + i]] = scount_arr[rix * b + i];
    /* tlca_map[enc_arr[rix * b + i]] = tlca_arr[rix * b + i]; */
  }
}

template<typename encT>
inline void
HTs<encT>::updateAfter(std::unordered_map<encT, scT> &scount_map, std::unordered_map<encT, tT> &tlca_map, uint32_t rix)
{
  uint8_t cix = ind_arr[rix];
  for (unsigned int i = 0; i < cix; ++i) {
    scount_arr[rix * b + i] = scount_map[enc_arr[rix * b + i]];
    /* tlca_arr[rix * b + i] = tlca_map[enc_arr[rix * b + i]]; */
  }
}

template<typename encT>
void HTs<encT>::makeUnique(bool update_size)
{
#ifdef DEBUG
  if (!HTs<encT>::areColumnsSorted()) {
    HTs<encT>::sortColumns();
    std::puts("HTs has unsorted columns before making rows unique.\n");
  }
#endif
#pragma omp parallel for schedule(static)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (ind_arr[rix] > 0) {
      std::unordered_map<encT, scT> scount_map{};
      std::unordered_map<encT, tT> tlca_map{};
      mapBefore(scount_map, tlca_map, rix);
      auto cix = std::unique(enc_arr + (rix * b), enc_arr + (rix * b) + ind_arr[rix]) - enc_arr - (rix * b);
      ind_arr[rix] = cix;
      updateAfter(scount_map, tlca_map, rix);
    }
  }
  if (update_size)
    HTs<encT>::updateSize();
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void HTd<encT>::sortColumns()
{
#pragma omp parallel for schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!enc_vvec[rix].empty()) {
      std::unordered_map<encT, scT> scount_map{};
      std::unordered_map<encT, tT> tlca_map{};
      mapBefore(scount_map, tlca_map, rix);
      std::sort(enc_vvec[rix].begin(), enc_vvec[rix].end());
      updateAfter(scount_map, tlca_map, rix);
    }
  }
}

template<typename encT>
void HTs<encT>::sortColumns()
{
  uint8_t b = HTs<encT>::b;
#pragma omp parallel for schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (ind_arr[rix] > 0) {
      std::unordered_map<encT, scT> scount_map{};
      std::unordered_map<encT, tT> tlca_map{};
      mapBefore(scount_map, tlca_map, rix);
      std::sort(enc_arr + (rix * b), enc_arr + (rix * b) + ind_arr[rix]);
      updateAfter(scount_map, tlca_map, rix);
    }
  }
}

template<typename encT>
bool HTd<encT>::areColumnsSorted()
{
  bool is_ok = true;
#pragma omp parallel for schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!enc_vvec[rix].empty() && !std::is_sorted(enc_vvec[rix].begin(), enc_vvec[rix].end())) {
#pragma omp atomic write
      is_ok = false;
    }
  }
  return is_ok;
}

template<typename encT>
bool HTs<encT>::areColumnsSorted()
{
  bool is_ok = true;
#pragma omp parallel for schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if ((ind_arr[rix] > 0) && !std::is_sorted(enc_arr + (rix * b), enc_arr + (rix * b) + ind_arr[rix])) {
#pragma omp atomic write
      is_ok = false;
    }
  }
  return is_ok;
}

template<typename encT>
void HTd<encT>::clearRows()
{
#ifdef DEBUG
  LOG(INFO) << "The current size, before clearing, of the table is " << num_kmers << std::endl;
#endif
  enc_vvec.clear();
  enc_vvec.resize(num_rows);
  scount_vvec.clear();
  scount_vvec.resize(num_rows);
  /* tlca_vvec.clear(); */
  /* tlca_vvec.resize(num_rows); */
  HTd<encT>::updateSize();
  HTd<encT>::num_basis = 0;
  HTd<encT>::tIDsBasis = {};
#ifdef DEBUG
  LOG(INFO) << "The current size, after clearing, of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void HTs<encT>::clearRows()
{
#ifdef DEBUG
  LOG(INFO) << "The current size, before clearing, of the table is " << num_kmers << std::endl;
#endif
#pragma omp parallel for schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    ind_arr[rix] = 0;
  }
  HTs<encT>::updateSize();
  num_basis = 0;
  tIDsBasis = {};
#ifdef DEBUG
  LOG(INFO) << "The current size, after clearing, of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void HTd<encT>::updateSize()
{
  uint64_t num_kmers_sum = 0;
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    num_kmers_sum += enc_vvec[rix].size();
  }
  num_kmers = num_kmers_sum;
  table_size = num_kmers * sizeof(encT) + num_rows * sizeof(std::vector<encT>);
  table_size += num_kmers * sizeof(scT) + num_rows * sizeof(std::vector<scT>);
  table_size += num_kmers * sizeof(tT) + num_rows * sizeof(std::vector<tT>);
}

template<typename encT>
void HTs<encT>::updateSize()
{
  uint64_t num_kmers_sum = 0;
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    num_kmers_sum += ind_arr[rix];
  }
  num_kmers = num_kmers_sum;
  table_size = table_size;
}

template<typename encT>
void HTd<encT>::initBasis(tT tID)
{
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
#pragma omp parallel for schedule(static)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!enc_vvec[rix].empty()) {
      scount_vvec[rix].resize(enc_vvec[rix].size());
      std::fill(scount_vvec[rix].begin(), scount_vvec[rix].end(), 1);
      /* tlca_vvec[rix].resize(enc_vvec[rix].size()); */
      /* std::fill(tlca_vvec[rix].begin(), tlca_vvec[rix].end(), tID); */
    } else {
      scount_vvec[rix].clear();
      /* tlca_vvec[rix].clear(); */
    }
  }
  num_basis = 1;
  tIDsBasis = {tID};
}

template<typename encT>
void HTd<encT>::trimColumns(size_t b_max)
{
#pragma omp parallel for schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (enc_vvec[rix].size() > b_max) {
      enc_vvec[rix].resize(b_max);
      scount_vvec[rix].resize(b_max);
      /* tlca_vvec[rix].resize(b_max); */
    }
  }
#ifdef DEBUG
  HTd<encT>::updateSize();
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void HTd<encT>::getScores(std::vector<float> &scores_vec, uint32_t rix)
{
  /* std::random_device rd_; */
  /* std::mt19937 gen_(rd_()); */
  /* std::beta_distribution<> beta_dist(1.0 / childrenHT.size(), 1.0 / childrenHT.size()); */
  /* float tsnorm = 0; */
  /* float rsnorm = 0; */
  /* float nbnorm = 0; */
  /* for (auto &ht : childrenHT) { */
  /*   tsnorm += ht.num_kmers; */
  /*   rsnorm += ht.scount_vvec[rix].size(); */
  /*   nbnorm += ht.num_basis; */
  /* } */
  std::unordered_map<encT, float> values_map{};
  for (auto &ht : childrenHT) {
    /* std::geometric_distribution<> rgdist(static_cast<float>(ht.enc_vvec[rix].size()) / rsnorm); */
    /* std::geometric_distribution<> tgdist(static_cast<float>(ht.num_kmers) / tsnorm); */
    /* std::geometric_distribution<> nbdist(static_cast<float>(ht.num_basis) / nbnorm); */
    for (unsigned int i = 0; i < ht.scount_vvec[rix].size(); ++i) {
      // Approach 1:
      /* std::binomial_distribution<> d( */
      /*   static_cast<unsigned int>(std::ceil(row_weight * (static_cast<float>(ht.num_kmers) / tsnorm))), beta_dist(gen_)); */
      /* values_map[ht.enc_vvec[rix][i]] += (d(gen_) + 1) * ht.scount_vvec[rix][i]; */
      // Approach 2:
      /* values_map[ht.enc_vvec[rix][i]] += log2(nbdist(gen_) + 2) * ht.scount_vvec[rix][i]; */
      // Approach 3:
      values_map[ht.enc_vvec[rix][i]] +=
        static_cast<float>(ht.scount_vvec[rix][i]) / static_cast<float>(ht.scount_vvec[rix].size());
    }
  }
  scores_vec.resize(enc_vvec[rix].size());
  for (unsigned int i = 0; i < enc_vvec[rix].size(); ++i) {
    scores_vec[i] = values_map[enc_vvec[rix][i]];
  }
}

template<typename encT>
void HTd<encT>::selectCoverage(std::vector<size_t> &ixs_r, uint32_t rix, size_t b_max)
{
  std::vector<size_t> ixs, ixs_k, ixs_a;
  ixs.resize(enc_vvec[rix].size());
  assert((ixs.size() > 0) && (b_max > 0));
  std::iota(ixs.begin(), ixs.end(), 0);
  std::random_shuffle(ixs.begin(), ixs.end());
  unsigned int nrm = (ixs.size() - b_max);
  auto &v = scount_vvec[rix];
  if (FL)
    std::sort(ixs.begin(), ixs.end(), [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });
  else
    std::sort(ixs.begin(), ixs.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
  ixs_k.push_back(ixs.back());

  for (unsigned int j = ixs.size() - 1; j-- != 0 && ixs_r.size() < nrm;) {
    bool uc = ixs_k.size() < b_max;
    bool am = false;
    for (unsigned int k = 0; k < ixs_k.size() && uc; ++k) {
      tT c_lca = taxonomy_record->getLowestCommonAncestor(tlca_vvec[rix][ixs[j]], tlca_vvec[rix][ixs_k[k]]);
      if ((c_lca == tlca_vvec[rix][ixs[j]])) {
        // Has lower score and someones parent (or the same taxon).
        uc = false;
      } else if (c_lca == tlca_vvec[rix][ixs_k[k]]) {
        // Ambiguous, k has higher score but up in the tree or they are coming from the same child.
        am = true;
      } else { // ((c_lca != tlca_vvec[rix][ixs_k[k]]) && (c_lca != tlca_vvec[rix][ixs[j]]))
        // We don't know, keep it. They are not coming from the same child.
        uc = true;
      }
    }
    if (uc && !am)
      ixs_k.push_back(ixs[j]);
    else if (am)
      ixs_a.push_back(ixs[j]);
    else
      ixs_r.push_back(ixs[j]);
  }
  for (unsigned int j = ixs_a.size(); j-- != 0 && ixs_r.size() < nrm;) {
    ixs_r.push_back(ixs_a[j]);
  }
  for (unsigned int j = ixs_k.size(); j-- != 0 && ixs_r.size() < nrm;) {
    ixs_r.push_back(ixs_k[j]);
  }
  std::sort(ixs_r.begin(), ixs_r.end());
}

template<typename encT>
void HTd<encT>::pruneColumns(size_t b_max)
{
#ifdef DEBUG
  LOG(INFO) << "The current size, before pruning, of the table is " << num_kmers << std::endl;
#endif
#pragma omp parallel for schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (enc_vvec[rix].size() > b_max) {
      std::vector<size_t> ixs_r;
      switch (ranking_method) {
        case random_kmer: {
          getIxsRandom(ixs_r, enc_vvec[rix].size(), enc_vvec[rix].size() - b_max);
          break;
        }
        case representative_kmer: {
          std::vector<float> scores_vec;
          HTd<encT>::getScores(scores_vec, rix);
          getIxsArgsort(ixs_r, scores_vec, enc_vvec[rix].size() - b_max, FL);
          /* HTd<encT>::selectCoverage(ixs_r, rix, b_max); */
          break;
        }
      }
      vecRemoveIxs(enc_vvec[rix], ixs_r);
      vecRemoveIxs(scount_vvec[rix], ixs_r);
      /* vecRemoveIxs(tlca_vvec[rix], ixs_r); */
    }
  }
  HTd<encT>::updateSize();
#ifdef DEBUG
  if (!HTd<encT>::areColumnsSorted()) {
    HTd<encT>::sortColumns();
    std::puts("HTd has unsorted columns after pruning.\n");
  }
  LOG(INFO) << "The current size, after pruning, of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void HTs<encT>::getScores(std::vector<float> &scores_vec, uint32_t rix)
{
  std::unordered_map<scT, scT> cmap{};
  std::unordered_map<encT, float> values_map{};
  for (auto &ht : childrenHT) {
    for (unsigned int i = 0; i < ht.ind_arr[rix]; ++i)
      cmap[ht.tID] += ht.scount_arr[rix * b + i];
  }
  for (auto &ht : childrenHT) {
    for (unsigned int i = 0; i < ht.ind_arr[rix]; ++i) {
      values_map[ht.enc_arr[rix * b + i]] += static_cast<float>(ht.scount_arr[rix * b + i]) / ht.ind_arr[rix];
    }
  }
  scores_vec.resize(ind_arr[rix]);
  for (unsigned int i = 0; i < ind_arr[rix]; ++i) {
    scores_vec[i] = values_map[enc_arr[rix * b + i]];
  }
}

template<typename encT>
void HTd<encT>::accumulateCounts(std::unordered_map<encT, scT> &scount_map, HTd<encT> &child, uint32_t rix)
{
  for (unsigned int i = 0; i < enc_vvec[rix].size(); ++i)
    scount_map[enc_vvec[rix][i]] = scount_vvec[rix][i];
  std::unordered_map<encT, scT> increment_map{};
  for (unsigned int i = 0; i < child.enc_vvec[rix].size(); ++i)
    increment_map[child.enc_vvec[rix][i]] = child.scount_vvec[rix][i];
  scount_map = std::accumulate(
    increment_map.begin(),
    increment_map.end(),
    scount_map,
    [](std::unordered_map<encT, scT> &p, const std::pair<const encT, scT> &c) { return (p[c.first] += c.second, p); });
}

template<typename encT>
void HTd<encT>::unionRows(HTd<encT> &child, bool update_size)
{
#ifdef DEBUG
  if (!HTd<encT>::areColumnsSorted()) {
    HTd<encT>::sortColumns();
    std::puts("HTd (parent) has unsorted columns before computing the union of rows.\n");
  }
  if (!child.areColumnsSorted()) {
    child.sortColumns();
    std::puts("HTd (child) has unsorted columns before computing the union of rows.\n");
  }
#endif
  assert(num_rows == child.num_rows);
#pragma omp parallel for schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!enc_vvec[rix].empty() && !child.enc_vvec[rix].empty()) {
      std::unordered_map<encT, scT> scount_map{};
      accumulateCounts(scount_map, child, rix);
      enc_vvec[rix].insert(enc_vvec[rix].end(), child.enc_vvec[rix].begin(), child.enc_vvec[rix].end());
      std::inplace_merge(enc_vvec[rix].begin(),
                         enc_vvec[rix].begin() + (enc_vvec[rix].size() - child.enc_vvec[rix].size()),
                         enc_vvec[rix].end());
      enc_vvec[rix].erase(std::unique(enc_vvec[rix].begin(), enc_vvec[rix].end()), enc_vvec[rix].end());
      scount_vvec[rix].resize(enc_vvec[rix].size());
      for (unsigned int i = 0; i < enc_vvec[rix].size(); ++i)
        scount_vvec[rix][i] = scount_map[enc_vvec[rix][i]];
    } else if (!child.enc_vvec[rix].empty()) {
      enc_vvec[rix] = child.enc_vvec[rix];
      scount_vvec[rix] = child.scount_vvec[rix];
      /* tlca_vvec[rix] = child.tlca_vvec[rix]; */
    }
  }
  num_basis += child.num_basis;
  tIDsBasis.insert(child.tIDsBasis.begin(), child.tIDsBasis.end());
  if (update_size)
    HTd<encT>::updateSize();
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void HTs<encT>::accumulateCounts(std::unordered_map<encT, scT> &scount_map, HTs<encT> &child, uint32_t rix)
{
  for (unsigned int i = 0; i < ind_arr[rix]; ++i)
    scount_map[enc_arr[rix * b + i]] = scount_arr[rix * b + i];
  std::unordered_map<encT, scT> increment_map{};
  for (unsigned int i = 0; i < child.ind_arr[rix]; ++i)
    increment_map[child.enc_arr[rix * b + i]] = child.scount_arr[rix * b + i];
  scount_map = std::accumulate(
    increment_map.begin(),
    increment_map.end(),
    scount_map,
    [](std::unordered_map<encT, scT> &p, const std::pair<const encT, scT> &c) { return (p[c.first] += c.second, p); });
}

template<typename encT>
void HTs<encT>::unionRows(HTs<encT> &child, bool update_size)
{
#ifdef DEBUG
  if (!HTs<encT>::areColumnsSorted()) {
    HTs<encT>::sortColumns();
    std::puts("HTs (parent) has unsorted columns before computing the union of rows.\n");
  }
  if (!child.areColumnsSorted()) {
    child.sortColumns();
    std::puts("HTs (child) has unsorted columns before computing the union of rows.\n");
  }
#endif
  assert(b == child.b);
  assert(num_rows == child.num_rows);
#pragma omp parallel for schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    unsigned int ix = rix * b;
    if ((child.ind_arr[rix] > 0) && (ind_arr[rix] > 0)) {
      std::unordered_map<encT, scT> scount_map{};
      accumulateCounts(scount_map, child, rix);
      std::vector<encT> n_enc_arr;
      std::set_union(enc_arr + ix,
                     enc_arr + (ix + ind_arr[rix]),
                     child.enc_arr + ix,
                     child.enc_arr + (ix + child.ind_arr[rix]),
                     std::back_inserter(n_enc_arr));
      if (n_enc_arr.size() > b) {
        std::vector<size_t> ixs_r;
        switch (ranking_method) {
          case random_kmer: {
            getIxsRandom(ixs_r, n_enc_arr.size(), n_enc_arr.size() - b);
            break;
          }
          case representative_kmer: {
            std::vector<float> scores_vec;
            HTs<encT>::getScores(scores_vec, rix);
            getIxsArgsort(ixs_r, scores_vec, n_enc_arr.size() - b, FL);
            break;
          }
        }
        vecRemoveIxs(n_enc_arr, ixs_r);
      }
      std::copy(n_enc_arr.begin(), n_enc_arr.end(), enc_arr + ix);
      ind_arr[rix] = n_enc_arr.size();
      for (unsigned int i = 0; i < ind_arr[rix]; ++i)
        scount_arr[ix + i] = scount_map[enc_arr[ix + i]];
    } else if (child.ind_arr[rix] > 0) {
      std::copy(child.enc_arr + ix, child.enc_arr + (ix + child.ind_arr[rix]), enc_arr + ix);
      std::copy(child.scount_arr + ix, child.scount_arr + (ix + child.ind_arr[rix]), scount_arr + ix);
      /* std::copy(child.tlca_arr + ix, child.tlca_arr + (ix + child.ind_arr[rix]), tlca_arr + ix); */
      ind_arr[rix] = child.ind_arr[rix];
    }
  }
  num_basis += child.num_basis;
  tIDsBasis.insert(child.tIDsBasis.begin(), child.tIDsBasis.end());
  if (update_size)
    HTs<encT>::updateSize();
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

template<typename encT>
void HTd<encT>::getRowOrdering(std::vector<unsigned int> &row_order, bool reverse)
{
  row_order.resize(enc_vvec.size());
  std::iota(row_order.begin(), row_order.end(), 0);
  if (reverse)
    std::sort(row_order.begin(), row_order.end(), [&](unsigned int i1, unsigned int i2) {
      return HTd<encT>::enc_vvec[i1].size() > HTd<encT>::enc_vvec[i2].size();
    });
  else
    std::sort(row_order.begin(), row_order.end(), [&](unsigned int i1, unsigned int i2) {
      return HTd<encT>::enc_vvec[i1].size() < HTd<encT>::enc_vvec[i2].size();
    });
}

template<typename encT>
void HTd<encT>::shrinkHT(uint64_t num_rm, size_t b_max)
{
  assert(num_rm < num_kmers);
  uint64_t init_num_kmers = num_kmers;
  HTd<encT>::pruneColumns(b_max);
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
    HTd<encT>::getRowOrdering(row_order, true);
    size_t n = static_cast<uint64_t>(to_rm) / num_rows + 1;
#pragma omp parallel for schedule(dynamic)
    for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
      uint32_t &rix = row_order[ix];
      if (enc_vvec[rix].size() >= n) {
        int64_t r_to_rm;
#pragma omp atomic read
        r_to_rm = to_rm;
        if (r_to_rm > 0) {
          std::vector<size_t> ixs_r;
          switch (ranking_method) {
            case random_kmer: {
              getIxsRandom(ixs_r, enc_vvec[rix].size(), n);
              break;
            }
            case representative_kmer: {
              std::vector<float> scores_vec;
              HTd<encT>::getScores(scores_vec, rix);
              getIxsArgsort(ixs_r, scores_vec, n, FL);
              /* if (enc_vvec[rix].size() == n) { */
              /*   ixs_r.resize(n); */
              /*   std::iota(ixs_r.begin(), ixs_r.end(), 0); */
              /* } else */
              /*   HTd<encT>::selectCoverage(ixs_r, rix, enc_vvec[rix].size() - n); */
              break;
            }
          }
          vecRemoveIxs(scount_vvec[rix], ixs_r);
          vecRemoveIxs(enc_vvec[rix], ixs_r);
          /* vecRemoveIxs(tlca_vvec[rix], ixs_r); */
#pragma omp atomic update
          to_rm = to_rm - ixs_r.size();
        }
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
void HTs<encT>::getRowOrdering(std::vector<unsigned int> &row_order, bool reverse)
{
  row_order.resize(num_rows);
  std::iota(row_order.begin(), row_order.end(), 0);
  if (reverse)
    std::sort(row_order.begin(), row_order.end(), [&](unsigned int i1, unsigned int i2) {
      return HTs<encT>::ind_arr[i1] > HTs<encT>::ind_arr[i2];
    });
  else
    std::sort(row_order.begin(), row_order.end(), [&](unsigned int i1, unsigned int i2) {
      return HTs<encT>::ind_arr[i1] < HTs<encT>::ind_arr[i2];
    });
}

template<typename encT>
void HTs<encT>::shrinkHT(uint64_t num_rm)
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
    HTs<encT>::getRowOrdering(row_order, true);
    size_t n = to_rm / num_rows + 1;
#pragma omp parallel for schedule(dynamic)
    for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
      uint32_t &rix = row_order[ix];
      if (ind_arr[rix] >= n) {
        int64_t r_to_rm;
#pragma omp atomic read
        r_to_rm = to_rm;
        if (r_to_rm > 0) {
          std::vector<size_t> ixs_r;
          switch (ranking_method) {
            case random_kmer: {
              getIxsRandom(ixs_r, ind_arr[rix], n);
              break;
            }
            case representative_kmer: {
              std::vector<float> scores_vec;
              HTs<encT>::getScores(scores_vec, rix);
              getIxsArgsort(ixs_r, scores_vec, n, FL);
              break;
            }
          }
          arrRemoveIxs(enc_arr + (rix * b), ind_arr[rix], ixs_r);
          arrRemoveIxs(scount_arr + (rix * b), ind_arr[rix], ixs_r);
          /* arrRemoveIxs(tlca_arr + (rix * b), ind_arr[rix], ixs_r); */
          ind_arr[rix] = ind_arr[rix] - ixs_r.size();
#pragma omp atomic update
          to_rm = to_rm - ixs_r.size();
        }
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
void HTd<encT>::updateLCA()
{
  if (childrenHT.size() > 1) {
#pragma omp parallel for schedule(dynamic)
    for (uint32_t rix = 0; rix < num_rows; ++rix) {
      if (enc_vvec[rix].size() > 0) {
        std::unordered_map<encT, tT> tlca_map{};
        for (auto &ht : childrenHT) {
          for (unsigned int i = 0; i < ht.enc_vvec[rix].size(); ++i) {
            if (tlca_map.find(ht.enc_vvec[rix][i]) == tlca_map.end())
              tlca_map[ht.enc_vvec[rix][i]] = ht.tlca_vvec[rix][i];
            else
              tlca_map[ht.enc_vvec[rix][i]] = tID;
          }
        }
        tlca_vvec[rix].resize(enc_vvec[rix].size());
        for (unsigned int i = 0; i < enc_vvec[rix].size(); ++i) {
          tlca_vvec[rix][i] = tlca_map[enc_vvec[rix][i]];
        }
      }
    }
  }
}

template<typename encT>
void HTs<encT>::updateLCA()
{
  if (childrenHT.size() > 1) {
#pragma omp parallel for schedule(dynamic)
    for (uint32_t rix = 0; rix < num_rows; ++rix) {
      if (ind_arr[rix] > 0) {
        std::unordered_map<encT, tT> tlca_map{};
        for (auto &ht : childrenHT) {
          for (unsigned int i = 0; i < ht.ind_arr[rix]; ++i) {
            if (tlca_map.find(ht.enc_arr[rix * ht.b + i]) == tlca_map.end())
              tlca_map[ht.enc_arr[rix * ht.b + i]] = ht.tlca_arr[rix * ht.b + i];
            else
              tlca_map[ht.enc_arr[rix * ht.b + i]] = tID;
          }
        }
        for (unsigned int i = 0; i < ind_arr[rix]; ++i) {
          tlca_arr[rix * b + i] = tlca_map[enc_arr[rix * b + i]];
        }
      }
    }
  }
}

template<typename encT>
void HTd<encT>::convertHTs(HTs<encT> *new_table)
{
  assert(new_table->k == k);
  assert(new_table->h == h);
  assert(new_table->num_rows == num_rows);
  assert(new_table->ptr_lsh_vg == ptr_lsh_vg);
  size_t b = new_table->b;
  uint32_t num_rows = new_table->num_rows;
  /* new_table->clearRows(); */
  HTd<encT>::pruneColumns(b);
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (enc_vvec[rix].size() < b) {
      new_table->ind_arr[rix] = enc_vvec[rix].size();
      if (new_table->ind_arr[rix] > 0) {
        std::copy(enc_vvec[rix].begin(), enc_vvec[rix].end(), new_table->enc_arr + (rix * b));
        std::copy(scount_vvec[rix].begin(), scount_vvec[rix].end(), new_table->scount_arr + (rix * b));
        /* std::copy(tlca_vvec[rix].begin(), tlca_vvec[rix].end(), new_table->tlca_arr + (rix * b)); */
      }
    } else {
      new_table->ind_arr[rix] = b;
      std::copy(enc_vvec[rix].begin(), enc_vvec[rix].begin() + b, new_table->enc_arr + (rix * b));
      std::copy(scount_vvec[rix].begin(), scount_vvec[rix].begin() + b, new_table->scount_arr + (rix * b));
      /* std::copy(tlca_vvec[rix].begin(), tlca_vvec[rix].begin() + b, new_table->tlca_arr + (rix * b)); */
    }
  }
  new_table->num_basis = num_basis;
  new_table->tIDsBasis = tIDsBasis;
  new_table->updateSize();
}

template<typename encT>
void HTd<encT>::filterLSR(std::vector<uint8_t> &depth_vec, uint8_t slr_depth)
{
#pragma omp parallel for schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (tlca_vvec[rix].size() >= 1) {
      std::vector<size_t> ixs_r;
      for (unsigned int j = 0; j < tlca_vvec[rix].size(); ++j) {
        if (depth_vec[tlca_vvec[rix][j]] <= slr_depth)
          ixs_r.push_back(j);
      }
      vecRemoveIxs(scount_vvec[rix], ixs_r);
      vecRemoveIxs(enc_vvec[rix], ixs_r);
      vecRemoveIxs(tlca_vvec[rix], ixs_r);
    }
  }
  HTd<encT>::updateSize();
#ifdef DEBUG
  LOG(INFO) << "The current size of the table is " << num_kmers << std::endl;
#endif
}

#include "tableins.cpp"
