#include "table.h"

#define SIZE_ESTIMATE 4194304

// TODO: Benchmark different thread schedules etc.
// TODO: Check what is sorted and what is not.
// TODO: Check if size is updated or not.
// TODO: Add macros for debugging to avoid redundant computation.

template<typename T>
inline void
sortColumns(vvec<T>& table)
{
  uint32_t num_rows = table.size();
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!table[rix].empty())
      std::sort(table[rix].begin(), table[rix].end());
  }
}

template<typename encT>
uint64_t
StreamIM<encT>::readBatch(const char* filepath, unsigned int batch_size)
{
  kseq_t* reader = IO::getReader(filepath);
  batch_size = IO::adjustBatchSize(batch_size, num_threads);
  std::vector<sseq_t> seqBatch = IO::readBatch(reader, batch_size);
  uint64_t total_kmer = 0;
  uint32_t max_rix = pow(2, 2 * StreamIM<encT>::h);
  StreamIM<encT>::lsh_enc_vec.reserve(SIZE_ESTIMATE);
  while (!(seqBatch.empty())) {
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
    for (uint32_t ix = 0; ix < seqBatch.size(); ++ix) {
      const char* kmer_seq;
      // TODO: Consider 32bit encodings if possible.
      uint64_t enc_bp;
      uint64_t enc_lr;
      uint32_t rix;
      kmer_seq = seqBatch[ix].nseq.c_str();
      kmerEncodingCompute(kmer_seq, enc_lr, enc_bp);
      rix = computeValueLSH(enc_bp, *(StreamIM<encT>::ptr_lsh_vg));
      assert(rix < max_rix);
#pragma omp critical
      {
        StreamIM<encT>::lsh_enc_vec.push_back(std::make_pair(rix, enc_lr));
      }
    }
    total_kmer += seqBatch.size();
    if (seqBatch.size() == batch_size)
      seqBatch = IO::readBatch(reader, batch_size);
    else
      seqBatch.clear();
  }
  StreamIM<encT>::lsh_enc_vec.shrink_to_fit();
  std::sort(StreamIM<encT>::lsh_enc_vec.begin(),
            StreamIM<encT>::lsh_enc_vec.end(),
            [](const std::pair<uint32_t, uint64_t>& l, const std::pair<uint32_t, uint64_t>& r) { return l.first < r.first; });
  return total_kmer;
}

template<typename encT>
std::unordered_map<uint8_t, uint64_t>
StreamIM<encT>::histRowSizes()
{
  std::unordered_map<uint32_t, uint8_t> row_sizes;
  for (auto kv : StreamIM<encT>::lsh_enc_vec) {
    row_sizes[kv.first]++;
  }
  std::unordered_map<uint8_t, uint64_t> hist_map;
  for (auto kv : row_sizes) {
    hist_map[kv.second]++;
  }
  return hist_map;
}

template<typename encT>
std::unordered_map<uint8_t, uint64_t>
HTd<encT>::histRowSizes()
{
  std::unordered_map<uint8_t, uint64_t> hist_map;
  for (auto& row : HTd<encT>::enc_vvec) {
    hist_map[row.size()]++;
  }
  return hist_map;
}

template<typename encT>
std::unordered_map<uint8_t, uint64_t>
HTs<encT>::histRowSizes()
{
  std::unordered_map<uint8_t, uint64_t> hist_map;
  uint32_t num_rows = HTs<encT>::num_rows;
  for (unsigned int rix = 0; rix < num_rows; ++rix) {
    hist_map[HTs<encT>::ind_arr[rix]]++;
  }
  return hist_map;
}

template<typename encT>
bool
StreamIM<encT>::save(const char* filepath)
{
  bool is_ok = false;
  if (StreamIM<encT>::lsh_enc_vec.empty()) {
    std::puts("The LSH-value and encoding pair vector is empty, nothing to save!");
    return is_ok;
  }
  FILE* vec_f = IO::open_file(filepath, is_ok, "wb");
  std::fwrite(StreamIM<encT>::lsh_enc_vec.data(), StreamIM<encT>::lsh_enc_vec.size(), sizeof(std::pair<uint32_t, encT>), vec_f);
  if (std::ferror(vec_f)) {
    std::puts("I/O error when writing LSH-value and encoding pairs.");
  } else if (std::feof(vec_f)) {
    is_ok = true;
  }
  std::fclose(vec_f);
  return is_ok;
}

template<typename encT>
bool
StreamIM<encT>::load(const char* filepath)
{
  bool is_ok = false;
  StreamIM<encT>::lsh_enc_vec.clear();
  std::ifstream vec_ifs = IO::open_ifstream(filepath, is_ok);
  while (!vec_ifs.eof() && vec_ifs.good()) {
    std::pair<uint32_t, encT> lsh_enc;
    vec_ifs.read((char*)&lsh_enc, sizeof(std::pair<uint32_t, encT>));
    if (!vec_ifs.eof())
      StreamIM<encT>::lsh_enc_vec.push_back(lsh_enc);
  }
  if (vec_ifs.fail() && !vec_ifs.eof()) {
    std::puts("I/Oerror when reading LSH-value and encoding pairs.");
  } else if (vec_ifs.eof()) {
    is_ok = true;
  }
  vec_ifs.close();
  return is_ok;
}

template<typename encT>
uint64_t
StreamIM<encT>::getBatch(vvec<encT>& batch_table, uint32_t batch_size)
{
  assert(batch_table.size() >= batch_size);
  uint64_t num_kmers = 0;
  while (StreamIM<encT>::curr_vix < StreamIM<encT>::lsh_enc_vec.size()) {
    std::pair<uint32_t, encT> lsh_enc = StreamIM<encT>::lsh_enc_vec[StreamIM<encT>::curr_vix];
    if (lsh_enc.first < (StreamIM<encT>::l_rix + batch_size)) {
      batch_table[lsh_enc.first - StreamIM<encT>::l_rix].push_back(lsh_enc.second);
      num_kmers++;
    } else {
      break;
    }
    StreamIM<encT>::curr_vix++;
  }
  StreamIM<encT>::l_rix += batch_size;
  sortColumns(batch_table);
  return num_kmers;
}

template<typename encT>
bool
StreamOD<encT>::openStream()
{
  bool is_ok = true;
  StreamOD<encT>::vec_ifs = IO::open_ifstream(StreamOD::filepath, is_ok);
  return is_ok;
}

template<typename encT>
uint64_t
StreamOD<encT>::getBatch(vvec<encT>& batch_table, uint32_t batch_size)
{
  assert(batch_table.size() >= batch_size);
  uint64_t num_kmers = 0;
  uint64_t row_ix = 0;
  uint32_t curr_rix = StreamOD<encT>::curr_rix;
  uint32_t f_rix = StreamOD<encT>::f_rix;
  while ((row_ix < batch_size) && (!StreamOD::vec_ifs.eof() && StreamOD::vec_ifs.good())) {
    std::pair<uint32_t, encT> lsh_enc;
    StreamOD::vec_ifs.read((char*)&lsh_enc, sizeof(std::pair<uint32_t, encT>));
    curr_rix = lsh_enc.first;
    row_ix = curr_rix - f_rix;
    if (!StreamOD::vec_ifs.eof() && (row_ix < batch_size)) {
      assert(row_ix < batch_table.size());
      batch_table[row_ix].push_back(lsh_enc.second);
      StreamOD<encT>::curr_rix = curr_rix + 1;
      StreamOD<encT>::curr_pos = StreamOD<encT>::vec_ifs.tellg();
      num_kmers++;
    }
  }
  if (StreamOD<encT>::vec_ifs.fail() && !StreamOD<encT>::vec_ifs.eof()) {
    std::puts("I/Oerror when reading LSH-value and encoding pairs.");
  } else if (StreamOD<encT>::vec_ifs.eof()) {
    StreamOD<encT>::vec_ifs.close();
  } else {
    StreamOD<encT>::vec_ifs.seekg(curr_pos);
    StreamOD<encT>::f_rix = f_rix + batch_size;
  }
  sortColumns(batch_table);
  return num_kmers;
}

template<typename encT>
void
HTd<encT>::makeUnique(bool update_size)
{
#ifdef DEBUG
  if (!HTd<encT>::areColumnsSorted()) {
    HTd<encT>::sortColumns();
    std::puts("HTd has unsorted columns before making columns unique.");
  }
#endif
  uint32_t num_rows = HTd<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!HTd<encT>::enc_vvec[rix].empty()) {
      std::unordered_map<encT, scT> count_map = mapCountsHTd(HTd<encT>::enc_vvec[rix], HTd<encT>::scount_vvec[rix]);
      HTd<encT>::enc_vvec[rix].erase(std::unique(HTd<encT>::enc_vvec[rix].begin(), HTd<encT>::enc_vvec[rix].end()), HTd<encT>::enc_vvec[rix].end());
      updateCountsHTd(HTd<encT>::enc_vvec[rix], HTd<encT>::scount_vvec[rix], count_map);
    }
  }
  if (update_size)
    HTd<encT>::updateSize();
}

template<typename encT>
void
HTs<encT>::makeUnique(bool update_size)
{
#ifdef DEBUG
  if (!HTs<encT>::areColumnsSorted()) {
    HTs<encT>::sortColumns();
    std::puts("HTs has unsorted columns before making columns unique.");
  }
#endif
  uint8_t b = HTs<encT>::b;
  uint32_t num_rows = HTs<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (HTs<encT>::ind_arr[rix] > 0) {
      std::unordered_map<encT, scT> count_map = mapCountsHTs(HTs<encT>::enc_arr + (rix * b), HTs<encT>::scount_arr + (rix * b), HTs<encT>::ind_arr[rix]);
      auto last = std::unique(HTs<encT>::enc_arr + (rix * b), HTs<encT>::enc_arr + (rix * b) + HTs<encT>::ind_arr[rix]) - HTs<encT>::enc_arr - (rix * b);
      HTs<encT>::ind_arr[rix] = last;
      updateCountsHTs(HTs<encT>::enc_arr + (rix * b), HTs<encT>::scount_arr + (rix * b), HTs<encT>::ind_arr[rix], count_map);
    }
  }
  if (update_size)
    HTs<encT>::updateSize();
}

template<typename encT>
void
HTd<encT>::sortColumns()
{
  uint32_t num_rows = HTd<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!HTd<encT>::enc_vvec[rix].empty()) {
      std::unordered_map<encT, scT> count_map = mapCountsHTd(HTd<encT>::enc_vvec[rix], HTd<encT>::scount_vvec[rix]);
      std::sort(HTd<encT>::enc_vvec[rix].begin(), HTd<encT>::enc_vvec[rix].end());
      updateCountsHTd(HTd<encT>::enc_vvec[rix], HTd<encT>::scount_vvec[rix], count_map);
    }
  }
}

template<typename encT>
void
HTs<encT>::sortColumns()
{
  uint8_t b = HTs<encT>::b;
  uint32_t num_rows = HTs<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (HTs<encT>::ind_arr[rix] > 0) {
      std::unordered_map<encT, scT> count_map = mapCountsHTs(HTs<encT>::enc_arr + (rix * b), HTs<encT>::scount_arr + (rix * b), HTs<encT>::ind_arr[rix]);
      std::sort(HTs<encT>::enc_arr + (rix * b), HTs<encT>::enc_arr + (rix * b) + HTs<encT>::ind_arr[rix]);
      updateCountsHTs(HTs<encT>::enc_arr + (rix * b), HTs<encT>::scount_arr + (rix * b), HTs<encT>::ind_arr[rix], count_map);
    }
  }
}

template<typename encT>
bool
HTd<encT>::areColumnsSorted()
{
  bool is_ok = true;
  uint32_t num_rows = HTd<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!HTd<encT>::enc_vvec[rix].empty() && !std::is_sorted(HTd<encT>::enc_vvec[rix].begin(), HTd<encT>::enc_vvec[rix].end())) {
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
  uint8_t b = HTs<encT>::b;
  uint32_t num_rows = HTs<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if ((HTs<encT>::ind_arr[rix] > 0) && !std::is_sorted(HTs<encT>::enc_arr + (rix * b), HTs<encT>::enc_arr + (rix * b) + HTs<encT>::ind_arr[rix])) {
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
  HTd<encT>::enc_vvec.clear();
  HTd<encT>::enc_vvec.resize(HTd<encT>::num_rows);
  HTd<encT>::updateSize();
  HTd<encT>::scount_vvec.clear();
  HTd<encT>::scount_vvec.resize(HTd<encT>::num_rows);
  HTd<encT>::num_species = 0;
}

template<typename encT>
void
HTs<encT>::clearRows()
{
  uint32_t num_rows = HTs<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    HTs<encT>::ind_arr[rix] = 0;
  }
  HTs<encT>::num_species = 0;
}

template<typename encT>
void
HTd<encT>::updateSize()
{
  uint64_t num_kmers = 0;
  for (auto& row : enc_vvec) {
    if (!row.empty())
      num_kmers += row.size();
  }
  HTd<encT>::num_kmers = num_kmers;
  HTd<encT>::table_size = HTd<encT>::num_kmers * sizeof(encT) + HTd<encT>::num_rows * sizeof(std::vector<encT>);
  HTd<encT>::table_size += HTd<encT>::num_kmers * sizeof(scT) + HTd<encT>::num_rows * sizeof(std::vector<scT>);
}

template<typename encT>
void
HTs<encT>::updateSize()
{
  uint64_t num_kmers = 0;
  uint32_t num_rows = HTs<encT>::num_rows;
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    num_kmers += HTs<encT>::ind_arr[rix];
  }
  HTs<encT>::num_kmers = num_kmers;
  HTs<encT>::table_size = HTs<encT>::table_size;
}

template<typename encT>
void
HTd<encT>::initCounts()
{
  uint32_t num_rows = HTd<encT>::num_rows;
  uint32_t value;
  if (HTd<encT>::num_kmers > 0)
    value = 1;
  else
    value = 0;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!HTd<encT>::enc_vvec[rix].empty()) {
      HTd<encT>::scount_vvec[rix].resize(HTd<encT>::enc_vvec[rix].size());
      std::fill(HTd<encT>::scount_vvec[rix].begin(), HTd<encT>::scount_vvec[rix].end(), value);
    } else {
      HTd<encT>::scount_vvec[rix].clear();
    }
  }
  HTd<encT>::num_species = value;
}

template<typename encT>
void
HTd<encT>::trimColumns(uint8_t b_max)
{
  uint32_t num_rows = HTd<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (HTd<encT>::enc_vvec[rix].size() > b_max) {
      HTd<encT>::enc_vvec[rix].resize(b_max);
      HTd<encT>::scount_vvec[rix].resize(b_max);
    }
  }
  HTd<encT>::updateSize();
}

template<typename encT>
void
HTd<encT>::pruneColumns(uint8_t b_max)
{
  assert(b_max > 0);
  uint32_t num_rows = HTd<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (HTd<encT>::enc_vvec[rix].size() > b_max) {
      std::vector<unsigned int> ixs;
      switch (HTd<encT>::rm_strategy) {
        case random_kmer:
          getIxsRandom(ixs, HTd<encT>::enc_vvec[rix].size(), HTd<encT>::enc_vvec[rix].size() - b_max);
          break;
        case small_scount:
          vecIxsNumber(ixs, HTd<encT>::scount_vvec[rix], HTd<encT>::enc_vvec[rix].size() - b_max);
          break;
      }
      vecRemoveIxs(HTd<encT>::enc_vvec[rix], ixs);
      vecRemoveIxs(HTd<encT>::scount_vvec[rix], ixs);
    }
  }
  HTd<encT>::updateSize();
#ifdef DEBUG
  if (!HTd<encT>::areColumnsSorted()) {
    HTd<encT>::sortColumns();
    std::puts("HTd has unsorted columns after pruning.");
  }
#endif
}

template<typename encT>
void
HTd<encT>::unionRows(HTd<encT>& sibling, bool update_size)
{
#ifdef DEBUG
  if (!HTd<encT>::areColumnsSorted()) {
    HTd<encT>::sortColumns();
    std::puts("HTd (main) has unsorted columns before computing the union of rows.");
  }
  if (!sibling.areColumnsSorted()) {
    sibling.sortColumns();
    std::puts("HTd (sibling) has unsorted columns before computing the union of rows.");
  }
#endif
  assert(HTd<encT>::num_rows == sibling.num_rows);
  uint32_t num_rows = HTd<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!HTd<encT>::enc_vvec[rix].empty() && !sibling.enc_vvec[rix].empty()) {
      std::unordered_map<encT, scT> count_map =
        mapTotalCountsHTd(HTd<encT>::enc_vvec[rix], sibling.enc_vvec[rix], HTd<encT>::scount_vvec[rix], sibling.scount_vvec[rix]);
      HTd<encT>::enc_vvec[rix].insert(HTd<encT>::enc_vvec[rix].end(), sibling.enc_vvec[rix].begin(), sibling.enc_vvec[rix].end());
      std::inplace_merge(HTd<encT>::enc_vvec[rix].begin(),
                         HTd<encT>::enc_vvec[rix].begin() + (HTd<encT>::enc_vvec[rix].size() - sibling.enc_vvec[rix].size()),
                         HTd<encT>::enc_vvec[rix].end());
      HTd<encT>::enc_vvec[rix].erase(std::unique(HTd<encT>::enc_vvec[rix].begin(), HTd<encT>::enc_vvec[rix].end()), HTd<encT>::enc_vvec[rix].end());
      updateCountsHTd(HTd<encT>::enc_vvec[rix], HTd<encT>::scount_vvec[rix], count_map);
    } else if (!sibling.enc_vvec[rix].empty()) {
      HTd<encT>::enc_vvec[rix] = sibling.enc_vvec[rix];
      HTd<encT>::scount_vvec[rix] = sibling.scount_vvec[rix];
    }
  }
  HTd<encT>::num_species += sibling.num_species;
  if (update_size)
    HTd<encT>::updateSize();
}

template<typename encT>
void
HTs<encT>::unionRows(HTs<encT>& sibling, bool update_size)
{
#ifdef DEBUG
  if (!HTs<encT>::areColumnsSorted()) {
    HTs<encT>::sortColumns();
    std::puts("HTs (main) has unsorted columns before computing the union of rows.");
  }
  if (!sibling.areColumnsSorted()) {
    sibling.sortColumns();
    std::puts("HTs (sibling) has unsorted columns before computing the union of rows.");
  }
#endif
  assert(HTs<encT>::b == sibling.b);
  assert(HTs<encT>::num_rows * HTs<encT>::b == sibling.num_rows * sibling.b);
  uint8_t b = HTs<encT>::b;
  uint32_t num_rows = HTs<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    unsigned int ix = rix * b;
    if ((sibling.ind_arr[rix] > 0) && (HTs<encT>::ind_arr[rix] > 0)) {
      std::unordered_map<encT, scT> count_map = mapTotalCountsHTs(
        HTs<encT>::enc_arr + ix, HTs<encT>::scount_arr + ix, HTs<encT>::ind_arr[rix], sibling.enc_arr + ix, sibling.scount_arr + ix, sibling.ind_arr[rix]);
      std::vector<encT> v;
      std::set_union(HTs<encT>::enc_arr + ix,
                     HTs<encT>::enc_arr + (ix + HTs<encT>::ind_arr[rix]),
                     sibling.enc_arr + ix,
                     sibling.enc_arr + (ix + sibling.ind_arr[rix]),
                     std::back_inserter(v));
      if (v.size() > b) {
        std::vector<unsigned int> ixs;
        switch (HTs<encT>::rm_strategy) {
          case random_kmer:
            getIxsRandom(ixs, v.size(), v.size() - b);
            break;
          case small_scount:
            std::vector<scT> s_v;
            s_v.resize(v.size());
            std::transform(v.begin(), v.end(), back_inserter(s_v), [&count_map](encT val) { return count_map[val]; });
            vecIxsNumber(ixs, s_v, v.size() - b);
            break;
        }
        vecRemoveIxs(v, ixs);
      }
      std::copy(v.begin(), v.end(), HTs<encT>::enc_arr + ix);
      HTs<encT>::ind_arr[rix] = v.size();
      updateCountsHTs(HTs<encT>::enc_arr + ix, HTs<encT>::scount_arr + ix, HTs<encT>::ind_arr[rix], count_map);
    } else if (sibling.ind_arr[rix] > 0) {
      std::copy(sibling.enc_arr + ix, sibling.enc_arr + (ix + sibling.ind_arr[rix]), HTs<encT>::enc_arr + ix);
      std::copy(sibling.scount_arr + ix, sibling.scount_arr + (ix + sibling.ind_arr[rix]), HTs<encT>::scount_arr + ix);
      HTs<encT>::ind_arr[rix] = sibling.ind_arr[rix];
    }
  }
  HTs<encT>::num_species += sibling.num_species;
  if (update_size)
    HTs<encT>::updateSize();
}

template<typename encT>
void
HTd<encT>::mergeRows(HTd<encT>& sibling, bool update_size)
{
#ifdef DEBUG
  if (!HTd<encT>::areColumnsSorted()) {
    HTd<encT>::sortColumns();
    std::puts("HTd (main) has unsorted columns before merging the rows.");
  }
  if (!sibling.areColumnsSorted()) {
    sibling.sortColumns();
    std::puts("HTd (sibling) has unsorted columns before merging the rows.");
  }
#endif
  assert(HTd<encT>::num_rows == sibling.num_rows);
  uint32_t num_rows = HTd<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!HTd<encT>::enc_vvec[rix].empty() && !sibling.enc_vvec[rix].empty()) {
      std::unordered_map<encT, scT> count_map =
        mapTotalCountsHTd(HTd<encT>::enc_vvec[rix], sibling.enc_vvec[rix], HTd<encT>::scount_vvec[rix], sibling.scount_vvec[rix]);
      HTd<encT>::enc_vvec[rix].insert(HTd<encT>::enc_vvec[rix].end(), sibling.enc_vvec[rix].begin(), sibling.enc_vvec[rix].end());
      std::inplace_merge(HTd<encT>::enc_vvec[rix].begin(),
                         HTd<encT>::enc_vvec[rix].begin() + (HTd<encT>::enc_vvec[rix].size() - sibling.enc_vvec[rix].size()),
                         HTd<encT>::enc_vvec[rix].end());
      updateCountsHTd(HTd<encT>::enc_vvec[rix], HTd<encT>::scount_vvec[rix], count_map);
    } else if (!sibling.enc_vvec[rix].empty()) {
      HTd<encT>::enc_vvec[rix] = sibling.enc_vvec[rix];
      HTd<encT>::scount_vvec[rix] = sibling.scount_vvec[rix];
    }
  }
  HTd<encT>::num_species += sibling.num_species;
  if (update_size)
    HTd<encT>::updateSize();
}

template<typename encT>
void
HTs<encT>::mergeRows(HTs<encT>& sibling, bool update_size)
{
#ifdef DEBUG
  if (!HTs<encT>::areColumnsSorted()) {
    HTs<encT>::sortColumns();
    std::puts("HTs (main) has unsorted columns before merging the rows.");
  }
  if (!sibling.areColumnsSorted()) {
    sibling.sortColumns();
    std::puts("HTs (sibling) has unsorted columns before merging the rows.");
  }
#endif
  assert(HTs<encT>::b == sibling.b);
  assert(HTs<encT>::num_rows * HTs<encT>::b == sibling.num_rows * sibling.b);
  uint8_t b = HTs<encT>::b;
  uint32_t num_rows = HTs<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    unsigned int ix = rix * b;
    if ((sibling.ind_arr[rix] > 0) && (HTs<encT>::ind_arr[rix] > 0)) {
      std::unordered_map<encT, scT> count_map = mapTotalCountsHTs(
        HTs<encT>::enc_arr + ix, HTs<encT>::scount_arr + ix, HTs<encT>::ind_arr[rix], sibling.enc_arr + ix, sibling.scount_arr + ix, sibling.ind_arr[rix]);
      std::vector<encT> v;
      std::merge(HTs<encT>::enc_arr + ix,
                 HTs<encT>::enc_arr + (ix + HTs<encT>::ind_arr[rix]),
                 sibling.enc_arr + ix,
                 sibling.enc_arr + (ix + sibling.ind_arr[rix]),
                 std::back_inserter(v));
      if (v.size() > b) {
        v.erase(std::unique(v.begin(), v.end()), v.end());
      }
      if (v.size() > b) {
        std::vector<unsigned int> ixs;
        switch (HTs<encT>::rm_strategy) {
          case random_kmer:
            getIxsRandom(ixs, v.size(), v.size() - b);
            break;
          case small_scount:
            std::vector<scT> s_v;
            s_v.resize(v.size());
            std::transform(v.begin(), v.end(), back_inserter(s_v), [&count_map](encT val) { return count_map[val]; });
            vecIxsNumber(ixs, s_v, v.size() - b);
            break;
        }
        vecRemoveIxs(v, ixs);
      }
      std::copy(v.begin(), v.end(), HTs<encT>::enc_arr + ix);
      HTs<encT>::ind_arr[rix] = v.size();
      updateCountsHTs(HTs<encT>::enc_arr + ix, HTs<encT>::scount_arr + ix, HTs<encT>::ind_arr[rix], count_map);
    } else if (sibling.ind_arr[rix] > 0) {
      std::copy(sibling.enc_arr + ix, sibling.enc_arr + (ix + sibling.ind_arr[rix]), HTs<encT>::enc_arr + ix);
      std::copy(sibling.scount_arr + ix, sibling.scount_arr + (ix + sibling.ind_arr[rix]), HTs<encT>::scount_arr + ix);
      HTs<encT>::ind_arr[rix] = sibling.ind_arr[rix];
    }
  }
  HTs<encT>::num_species += sibling.num_species;
  if (update_size)
    HTs<encT>::updateSize();
}

template<typename encT>
void
HTd<encT>::shrinkHT(uint64_t num_rm, uint8_t b_max)
{
  uint64_t init_num_kmers = HTd<encT>::num_kmers;
  assert(num_rm < HTd<encT>::num_kmers);
  HTd<encT>::pruneColumns(b_max);
  num_rm = num_rm - (init_num_kmers - HTd<encT>::num_kmers);
  assert(num_rm < std::numeric_limits<int64_t>::max());
  volatile int64_t to_rm = static_cast<uint64_t>(num_rm);
  while (to_rm > 0) {
    std::vector<unsigned int> row_order;
    vvecSizeOrder(row_order, HTd<encT>::scount_vvec, true);
    switch (HTd<encT>::rm_strategy) {
      case random_kmer:
#pragma omp parallel for num_threads(num_threads) schedule(static, 1) shared(to_rm)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          int64_t r_to_rm;
#pragma omp atomic read
          r_to_rm = to_rm;
          if (r_to_rm > 0) {
            std::vector<unsigned int> ixs;
            uint8_t n = r_to_rm / HTd<encT>::num_rows + 1;
            getIxsRandom(ixs, HTd<encT>::scount_vvec[row_order[ix]].size(), n);
            vecRemoveIxs(HTd<encT>::scount_vvec[row_order[ix]], ixs);
            vecRemoveIxs(HTd<encT>::enc_vvec[row_order[ix]], ixs);
#pragma omp atomic update
            to_rm = to_rm - ixs.size();
          } else {
            continue;
          }
        }
        break;
      case small_scount:
        scT threshold = vvecArgmax2D(HTd<encT>::scount_vvec, static_cast<uint64_t>(to_rm));
#pragma omp parallel for num_threads(num_threads) schedule(static, 1) shared(to_rm)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          int64_t r_to_rm;
#pragma omp atomic read
          r_to_rm = to_rm;
          if (r_to_rm > 0) {
            std::vector<unsigned int> ixs;
            vecIxsThreshold(ixs, HTd<encT>::scount_vvec[row_order[ix]], threshold);
            vecRemoveIxs(HTd<encT>::scount_vvec[row_order[ix]], ixs);
            vecRemoveIxs(HTd<encT>::enc_vvec[row_order[ix]], ixs);
#pragma omp atomic update
            to_rm = to_rm - ixs.size();
          } else {
            continue;
          }
        }
        break;
    }
  }
}

template<typename encT>
void
HTs<encT>::shrinkHT(uint64_t num_rm)
{
  uint8_t b = HTs<encT>::b;
  assert(num_rm <= HTs<encT>::num_kmers);
  assert(num_rm < std::numeric_limits<int64_t>::max());
  volatile int64_t to_rm = static_cast<int64_t>(num_rm);
  while (to_rm > 0) {
    std::vector<unsigned int> row_order;
    arrSizeOrder(row_order, HTs<encT>::ind_arr, HTs<encT>::num_rows, true);
    switch (HTs<encT>::rm_strategy) {
      case random_kmer:
#pragma omp parallel for num_threads(num_threads) schedule(static, 1) shared(to_rm)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          int64_t r_to_rm;
#pragma omp atomic read
          r_to_rm = to_rm;
          if (r_to_rm > 0) {
            std::vector<unsigned int> ixs;
            uint8_t n = r_to_rm / HTs<encT>::num_rows + 1;
            getIxsRandom(ixs, HTs<encT>::ind_arr[row_order[ix]], n);
            arrRemoveIxs(HTs<encT>::enc_arr + (row_order[ix] * b), HTs<encT>::ind_arr[row_order[ix]], ixs);
            arrRemoveIxs(HTs<encT>::scount_arr + (row_order[ix] * b), HTs<encT>::ind_arr[row_order[ix]], ixs);
#pragma omp atomic update
            to_rm = to_rm - ixs.size();
          } else {
            continue;
          }
        }
        break;
      case small_scount:
        scT threshold = arrArgmax2D(HTs<encT>::scount_arr, HTs<encT>::ind_arr, HTs<encT>::num_rows, HTs<encT>::b, static_cast<uint64_t>(num_rm));
#pragma omp parallel for num_threads(num_threads) schedule(static, 1) shared(to_rm)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          int64_t r_to_rm;
#pragma omp atomic read
          r_to_rm = to_rm;
          if (r_to_rm > 0) {
            std::vector<unsigned int> ixs;
            arrIxsThreshold(ixs, HTs<encT>::scount_arr + (row_order[ix] * b), HTs<encT>::ind_arr[row_order[ix]], threshold);
            arrRemoveIxs(HTs<encT>::enc_arr + (row_order[ix] * b), HTs<encT>::ind_arr[row_order[ix]], ixs);
            arrRemoveIxs(HTs<encT>::scount_arr + (row_order[ix] * b), HTs<encT>::ind_arr[row_order[ix]], ixs);
#pragma omp atomic update
            to_rm = to_rm - ixs.size();
          } else {
            continue;
          }
        }
        break;
    }
  }
}

template<typename encT>
void
HTd<encT>::convertHTs(HTs<encT>& new_table)
{
  assert(new_table.k == HTd<encT>::k);
  assert(new_table.h == HTd<encT>::h);
  assert(new_table.num_rows == HTd<encT>::num_rows);
  assert(new_table.ptr_lsh_vg == HTd<encT>::ptr_lsh_vg);
  uint8_t b = new_table.b;
  uint32_t num_rows = new_table.num_rows;
  new_table.clearRows();
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (HTd<encT>::enc_vvec[rix].size() <= b) {
      new_table.ind_arr[rix] = HTd<encT>::enc_vvec[rix].size();
      if (new_table.ind_arr[rix]) {
        std::copy(HTd<encT>::enc_vvec[rix].begin(), HTd<encT>::enc_vvec[rix].end(), new_table.enc_arr + (rix * b));
        std::copy(HTd<encT>::scount_vvec[rix].begin(), HTd<encT>::scount_vvec[rix].end(), new_table.scount_arr + (rix * b));
      }
    } else {
      new_table.ind_arr[rix] = b;
      std::copy(HTd<encT>::enc_vvec[rix].begin(), HTd<encT>::enc_vvec[rix].begin() + b, new_table.enc_arr + (rix * b));
      std::copy(HTd<encT>::scount_vvec[rix].begin(), HTd<encT>::scount_vvec[rix].begin() + b, new_table.scount_arr + (rix * b));
    }
  }
  new_table.num_species = HTd<encT>::num_species;
  new_table.updateSize();
}

#include "tableins.cpp"
