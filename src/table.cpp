#include "table.h"

// TODO: Check what is sorted and what is not.
// TODO: Check if size is updated or not.
// TODO: Add macros for debugging to avoid redundant computation.

template<typename T>
inline void
sortColumns(vvec<T>& table)
{
  uint32_t num_rows = table.size();
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!table[rix].empty())
      std::sort(table[rix].begin(), table[rix].end());
  }
}

template<typename encT>
uint64_t
StreamIM<encT>::processInput(uint64_t rbatch_size)
{
  kseq_t* reader = IO::getReader(filepath);
  rbatch_size = IO::adjustBatchSize(rbatch_size, num_threads);
  std::vector<sseq_t> seqBatch = IO::readBatch(reader, rbatch_size);
  uint64_t tnum_kmers_sum = 0;
  uint32_t max_rix = std::numeric_limits<uint32_t>::max();
  max_rix = max_rix >> (32 - 2 * h);
  while (!(seqBatch.empty())) {
    uint32_t pix = lsh_enc_vec.size();
    lsh_enc_vec.resize(lsh_enc_vec.size() + seqBatch.size());
    for (uint32_t ix = 0; ix < seqBatch.size(); ++ix) {
      const char* kmer_seq;
      uint64_t enc64_bp;
      uint64_t enc64_lr;
      uint32_t enc32_lr;
      uint32_t enc32_bp;
      uint32_t rix;
      kmer_seq = seqBatch[ix].nseq.c_str();
      kmerEncodingCompute(kmer_seq, enc64_lr, enc64_bp);
      rix = computeValueLSH(enc64_bp, *(ptr_lsh_vg));
      assert(rix <= max_rix);
      if (std::is_same<encT, uint64_t>::value) {
        lsh_enc_vec[pix + ix] = std::make_pair(rix, enc64_lr);
      } else if (std::is_same<encT, uint32_t>::value) {
        drop64Encoding32(*ptr_npositions, enc64_bp, enc64_lr, enc32_bp, enc32_lr);
        lsh_enc_vec[pix + ix] = std::make_pair(rix, enc32_lr);
      } else {
        std::puts("Available encoding types are 'uint64_t' and 'uint32_t'\n.");
        exit(EXIT_FAILURE);
      }
    }
    tnum_kmers_sum += seqBatch.size();
    if (seqBatch.size() == rbatch_size)
      seqBatch = IO::readBatch(reader, rbatch_size);
    else
      seqBatch.clear();
  }
  lsh_enc_vec.shrink_to_fit();
  std::sort(lsh_enc_vec.begin(),
            lsh_enc_vec.end(),
            [](const std::pair<uint32_t, uint64_t>& l, const std::pair<uint32_t, uint64_t>& r) {
              return l.first < r.first;
            });
  tnum_kmers = tnum_kmers_sum;
  return tnum_kmers;
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
bool
StreamIM<encT>::save(const char* filepath)
{
  bool is_ok = true;
  if (StreamIM<encT>::lsh_enc_vec.empty()) {
    std::puts("The LSH-value and encoding pair vector is empty, nothing to save!\n");
    is_ok = false;
    return is_ok;
  }
  FILE* vec_f = IO::open_file(filepath, is_ok, "wb");
  std::fwrite(StreamIM<encT>::lsh_enc_vec.data(),
              sizeof(std::pair<uint32_t, encT>),
              StreamIM<encT>::lsh_enc_vec.size(),
              vec_f);
  if (std::ferror(vec_f)) {
    std::puts("I/O error when writing LSH-value and encoding pairs.\n");
    is_ok = false;
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
    std::puts("I/O rror when reading LSH-value and encoding pairs.\n");
  } else if (vec_ifs.eof()) {
    is_ok = true;
  }
  vec_ifs.close();
  return is_ok;
}

template<typename encT>
void
StreamIM<encT>::clear()
{
  StreamIM<encT>::lsh_enc_vec.clear();
  l_rix = 0;
  curr_vix = 0;
}

template<typename encT>
uint64_t
StreamIM<encT>::getBatch(vvec<encT>& batch_table, uint32_t tbatch_size)
{
  assert(batch_table.size() >= tbatch_size);
  uint64_t num_kmers = 0;
  while (StreamIM<encT>::curr_vix < StreamIM<encT>::lsh_enc_vec.size()) {
    std::pair<uint32_t, encT> lsh_enc = StreamIM<encT>::lsh_enc_vec[StreamIM<encT>::curr_vix];
    if (lsh_enc.first < (StreamIM<encT>::l_rix + tbatch_size)) {
      batch_table[lsh_enc.first - StreamIM<encT>::l_rix].push_back(lsh_enc.second);
      num_kmers++;
    } else {
      break;
    }
    StreamIM<encT>::curr_vix++;
  }
  StreamIM<encT>::l_rix += tbatch_size;
  sortColumns(batch_table);
  return num_kmers;
}

template<typename encT>
void
StreamOD<encT>::openStream()
{
  is_open = true;
  StreamOD<encT>::vec_ifs = IO::open_ifstream(StreamOD::filepath, is_open);
  if (is_open)
    curr_pos = vec_ifs.tellg();
}

template<typename encT>
uint64_t
StreamOD<encT>::getBatch(vvec<encT>& batch_table, uint32_t tbatch_size, bool contd)
{
  if (!contd)
    assert(batch_table.size() >= tbatch_size);
  uint64_t num_kmers = 0;
  uint64_t row_ix = 0;
  uint32_t curr_rix = StreamOD<encT>::curr_rix;
  uint32_t f_rix = StreamOD<encT>::f_rix;
  const size_t bufsize = 1024 * 1024;
  char buf[bufsize];
  vec_ifs.rdbuf()->pubsetbuf(buf, bufsize);
  while ((row_ix < tbatch_size) && (!StreamOD::vec_ifs.eof() && StreamOD::vec_ifs.good())) {
    std::pair<uint32_t, encT> lsh_enc;
    StreamOD::vec_ifs.read((char*)&lsh_enc, sizeof(std::pair<uint32_t, encT>));
    curr_rix = lsh_enc.first;
    row_ix = curr_rix - f_rix;
    if (!StreamOD::vec_ifs.eof() && (row_ix < tbatch_size)) {
      if (!contd) {
        assert(row_ix < batch_table.size());
        batch_table[row_ix].push_back(lsh_enc.second);
      }
      StreamOD<encT>::curr_rix = curr_rix + 1;
      StreamOD<encT>::curr_pos = StreamOD<encT>::vec_ifs.tellg();
      num_kmers++;
    }
  }
  if (StreamOD<encT>::vec_ifs.fail() && !StreamOD<encT>::vec_ifs.eof()) {
    std::puts("I/O rror when reading LSH-value and encoding pairs.\n");
  } else if (StreamOD<encT>::vec_ifs.eof()) {
    StreamOD<encT>::vec_ifs.close();
  } else {
    StreamOD<encT>::vec_ifs.seekg(curr_pos);
    StreamOD<encT>::f_rix = f_rix + tbatch_size;
  }
  if (!contd)
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
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!enc_vvec[rix].empty()) {
#pragma omp atomic update
      num_kmers_sum += enc_vvec[rix].size();
    }
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
#pragma omp parallel for num_threads(num_threads), schedule(dynamic)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
#pragma omp atomic update
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
    } else {
      scount_vvec[rix].clear();
    }
    tlca_vvec[rix].resize(enc_vvec[rix].size());
    std::fill(tlca_vvec[rix].begin(), tlca_vvec[rix].end(), tID);
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
  HTd<encT>::updateSize();
#ifdef DEBUG
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
      switch (kmer_ranking) {
        case random_kmer:
          getIxsRandom(ixs, enc_vvec[rix].size(), enc_vvec[rix].size() - b_max);
          break;
        case large_scount:
          vecIxsNumber(ixs, scount_vvec[rix], enc_vvec[rix].size() - b_max, true);
          break;
        case information_score:
          std::unordered_map<encT, std::vector<scT>> values_map =
            mapValuesCountsHTd(childrenHT, rix);
          std::vector<scT> scores_vec;
          vecInformationScores(scores_vec, enc_vvec[rix], values_map);
          vecIxsNumber(ixs, scores_vec, enc_vvec[rix].size() - b_max, true);
          break;
      }
      vecRemoveIxs(enc_vvec[rix], ixs);
      vecRemoveIxs(scount_vvec[rix], ixs);
      vecRemoveIxs(tlca_vvec[rix], ixs);
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
void
HTd<encT>::unionRows(HTd<encT>& sibling, bool update_size)
{
#ifdef DEBUG
  if (!HTd<encT>::areColumnsSorted()) {
    HTd<encT>::sortColumns();
    std::puts("HTd (main) has unsorted columns before computing the union of rows.\n");
  }
  if (!sibling.areColumnsSorted()) {
    sibling.sortColumns();
    std::puts("HTd (sibling) has unsorted columns before computing the union of rows.\n");
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
    std::puts("HTs (main) has unsorted columns before computing the union of rows.\n");
  }
  if (!sibling.areColumnsSorted()) {
    sibling.sortColumns();
    std::puts("HTs (sibling) has unsorted columns before computing the union of rows.\n");
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
        switch (kmer_ranking) {
          case random_kmer: {
            getIxsRandom(ixs, v.size(), v.size() - b);
            break;
          }
          case large_scount: {
            std::vector<scT> s_v;
            std::transform(v.begin(), v.end(), back_inserter(s_v), [&count_map](encT val) {
              return count_map[val];
            });
            vecIxsNumber(ixs, s_v, v.size() - b, true);
            break;
          }
          case information_score: {
            std::unordered_map<encT, std::vector<scT>> values_map =
              mapValuesCountsHTs(childrenHT, rix);
            std::vector<scT> scores_vec;
            vecInformationScores(scores_vec, v, values_map);
            vecIxsNumber(ixs, scores_vec, v.size() - b, true);
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
        switch (kmer_ranking) {
          case random_kmer: {
            getIxsRandom(ixs, v.size(), v.size() - b);
            break;
          }
          case large_scount: {
            std::vector<scT> s_v;
            std::transform(v.begin(), v.end(), back_inserter(s_v), [&count_map](encT val) {
              return count_map[val];
            });
            vecIxsNumber(ixs, s_v, v.size() - b, true);
            break;
          }
          case information_score: {
            std::unordered_map<encT, std::vector<scT>> values_map =
              mapValuesCountsHTs(childrenHT, rix);
            std::vector<scT> scores_vec;
            vecInformationScores(scores_vec, v, values_map);
            vecIxsNumber(ixs, scores_vec, v.size() - b, true);
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
    switch (kmer_ranking) {
      case random_kmer: {
        uint8_t n = static_cast<uint64_t>(to_rm) / num_rows + 1;
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          uint32_t& rix = row_order[ix];
          if (scount_vvec[rix].size() >= n) {
            int64_t r_to_rm;
#pragma omp atomic read
            r_to_rm = to_rm;
            if (r_to_rm > 0) {
              std::vector<unsigned int> ixs;
              getIxsRandom(ixs, scount_vvec[rix].size(), n);
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
      case large_scount: {
        scT threshold = vvecArgmax2D(scount_vvec, static_cast<uint64_t>(to_rm), true);
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          uint32_t& rix = row_order[ix];
          if (scount_vvec[rix].size() > 0) {
            int64_t r_to_rm;
#pragma omp atomic read
            r_to_rm = to_rm;
            if (r_to_rm > 0) {
              std::vector<unsigned int> ixs;
              vecIxsThreshold(ixs, scount_vvec[rix], threshold, true);
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
        std::map<scT, uint64_t> val_counts{};
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t rix = 0; rix < num_rows; ++rix) {
          if (enc_vvec[rix].size() > 0) {
            std::unordered_map<encT, std::vector<scT>> values_map =
              mapValuesCountsHTd(childrenHT, rix);
            std::vector<scT> scores_vec;
            vecInformationScores(scores_vec, enc_vvec[rix], values_map);
#pragma omp critical
            {
              for (auto& val : scores_vec)
                val_counts[val]++;
            }
          }
        }
        scT threshold = mapArgmax(val_counts, static_cast<uint64_t>(to_rm), true);
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          uint32_t& rix = row_order[ix];
          if (enc_vvec[rix].size() > 0) {
            int64_t r_to_rm;
#pragma omp atomic read
            r_to_rm = to_rm;
            if (r_to_rm > 0) {
              std::unordered_map<encT, std::vector<scT>> values_map =
                mapValuesCountsHTd(childrenHT, rix);
              std::vector<scT> scores_vec;
              vecInformationScores(scores_vec, enc_vvec[rix], values_map);
              std::vector<unsigned int> ixs;
              vecIxsThreshold(ixs, scores_vec, threshold, true);
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
    switch (kmer_ranking) {
      case random_kmer: {
        uint8_t n = to_rm / num_rows + 1;
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
      case large_scount: {
        scT threshold =
          arrArgmax2D(scount_arr, ind_arr, num_rows, b, static_cast<uint64_t>(to_rm), true);
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          uint32_t& rix = row_order[ix];
          if (ind_arr[rix] > 0) {
            int64_t r_to_rm;
#pragma omp atomic read
            r_to_rm = to_rm;
            if (r_to_rm > 0) {
              std::vector<unsigned int> ixs;
              arrIxsThreshold(ixs, scount_arr + (rix * b), ind_arr[rix], threshold, true);
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
        std::map<scT, uint64_t> val_counts{};
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t rix = 0; rix < num_rows; ++rix) {
          if (ind_arr[rix] > 0) {
            std::unordered_map<encT, std::vector<scT>> values_map =
              mapValuesCountsHTs(childrenHT, rix);
            std::vector<scT> scores_vec;
            arrInformationScores(scores_vec, enc_arr + rix * b, ind_arr[rix], values_map);
#pragma omp critical
            {
              for (auto& val : scores_vec)
                val_counts[val]++;
            }
          }
        }
        scT threshold = mapArgmax(val_counts, static_cast<uint64_t>(to_rm), true);
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint32_t ix = 0; ix < row_order.size(); ++ix) {
          uint32_t& rix = row_order[ix];
          if (ind_arr[rix] > 0) {
            int64_t r_to_rm;
#pragma omp atomic read
            r_to_rm = to_rm;
            if (r_to_rm > 0) {
              std::unordered_map<encT, std::vector<scT>> values_map =
                mapValuesCountsHTs(childrenHT, rix);
              std::vector<scT> scores_vec;
              arrInformationScores(scores_vec, enc_arr + rix * b, ind_arr[rix], values_map);
              std::vector<unsigned int> ixs;
              vecIxsThreshold(ixs, scores_vec, threshold, true);
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
          /* std::unordered_map<encT, std::vector<float>> select_map =
           * mapValuesSelectHTd(childrenHT, rix); */
          /* std::bernoulli_distribution d(updateLCAtProbabilityAPN(prob_map[enc_vvec[rix][i]],
           * select_map[enc_vvec[rix][i]])); */
          std::bernoulli_distribution d(updateLCAtProbabilityC2N(prob_map[enc_vvec[rix][i]]));
          if (d(gen)) {
            tlca_vvec[rix][i] = tID;
          } else {
            std::discrete_distribution<> d(prob_map[enc_vvec[rix][i]].begin(),
                                           prob_map[enc_vvec[rix][i]].end());
            for (auto v : tlca_map[enc_vvec[rix][i]])
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
          /* std::unordered_map<encT, std::vector<float>> select_map =
           * mapValuesSelectHTs(childrenHT, rix); */
          /* std::bernoulli_distribution d(updateLCAtProbabilityAPN(prob_map[enc_arr[rix * b + i]],
           * select_map[enc_arr[rix * b + i]])); */
          std::bernoulli_distribution d(updateLCAtProbabilityC2N(prob_map[enc_arr[rix * b + i]]));
          if (d(gen)) {
            tlca_arr[rix * b + i] = tID;
          } else {
            std::discrete_distribution<> d(prob_map[enc_arr[rix * b + i]].begin(),
                                           prob_map[enc_arr[rix * b + i]].end());
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
