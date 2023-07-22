#include "table.h"

#define SIZE_ESTIMATE 4194304

// TODO: Benchmark different thread schedules etc.
// TODO: Check what is sorted and what is not.
// TODO: Check if size is updated or not.
// TODO: Add macros for debugging to avoid redundant computation.

template<typename T>
inline void
sortColumns(std::vector<std::vector<T>>& table)
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
tableC<encT>::fill(const char* filepath, unsigned int batch_size)
{
  kseq_t* reader = getReader(filepath);
  batch_size = adjustBatchSize(batch_size, num_threads);
  std::vector<sseq_t> seqBatch = readBatch(reader, batch_size);
  uint64_t total_kmer = 0;
  uint32_t max_rix = pow(2, 2 * tableC<encT>::h);
  tableC<encT>::rix_enc_vec.reserve(SIZE_ESTIMATE);
  while (!(seqBatch.empty())) {
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
    for (uint32_t ix = 0; ix < (uint32_t)seqBatch.size(); ++ix) {
      const char* kmer_seq;
      uint64_t enc_bp;
      uint64_t enc_lr;
      uint32_t rix;
      kmer_seq = seqBatch[ix].nseq.c_str();
      kmerEncodingCompute(kmer_seq, enc_lr, enc_bp);
      rix = computeValueLSH(enc_bp, *(tableC<encT>::ptr_lsh_vg));
      assert(rix < max_rix);
#pragma omp critical
      {
        tableC<encT>::rix_enc_vec.push_back(std::make_pair(rix, enc_lr));
      }
    }
    total_kmer += seqBatch.size();
    if (seqBatch.size() == batch_size)
      seqBatch = readBatch(reader, batch_size);
    else
      seqBatch.clear();
  }
  tableC<encT>::rix_enc_vec.shrink_to_fit();
  std::sort(tableC<encT>::rix_enc_vec.begin(),
            tableC<encT>::rix_enc_vec.end(),
            [](const std::pair<uint32_t, uint64_t>& left, const std::pair<uint32_t, uint64_t>& right) { return left.first < right.first; });
  return total_kmer;
}

template<typename encT>
std::unordered_map<uint8_t, uint64_t>
tableC<encT>::histNumCols()
{
  std::unordered_map<uint32_t, uint8_t> row_sizes;
  for (auto kv : tableC<encT>::rix_enc_vec) {
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
tableD<encT>::histNumCols()
{
  std::unordered_map<uint8_t, uint64_t> hist_map;
  for (auto& row : tableD<encT>::enc_vvec) {
    hist_map[row.size()]++;
  }
  return hist_map;
}

template<typename encT>
std::unordered_map<uint8_t, uint64_t>
tableF<encT>::histNumCols()
{
  std::unordered_map<uint8_t, uint64_t> hist_map;
  uint32_t num_rows = tableF<encT>::num_rows;
  for (unsigned int rix = 0; rix < num_rows; ++rix) {
    hist_map[tableF<encT>::ind_arr[rix]]++;
  }
  return hist_map;
}

template<typename encT>
bool
tableC<encT>::save(const char* filepath)
{
  bool is_ok = false;
  if (tableC<encT>::rix_enc_vec.empty()) {
    std::puts("The LSH-value and encoding pair vector is empty, nothing to save!");
    return is_ok;
  }
  FILE* vec_f = open_file(filepath, is_ok, "wb");
  std::fwrite(tableC<encT>::rix_enc_vec.data(), tableC<encT>::rix_enc_vec.size(), sizeof(std::pair<uint32_t, encT>), vec_f);
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
tableC<encT>::load(const char* filepath)
{
  bool is_ok = false;
  tableC<encT>::rix_enc_vec.clear();
  std::ifstream vec_ifs = open_ifstream(filepath, is_ok);
  while (!vec_ifs.eof() && vec_ifs.good()) {
    std::pair<uint32_t, encT> lsh_enc;
    vec_ifs.read((char*)&lsh_enc, sizeof(std::pair<uint32_t, encT>));
    if (!vec_ifs.eof())
      tableC<encT>::rix_enc_vec.push_back(lsh_enc);
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
tableC<encT>::getBatch(std::vector<std::vector<encT>>& batch_table, uint32_t batch_size)
{
  assert(batch_table.size() >= batch_size);
  uint64_t num_kmers = 0;
  while (tableC<encT>::curr_vix < tableC<encT>::rix_enc_vec.size()) {
    std::pair<uint32_t, uint64_t> lsh_enc = tableC<encT>::rix_enc_vec[tableC<encT>::curr_vix];
    if (lsh_enc.first < (tableC<encT>::l_rix + batch_size)) {
      batch_table[lsh_enc.first - tableC<encT>::l_rix].push_back(lsh_enc.second);
      num_kmers++;
    } else {
      break;
    }
    tableC<encT>::curr_vix++;
  }
  tableC<encT>::l_rix += batch_size;
  sortColumns(batch_table);
  return num_kmers;
}

template<typename encT>
bool
tableS<encT>::openStream()
{
  bool is_ok = true;
  tableS<encT>::vec_ifs = open_ifstream(tableS::filepath, is_ok);
  return is_ok;
}

template<typename encT>
uint64_t
tableS<encT>::getBatch(std::vector<std::vector<encT>>& batch_table, uint32_t batch_size)
{
  assert(batch_table.size() >= batch_size);
  uint64_t num_kmers = 0;
  uint64_t row_ix = 0;
  uint32_t curr_rix = tableS<encT>::curr_rix;
  uint32_t f_rix = tableS<encT>::f_rix;
  while ((row_ix < batch_size) && (!tableS::vec_ifs.eof() && tableS::vec_ifs.good())) {
    std::pair<uint32_t, encT> lsh_enc;
    tableS::vec_ifs.read((char*)&lsh_enc, sizeof(std::pair<uint32_t, encT>));
    curr_rix = lsh_enc.first;
    row_ix = curr_rix - f_rix;
    if (!tableS::vec_ifs.eof() && (row_ix < batch_size)) {
      assert(row_ix < batch_table.size());
      batch_table[row_ix].push_back(lsh_enc.second);
      tableS<encT>::curr_rix = curr_rix + 1;
      tableS<encT>::curr_pos = tableS<encT>::vec_ifs.tellg();
      num_kmers++;
    }
  }
  if (tableS<encT>::vec_ifs.fail() && !tableS<encT>::vec_ifs.eof()) {
    std::puts("I/Oerror when reading LSH-value and encoding pairs.");
  } else if (tableS<encT>::vec_ifs.eof()) {
    tableS<encT>::vec_ifs.close();
  } else {
    tableS<encT>::vec_ifs.seekg(curr_pos);
    tableS<encT>::f_rix = f_rix + batch_size;
  }
  sortColumns(batch_table);
  return num_kmers;
}

template<typename encT>
void
tableD<encT>::makeUnique(bool update_size)
{
  if (!tableD<encT>::areColumnsSorted())
    tableD<encT>::sortColumns();
  uint32_t num_rows = tableD<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!tableD<encT>::enc_vvec[rix].empty())
      tableD<encT>::enc_vvec[rix].erase(std::unique(tableD<encT>::enc_vvec[rix].begin(), tableD<encT>::enc_vvec[rix].end()), tableD<encT>::enc_vvec[rix].end());
  }
  if (update_size)
    tableD<encT>::updateSize();
}

template<typename encT>
void
tableF<encT>::makeUnique(bool update_size)
{
  if (!tableF<encT>::areColumnsSorted())
    tableF<encT>::sortColumns();
  uint8_t b = tableF<encT>::b;
  uint32_t num_rows = tableF<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (tableF<encT>::ind_arr[rix] > 0) {
      auto last =
        std::unique(tableF<encT>::enc_arr + (rix * b), tableF<encT>::enc_arr + (rix * b) + tableF<encT>::ind_arr[rix]) - tableF<encT>::enc_arr - (rix * b);
      tableF<encT>::ind_arr[rix] = last;
    }
  }
  if (update_size)
    tableF<encT>::updateSize();
}

template<typename encT>
void
tableD<encT>::sortColumns()
{
  uint32_t num_rows = tableD<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!tableD<encT>::enc_vvec[rix].empty())
      std::sort(tableD<encT>::enc_vvec[rix].begin(), tableD<encT>::enc_vvec[rix].end());
  }
}

template<typename encT>
void
tableF<encT>::sortColumns()
{
  uint8_t b = tableF<encT>::b;
  uint32_t num_rows = tableF<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (tableF<encT>::ind_arr[rix] > 0) {
      std::sort(tableF<encT>::enc_arr + (rix * b), tableF<encT>::enc_arr + (rix * b) + tableF<encT>::ind_arr[rix]);
    }
  }
}

template<typename encT>
bool
tableD<encT>::areColumnsSorted()
{
  bool is_ok = true;
  uint32_t num_rows = tableD<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!tableD<encT>::enc_vvec[rix].empty() && !std::is_sorted(tableD<encT>::enc_vvec[rix].begin(), tableD<encT>::enc_vvec[rix].end())) {
#pragma omp atomic write
      is_ok = false;
    }
  }
  return is_ok;
}

template<typename encT>
bool
tableF<encT>::areColumnsSorted()
{
  bool is_ok = true;
  uint8_t b = tableF<encT>::b;
  uint32_t num_rows = tableF<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if ((tableF<encT>::ind_arr[rix] > 0) &&
        !std::is_sorted(tableF<encT>::enc_arr + (rix * b), tableF<encT>::enc_arr + (rix * b) + tableF<encT>::ind_arr[rix])) {
#pragma omp atomic write
      is_ok = false;
    }
  }
  return is_ok;
}

template<typename encT>
void
tableD<encT>::clearRows()
{
  tableD<encT>::enc_vvec.clear();
  tableD<encT>::enc_vvec.resize(tableD<encT>::num_rows);
  tableD<encT>::updateSize();
}

template<typename encT>
void
tableF<encT>::clearRows()
{
  uint32_t num_rows = tableF<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    tableF<encT>::ind_arr[rix] = 0;
  }
}

template<typename encT>
void
tableD<encT>::updateSize()
{
  uint64_t num_kmers = 0;
  for (auto& row : enc_vvec) {
    if (!row.empty())
      num_kmers += row.size();
  }
  tableD<encT>::num_kmers = num_kmers;
  tableD<encT>::table_size = tableD<encT>::num_kmers * sizeof(encT) + tableD<encT>::num_rows * sizeof(std::vector<encT>);
}

template<typename encT>
void
tableF<encT>::updateSize()
{
  uint64_t num_kmers = 0;
  uint32_t num_rows = tableF<encT>::num_rows;
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    num_kmers += tableF<encT>::ind_arr[rix];
  }
  tableF<encT>::num_kmers = num_kmers;
  tableF<encT>::table_size = tableF<encT>::table_size;
}

template<typename encT>
void
tableD<encT>::trimColumns(uint8_t b)
{
  uint32_t num_rows = tableD<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (tableD<encT>::enc_vvec[rix].size() > b)
      tableD<encT>::enc_vvec[rix].resize(b);
  }
}

template<typename encT>
void
tableD<encT>::pruneColumns(uint8_t b)
{
  uint32_t num_rows = tableD<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (tableD<encT>::enc_vvec[rix].size() > b)
      pruneVectorPseudorandom(tableD<encT>::enc_vvec[rix], b);
  }
}

template<typename encT>
void
tableD<encT>::unionRows(tableD<encT>& sibling, bool update_size)
{
  if (!tableD<encT>::areColumnsSorted())
    tableD<encT>::sortColumns();
  if (!sibling.areColumnsSorted())
    sibling.sortColumns();
  assert(tableD<encT>::num_rows == sibling.num_rows);
  uint32_t num_rows = tableD<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!tableD<encT>::enc_vvec[rix].empty() && !sibling.enc_vvec[rix].empty()) {
      tableD<encT>::enc_vvec[rix].insert(tableD<encT>::enc_vvec[rix].end(), sibling.enc_vvec[rix].begin(), sibling.enc_vvec[rix].end());
      std::inplace_merge(tableD<encT>::enc_vvec[rix].begin(),
                         tableD<encT>::enc_vvec[rix].begin() + (tableD<encT>::enc_vvec[rix].size() - sibling.enc_vvec[rix].size()),
                         tableD<encT>::enc_vvec[rix].end());
      tableD<encT>::enc_vvec[rix].erase(std::unique(tableD<encT>::enc_vvec[rix].begin(), tableD<encT>::enc_vvec[rix].end()), tableD<encT>::enc_vvec[rix].end());
    } else if (!sibling.enc_vvec[rix].empty()) {
      tableD<encT>::enc_vvec[rix] = sibling.enc_vvec[rix];
    }
  }
  if (update_size)
    tableD<encT>::updateSize();
}

template<typename encT>
void
tableF<encT>::unionRows(tableF<encT>& sibling, bool update_size)
{
  if (!tableF<encT>::areColumnsSorted())
    tableF<encT>::sortColumns();
  if (!sibling.areColumnsSorted())
    sibling.sortColumns();
  assert(tableF<encT>::b == sibling.b);
  assert(tableF<encT>::num_rows * tableF<encT>::b == sibling.num_rows * sibling.b);
  uint32_t b = tableF<encT>::b;
  uint32_t num_rows = tableF<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    unsigned int ix = rix * b;
    if ((sibling.ind_arr[rix] > 0) && (tableF<encT>::ind_arr[rix] > 0)) {
      std::list<encT> l;
      std::set_union(tableF<encT>::enc_arr + ix,
                     tableF<encT>::enc_arr + (ix + tableF<encT>::ind_arr[rix]),
                     sibling.enc_arr + ix,
                     sibling.enc_arr + (ix + sibling.ind_arr[rix]),
                     std::back_inserter(l));
      if (l.size() > b)
        pruneListPseudorandom(l, b);
      std::copy(l.begin(), l.end(), tableF<encT>::enc_arr + ix);
      tableF<encT>::ind_arr[rix] = l.size();
    } else if (sibling.ind_arr[rix] > 0) {
      std::copy(sibling.enc_arr + ix, sibling.enc_arr + (ix + sibling.ind_arr[rix]), tableF<encT>::enc_arr + ix);
      tableF<encT>::ind_arr[rix] = sibling.ind_arr[rix];
    }
  }
  if (update_size)
    tableF<encT>::updateSize();
}

template<typename encT>
void
tableD<encT>::mergeRows(tableD<encT>& sibling, bool update_size)
{
  if (!tableD<encT>::areColumnsSorted())
    tableD<encT>::sortColumns();
  if (!sibling.areColumnsSorted())
    sibling.sortColumns();
  assert(tableD<encT>::num_rows == sibling.num_rows);
  uint32_t num_rows = tableD<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (!tableD<encT>::enc_vvec[rix].empty() && !sibling.enc_vvec[rix].empty()) {
      tableD<encT>::enc_vvec[rix].insert(tableD<encT>::enc_vvec[rix].end(), sibling.enc_vvec[rix].begin(), sibling.enc_vvec[rix].end());
      std::inplace_merge(tableD<encT>::enc_vvec[rix].begin(),
                         tableD<encT>::enc_vvec[rix].begin() + (tableD<encT>::enc_vvec[rix].size() - sibling.enc_vvec[rix].size()),
                         tableD<encT>::enc_vvec[rix].end());
    } else if (!sibling.enc_vvec[rix].empty()) {
      tableD<encT>::enc_vvec[rix] = sibling.enc_vvec[rix];
    }
  }
  if (update_size)
    tableD<encT>::updateSize();
}

template<typename encT>
void
tableF<encT>::mergeRows(tableF<encT>& sibling, bool update_size)
{
  if (!tableF<encT>::areColumnsSorted())
    tableF<encT>::sortColumns();
  if (!sibling.areColumnsSorted())
    sibling.sortColumns();
  assert(tableF<encT>::b == sibling.b);
  assert(tableF<encT>::num_rows * tableF<encT>::b == sibling.num_rows * sibling.b);
  uint32_t b = tableF<encT>::b;
  uint32_t num_rows = tableF<encT>::num_rows;
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    unsigned int ix = rix * b;
    if ((sibling.ind_arr[rix] > 0) && (tableF<encT>::ind_arr[rix] > 0)) {
      std::list<encT> l;
      std::merge(tableF<encT>::enc_arr + ix,
                 tableF<encT>::enc_arr + (ix + tableF<encT>::ind_arr[rix]),
                 sibling.enc_arr + ix,
                 sibling.enc_arr + (ix + sibling.ind_arr[rix]),
                 std::back_inserter(l));
      if (l.size() > b)
        l.erase(std::unique(l.begin(), l.end()), l.end());
      if (l.size() > b)
        pruneListPseudorandom(l, b);
      std::copy(l.begin(), l.end(), tableF<encT>::enc_arr + ix);
      tableF<encT>::ind_arr[rix] = l.size();
    } else if (sibling.ind_arr[rix] > 0) {
      std::copy(sibling.enc_arr + ix, sibling.enc_arr + (ix + sibling.ind_arr[rix]), tableF<encT>::enc_arr + ix);
      tableF<encT>::ind_arr[rix] = sibling.ind_arr[rix];
    }
  }
  if (update_size)
    tableF<encT>::updateSize();
}

template<typename encT>
void
tableD<encT>::removeIndices(std::vector<std::pair<uint32_t, uint8_t>>& indices_vec)
{
  // TODO: PARALLEL
  // TODO: USE A MAP INSTEAD OF VECTOR OF PAIRS
  for (auto kv : indices_vec) {
    if (!tableD<encT>::enc_vvec[kv.first].empty() && (kv.second < tableD<encT>::enc_vvec[kv.first].size()))
      tableD<encT>::enc_vvec[kv.first].erase(tableD<encT>::enc_vvec[kv.first].begin() + kv.second);
    else
      std::puts("The attempt of removing non-existent index has been ignored.");
  }
}

template<typename encT>
void
tableF<encT>::removeIndices(std::vector<std::pair<uint32_t, uint8_t>>& indices_vec)
{
  uint8_t b = tableF<encT>::b;
  uint32_t num_rows = tableF<encT>::num_rows;
  // TODO: PARALLEL
  // TODO: USE A MAP INSTEAD OF VECTOR OF PAIRS
  for (auto kv : indices_vec) {
    if (!tableF<encT>::ind_arr[kv.first] && (kv.second < tableF<encT>::ind_arr[kv.first])) {
      auto ix = tableF<encT>::enc_arr + (kv.first * b);
      std::copy(ix + kv.second + 1, ix + tableF<encT>::ind_arr[kv.first], ix + kv.second);
      tableF<encT>::ind_arr[kv.first]--;
    } else {
      std::puts("The attempt of removing non-existent index has been ignored.");
    }
  }
}

template<typename encT>
void
tableD<encT>::transformF(tableF<encT>& new_table)
{
  new_table.clearRows();
  assert(new_table.k == tableD<encT>::k);
  assert(new_table.h == tableD<encT>::h);
  assert(new_table.num_rows == tableD<encT>::num_rows);
  assert(new_table.ptr_lsh_vg == tableD<encT>::ptr_lsh_vg);
  uint8_t b = new_table.b;
  uint32_t num_rows = new_table.num_rows;
  for (uint32_t rix = 0; rix < num_rows; ++rix) {
    if (tableD<encT>::enc_vvec[rix].size() <= b) {
      new_table.ind_arr[rix] = tableD<encT>::enc_vvec[rix].size();
      if (new_table.ind_arr[rix])
        std::copy(tableD<encT>::enc_vvec[rix].begin(), tableD<encT>::enc_vvec[rix].end(), new_table.enc_arr + (rix * b));
    } else {
      std::vector<encT> row = tableD<encT>::enc_vvec[rix];
      pruneVectorPseudorandom(row, b);
      new_table.ind_arr[rix] = b;
      std::copy(row.begin(), row.end(), new_table.enc_arr + (rix * b));
    }
  }
  new_table.updateSize();
}

template uint64_t
tableC<uint32_t>::fill(const char* filepath, unsigned int batch_size);

template uint64_t
tableC<uint64_t>::fill(const char* filepath, unsigned int batch_size);

template bool
tableC<uint64_t>::save(const char* filepath);

template bool
tableC<uint32_t>::save(const char* filepath);

template bool
tableC<uint64_t>::load(const char* filepath);

template bool
tableC<uint32_t>::load(const char* filepath);

template uint64_t
tableC<uint64_t>::getBatch(std::vector<std::vector<u_int64_t>>& batch_table, uint32_t batch_size);

template uint64_t
tableC<uint32_t>::getBatch(std::vector<std::vector<u_int32_t>>& batch_table, uint32_t batch_size);

template bool
tableS<uint64_t>::openStream();

template bool
tableS<uint32_t>::openStream();

template uint64_t
tableS<uint64_t>::getBatch(std::vector<std::vector<u_int64_t>>& batch_table, uint32_t batch_size);

template uint64_t
tableS<uint32_t>::getBatch(std::vector<std::vector<u_int32_t>>& batch_table, uint32_t batch_size);

template std::unordered_map<uint8_t, uint64_t>
tableC<uint64_t>::histNumCols();

template std::unordered_map<uint8_t, uint64_t>
tableC<uint32_t>::histNumCols();

template void
tableD<uint32_t>::makeUnique(bool update_size);

template void
tableD<uint64_t>::makeUnique(bool update_size);

template void
tableF<uint32_t>::makeUnique(bool update_size);

template void
tableF<uint64_t>::makeUnique(bool update_size);

template void
tableD<uint32_t>::trimColumns(uint8_t b);

template void
tableD<uint64_t>::trimColumns(uint8_t b);

template void
tableD<uint32_t>::pruneColumns(uint8_t b);

template void
tableD<uint64_t>::pruneColumns(uint8_t b);

template void
tableD<uint64_t>::sortColumns();

template void
tableD<uint32_t>::sortColumns();

template void
tableF<uint64_t>::sortColumns();

template void
tableF<uint32_t>::sortColumns();

template bool
tableD<uint64_t>::areColumnsSorted();

template bool
tableD<uint32_t>::areColumnsSorted();

template bool
tableF<uint64_t>::areColumnsSorted();

template bool
tableF<uint32_t>::areColumnsSorted();

template void
tableD<uint32_t>::clearRows();

template void
tableD<uint64_t>::clearRows();

template void
tableF<uint32_t>::clearRows();

template void
tableF<uint64_t>::clearRows();

template void
tableD<uint64_t>::updateSize();

template void
tableD<uint32_t>::updateSize();

template void
tableF<uint64_t>::updateSize();

template void
tableF<uint32_t>::updateSize();

template void
tableD<uint64_t>::unionRows(tableD<uint64_t>& sibling, bool update_size);

template void
tableD<uint32_t>::unionRows(tableD<uint32_t>& sibling, bool update_size);

template void
tableF<uint64_t>::unionRows(tableF<uint64_t>& sibling, bool update_size);

template void
tableF<uint32_t>::unionRows(tableF<uint32_t>& sibling, bool update_size);

template void
tableD<uint64_t>::mergeRows(tableD<uint64_t>& sibling, bool update_size);

template void
tableD<uint32_t>::mergeRows(tableD<uint32_t>& sibling, bool update_size);

template void
tableF<uint64_t>::mergeRows(tableF<uint64_t>& sibling, bool update_size);

template void
tableF<uint32_t>::mergeRows(tableF<uint32_t>& sibling, bool update_size);

template void
tableD<uint64_t>::removeIndices(std::vector<std::pair<uint32_t, uint8_t>>& indices_vec);

template void
tableD<uint32_t>::removeIndices(std::vector<std::pair<uint32_t, uint8_t>>& indices_vec);

template void
tableF<uint64_t>::removeIndices(std::vector<std::pair<uint32_t, uint8_t>>& indices_vec);

template void
tableF<uint32_t>::removeIndices(std::vector<std::pair<uint32_t, uint8_t>>& indices_vec);

template std::unordered_map<uint8_t, uint64_t>
tableD<uint64_t>::histNumCols();

template std::unordered_map<uint8_t, uint64_t>
tableD<uint32_t>::histNumCols();

template std::unordered_map<uint8_t, uint64_t>
tableF<uint64_t>::histNumCols();

template std::unordered_map<uint8_t, uint64_t>
tableF<uint32_t>::histNumCols();

template void
tableD<uint32_t>::transformF(tableF<uint32_t>& table);

template void
tableD<uint64_t>::transformF(tableF<uint64_t>& table);
