#include "table.h"
#include <cstdint>
#include <sys/types.h>
#include <utility>

template<typename encT>
uint64_t
tableC<encT>::fillVec(const char* filepath, unsigned int batch_size, unsigned int num_threads)
{
  kseq_t* reader = getReader(filepath);
  batch_size = adjustBatchSize(batch_size, num_threads);
  std::vector<sseq_t> seqBatch = readBatch(reader, batch_size);
  uint64_t total_kmer = 0;

  /* tableC<encT>::rix_enc_vec.reserve(pow(2, 22)); */

  while (!(seqBatch.empty())) {
    /* tableC<encT>::rix_enc_vec.reserve(tableC<encT>::rix_enc_vec.capacity() + seqBatch.size()); */
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
    for (uint32_t ix = 0; ix < (uint32_t)seqBatch.size(); ++ix) {
      const char* kmer_seq;
      uint64_t enc_bp;
      uint64_t enc_lr;
      uint32_t rix;

      kmer_seq = seqBatch[ix].nseq.c_str();
      kmerEncodingCompute(kmer_seq, enc_lr, enc_bp);
      rix = computeValueLSH(enc_bp, *(tableC<encT>::ptr_lsh_vg));
      assert(rix < tableC<encT>::num_rows);
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

  /* tableC<encT>::rix_enc_vec.shrink_to_fit(); */
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
  std::unordered_map<uint8_t, uint64_t> hist_vec;
  for (auto kv : row_sizes) {
    hist_vec[kv.second]++;
  }
  return hist_vec;
}

template<typename encT>
bool
tableC<encT>::saveVec(const char* filepath)
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
tableC<encT>::loadVec(const char* filepath)
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
  uint64_t num_kmers = 0;

  assert(batch_table.size() >= batch_size);

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
  uint64_t num_kmers = 0;
  uint64_t num_rows = 0;
  uint32_t curr_rix = tableS<encT>::curr_rix;
  uint32_t f_rix = tableS<encT>::f_rix;

  assert(batch_table.size() >= batch_size);

  while ((num_rows < batch_size) && (!tableS::vec_ifs.eof() && tableS::vec_ifs.good())) {
    std::pair<uint32_t, encT> lsh_enc;
    tableS::vec_ifs.read((char*)&lsh_enc, sizeof(std::pair<uint32_t, encT>));
    curr_rix = lsh_enc.first;
    num_rows = curr_rix - f_rix;
    if (!tableS::vec_ifs.eof() && (num_rows < batch_size)) {
      assert(num_rows < batch_table.size());
      batch_table[num_rows].push_back(lsh_enc.second);
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

  return num_kmers;
}

template<typename encT>
void
tableD<encT>::makeUnique()
{
  for (auto& row : dynamic_table) {
    {
      std::sort(row.begin(), row.end());
      row.erase(std::unique(row.begin(), row.end()), row.end());
    }
  }
}

template<typename encT>
void
tableD<encT>::pruneColumns(uint8_t b)
{
  for (auto& row : dynamic_table) {
    if (row.size() > b)
      row.resize(b);
  }
}

template uint64_t
tableC<uint32_t>::fillVec(const char* filepath, unsigned int batch_size, unsigned int num_threads);

template uint64_t
tableC<uint64_t>::fillVec(const char* filepath, unsigned int batch_size, unsigned int num_threads);

template bool
tableC<uint64_t>::saveVec(const char* filepath);

template bool
tableC<uint32_t>::saveVec(const char* filepath);

template bool
tableC<uint64_t>::loadVec(const char* filepath);

template bool
tableC<uint32_t>::loadVec(const char* filepath);

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
tableD<uint32_t>::makeUnique();

template void
tableD<uint64_t>::makeUnique();

template void
tableD<uint32_t>::pruneColumns(uint8_t b);

template void
tableD<uint64_t>::pruneColumns(uint8_t b);
