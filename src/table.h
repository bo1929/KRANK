#ifndef _TABLE_H
#define _TABLE_H

#include "assess.h"
#include "common.h"
#include "encode.h"
#include "io.h"
#include "lsh.h"

template<typename encT>
struct tableC
{
  uint8_t k, h;
  uint32_t l_rix;
  uint32_t curr_vix;
  maskLSH* ptr_lsh_vg;
  std::vector<std::pair<uint32_t, encT>> rix_enc_vec;

  tableC(uint8_t k, uint8_t h, maskLSH* ptr_lsh_vg)
    : k(h)
    , h(h)
    , l_rix(0)
    , curr_vix(0)
    , ptr_lsh_vg(ptr_lsh_vg)
  {
  }
  bool saveVec(const char* filepath);
  bool loadVec(const char* filepath);
  uint64_t fillVec(const char* filepath, unsigned int batch_size, unsigned int num_threads);
  uint64_t getBatch(std::vector<std::vector<encT>>& batch_table, uint32_t batch_size);
  std::unordered_map<uint8_t, uint64_t> histNumCols();
};

template<typename encT>
struct tableS
{
  uint32_t f_rix;
  uint32_t curr_rix;
  const char* filepath;
  std::ifstream vec_ifs;
  std::streampos curr_pos;
  bool is_open;

  tableS(const char* filepath)
    : filepath(filepath)
    , curr_rix(0)
    , f_rix(0)
  {
    is_open = this->openStream();
    if (is_open)
      curr_pos = vec_ifs.tellg();
  }
  bool openStream();
  uint64_t getBatch(std::vector<std::vector<encT>>& batch_table, uint32_t batch_size);
};

template<typename encT>
struct tableF
{
  uint8_t k, h, b;
  uint64_t num_rows;
  uint64_t num_kmers;
  size_t table_size;
  maskLSH* ptr_lsh_vg;
  encT* enc_arr;
  uint8_t* ind_arr;

  tableF(uint8_t k, uint8_t h, uint8_t b, uint32_t num_rows, maskLSH* ptr_lsh_vg)
    : k(h)
    , h(h)
    , b(b)
    , num_rows(num_rows) // i.e., batch_size
    , num_kmers(0)
    , ptr_lsh_vg(ptr_lsh_vg)
  {
    enc_arr = new encT[num_rows * b];
    std::fill(enc_arr, enc_arr + num_rows * b, 0);
    ind_arr = new uint8_t[num_rows];
    std::fill(ind_arr, ind_arr + num_rows, 0);
    table_size = num_rows * b * sizeof(encT) + num_rows * sizeof(uint8_t);
  }
  ~tableF()
  {
    delete[] enc_arr;
    delete[] ind_arr;
  }

  void clearRows();
  void makeUnique();
  bool areRowsSorted();
  void sortColumns();
  void updateSize();
  void uniounRows(tableF<encT>& sibling);
  void mergeRows(tableF<encT>& sibling);
  void removeIndices(std::vector<std::pair<uint32_t, uint8_t>>& indices_vec);
  std::unordered_map<uint8_t, uint64_t> histNumCols();
};

template<typename encT>
struct tableD
{
  uint8_t k, h;
  uint64_t num_rows;
  uint64_t num_kmers;
  size_t table_size;
  maskLSH* ptr_lsh_vg;
  std::vector<std::vector<encT>> enc_vvec;

  tableD(uint8_t k, uint8_t h, uint32_t num_rows, maskLSH* ptr_lsh_vg)
    : k(h)
    , h(h)
    , num_rows(num_rows) // i.e., batch_size
    , num_kmers(0)
    , ptr_lsh_vg(ptr_lsh_vg)
  {
    enc_vvec.resize(num_rows);
  }

  void clearRows();
  bool areRowsSorted();
  void sortColumns();
  void makeUnique();
  void updateSize();
  void trimColumns(uint8_t b);
  void pruneColumns(uint8_t b);
  void uniounRows(tableD<encT>& sibling);
  void mergeRows(tableD<encT>& sibling);
  void removeIndices(std::vector<std::pair<uint32_t, uint8_t>>& indices_vec);
  std::unordered_map<uint8_t, uint64_t> histNumCols();
  void transformTableF(tableF<encT>& table);
};

#endif
