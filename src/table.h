#ifndef _TABLE_H
#define _TABLE_H

#include "assess.h"
#include "common.h"
#include "encode.h"
#include "io.h"
#include "lsh.h"

typedef uint32_t scT;
template<typename encT>
using vvec = std::vector<std::vector<encT>>;

template<typename encT>
struct StreamIM
{
  uint8_t k, h;
  uint32_t l_rix;
  uint32_t curr_vix;
  maskLSH* ptr_lsh_vg;
  std::vector<std::pair<uint32_t, encT>> lsh_enc_vec;

  StreamIM(uint8_t k, uint8_t h, maskLSH* ptr_lsh_vg)
    : k(h)
    , h(h)
    , l_rix(0)
    , curr_vix(0)
    , ptr_lsh_vg(ptr_lsh_vg)
  {}
  bool save(const char* filepath);
  bool load(const char* filepath);
  uint64_t readBatch(const char* filepath, unsigned int batch_size);
  uint64_t getBatch(vvec<encT>& batch_table, uint32_t batch_size);
  std::unordered_map<uint8_t, uint64_t> histRowSizes();
};

template<typename encT>
struct StreamOD
{
  uint32_t f_rix;
  uint32_t curr_rix;
  const char* filepath;
  std::ifstream vec_ifs;
  std::streampos curr_pos;
  bool is_open;

  StreamOD(const char* filepath)
    : filepath(filepath)
    , curr_rix(0)
    , f_rix(0)
  {
    is_open = this->openStream();
    if (is_open)
      curr_pos = vec_ifs.tellg();
  }
  bool openStream();
  uint64_t getBatch(vvec<encT>& batch_table, uint32_t batch_size);
};

template<typename encT>
struct HTs
{
  uint8_t k, h, b;
  uint32_t num_rows;
  uint64_t num_kmers;
  uint64_t num_species;
  size_t table_size;
  maskLSH* ptr_lsh_vg;
  encT* enc_arr;
  scT* scounts_arr;
  uint8_t* ind_arr;

  HTs(uint8_t k, uint8_t h, uint8_t b, uint32_t num_rows, maskLSH* ptr_lsh_vg)
    : k(h)
    , h(h)
    , b(b)
    , num_rows(num_rows) // i.e., batch_size
    , num_kmers(0)
    , num_species(0)
    , ptr_lsh_vg(ptr_lsh_vg)
  {
    enc_arr = new encT[num_rows * b];
    std::fill(enc_arr, enc_arr + num_rows * b, 0);
    scounts_arr = new scT[num_rows * b];
    std::fill(scounts_arr, scounts_arr + num_rows * b, 0);
    ind_arr = new uint8_t[num_rows];
    std::fill(ind_arr, ind_arr + num_rows, 0);
    table_size = num_rows * b * sizeof(encT) + num_rows * sizeof(uint8_t) + num_rows * b * sizeof(scT);
  }
  ~HTs()
  {
    delete[] enc_arr;
    delete[] ind_arr;
  }

  void clearRows();
  void updateSize();
  void sortColumns();
  bool areColumnsSorted();
  void makeUnique(bool update_size = true);
  void unionRows(HTs<encT>& sibling, bool update_size = true);
  void mergeRows(HTs<encT>& sibling, bool update_size = true);
  void removeIndices(std::vector<std::pair<uint32_t, uint8_t>>& indices_vec); // TODO
  std::unordered_map<uint8_t, uint64_t> histRowSizes();
};

template<typename encT>
struct HTd
{
  uint8_t k, h;
  uint32_t num_rows;
  uint64_t num_kmers;
  uint64_t num_species;
  size_t table_size;
  maskLSH* ptr_lsh_vg;
  vvec<encT> enc_vvec;
  vvec<scT> scounts_vvec;

  HTd(uint8_t k, uint8_t h, uint32_t num_rows, maskLSH* ptr_lsh_vg)
    : k(h)
    , h(h)
    , num_rows(num_rows) // i.e., batch_size
    , num_kmers(0)
    , num_species(0)
    , ptr_lsh_vg(ptr_lsh_vg)
  {
    enc_vvec.resize(num_rows);
    scounts_vvec.resize(num_rows);
    table_size = num_rows * sizeof(std::vector<encT>) + num_rows * sizeof(std::vector<scT>);
  }

  void clearRows();
  void updateSize();
  void initCounts();
  void sortColumns();
  bool areColumnsSorted();
  void makeUnique(bool update_size = true);
  void unionRows(HTd<encT>& sibling, bool update_size = true);
  void mergeRows(HTd<encT>& sibling, bool update_size = true);
  void trimColumns(uint8_t b);
  void pruneColumns(uint8_t b);                                               // TODO
  void removeIndices(std::vector<std::pair<uint32_t, uint8_t>>& indices_vec); // TODO
  std::unordered_map<uint8_t, uint64_t> histRowSizes();
  void transformHTs(HTs<encT>& table);
};

#include "scounts.h"
#endif