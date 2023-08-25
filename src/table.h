#ifndef _TABLE_H
#define _TABLE_H

#include "assess.h"
#include "common.h"
#include "encode.h"
#include "io.h"
#include "lsh.h"

enum PriorityMethod
{
  random_kmer,
  large_scount,
  nice_scount
};

template<typename encT>
struct StreamIM
{
  uint8_t k, h;
  uint32_t l_rix;
  uint64_t tnum_kmers;
  uint32_t curr_vix;
  const char* filepath;
  maskLSH* ptr_lsh_vg;
  std::vector<std::pair<uint32_t, encT>> lsh_enc_vec;

  StreamIM(const char* filepath, uint8_t k, uint8_t h, maskLSH* ptr_lsh_vg)
    : l_rix(0)
    , curr_vix(0)
    , tnum_kmers(0)
    , filepath(filepath)
    , k(h)
    , h(h)
    , ptr_lsh_vg(ptr_lsh_vg)
  {
  }
  bool save(const char* filepath);
  bool load(const char* filepath);
  void clear();
  uint64_t processInput(uint32_t rbatch_size);
  uint64_t getBatch(vvec<encT>& batch_table, uint32_t tbatch_size);
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
  }
  void openStream();
  uint64_t getBatch(vvec<encT>& batch_table, uint32_t tbatch_size);
};

template<typename encT>
struct HTs
{
  tT tID;
  uint8_t k, h, b;
  uint32_t num_rows;
  uint64_t num_kmers;
  uint64_t num_species;
  size_t table_size;
  maskLSH* ptr_lsh_vg;
  encT* enc_arr;
  scT* scount_arr;
  uint8_t* ind_arr;
  PriorityMethod kmer_priority;
  std::vector<HTs<encT>> childrenHT;
  std::set<tT> tIDsBasis;

  HTs(tT tID, uint8_t k, uint8_t h, uint8_t b, uint32_t num_rows, maskLSH* ptr_lsh_vg, PriorityMethod kmer_priority)
    : tID(tID)
    , k(h)
    , h(h)
    , b(b)
    , num_rows(num_rows) // i.e., tbatch_size
    , num_kmers(0)
    , num_species(0)
    , ptr_lsh_vg(ptr_lsh_vg)
    , kmer_priority(kmer_priority)
  {
    enc_arr = new encT[num_rows * b];
    std::fill(enc_arr, enc_arr + num_rows * b, 0);
    scount_arr = new scT[num_rows * b];
    std::fill(scount_arr, scount_arr + num_rows * b, 0);
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
  void shrinkHT(uint64_t num_rm);
  std::unordered_map<uint8_t, uint64_t> histRowSizes();
};

template<typename encT>
struct HTd
{
  tT tID;
  uint8_t k, h;
  uint32_t num_rows;
  uint64_t num_kmers;
  uint64_t num_species;
  size_t table_size;
  maskLSH* ptr_lsh_vg;
  vvec<encT> enc_vvec;
  vvec<scT> scount_vvec;
  PriorityMethod kmer_priority;
  std::vector<HTd<encT>> childrenHT;
  std::set<tT> tIDsBasis;

  HTd(tT tID, uint8_t k, uint8_t h, uint32_t num_rows, maskLSH* ptr_lsh_vg, PriorityMethod kmer_priority)
    : tID(tID)
    , k(h)
    , h(h)
    , num_rows(num_rows) // i.e., tbatch_size
    , num_kmers(0)
    , num_species(0)
    , ptr_lsh_vg(ptr_lsh_vg)
    , kmer_priority(kmer_priority)
  {
    enc_vvec.resize(num_rows);
    scount_vvec.resize(num_rows);
    table_size = num_rows * sizeof(std::vector<encT>) + num_rows * sizeof(std::vector<scT>);
  }

  void clearRows();
  void updateSize();
  void initBasis(tT tID);
  void sortColumns();
  bool areColumnsSorted();
  void makeUnique(bool update_size = true);
  void unionRows(HTd<encT>& sibling, bool update_size = true);
  void mergeRows(HTd<encT>& sibling, bool update_size = true);
  void trimColumns(uint8_t b_max);
  void pruneColumns(uint8_t b_max);
  void shrinkHT(uint64_t num_rm, uint8_t b_max);
  void convertHTs(HTs<encT>& table);
  std::unordered_map<uint8_t, uint64_t> histRowSizes();
};

#include "scount.h"
#endif
