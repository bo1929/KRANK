#ifndef _TABLE_H
#define _TABLE_H

#include "assess.h"
#include "common.h"
#include "encode.h"
#include "io.h"
#include "lsh.h"

enum RankingMethod
{
  random_kmer,
  species_count,
  information_score,
  taxa_count
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
  std::vector<uint8_t>* ptr_npositions;
  std::vector<std::pair<uint32_t, encT>> lsh_enc_vec;

  StreamIM(const char* filepath,
           uint8_t k,
           uint8_t h,
           maskLSH* ptr_lsh_vg,
           std::vector<uint8_t>* ptr_npositions)
    : l_rix(0)
    , curr_vix(0)
    , tnum_kmers(0)
    , filepath(filepath)
    , k(k)
    , h(h)
    , ptr_lsh_vg(ptr_lsh_vg)
    , ptr_npositions(ptr_npositions)
  {
  }
  bool save(const char* filepath);
  bool load(const char* filepath);
  void clearStream();
  void resetStream();
  uint64_t processInput(uint64_t rbatch_size);
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
  void closeStream();
  uint64_t getBatch(vvec<encT>& batch_table, uint32_t tbatch_size, bool contd = false);
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
  tT* tlca_arr;
  uint8_t* ind_arr;
  RankingMethod ranking_method;
  std::vector<HTs<encT>> childrenHT;
  std::set<tT> tIDsBasis;

  HTs(tT tID,
      uint8_t k,
      uint8_t h,
      uint8_t b,
      uint32_t num_rows,
      maskLSH* ptr_lsh_vg,
      RankingMethod ranking_method)
    : tID(tID)
    , k(k)
    , h(h)
    , b(b)
    , num_rows(num_rows) // i.e., tbatch_size
    , num_kmers(0)
    , num_species(0)
    , ptr_lsh_vg(ptr_lsh_vg)
    , ranking_method(ranking_method)
  {
    enc_arr = new encT[num_rows * b];
    std::fill(enc_arr, enc_arr + num_rows * b, 0);
    scount_arr = new scT[num_rows * b];
    std::fill(scount_arr, scount_arr + num_rows * b, 0);
    tlca_arr = new tT[num_rows * b];
    std::fill(tlca_arr, tlca_arr + num_rows * b, 0);
    ind_arr = new uint8_t[num_rows];
    std::fill(ind_arr, ind_arr + num_rows, 0);
    table_size = num_rows * b * sizeof(encT) + num_rows * sizeof(uint8_t) +
                 num_rows * b * sizeof(scT) + num_rows * b * sizeof(tT);
  }
  HTs(const HTs& ts)
    : tID(ts.tID)
    , k(ts.k)
    , h(ts.h)
    , b(ts.b)
    , num_rows(ts.num_rows) // i.e., tbatch_size
    , num_kmers(ts.num_kmers)
    , num_species(ts.num_species)
    , ptr_lsh_vg(ts.ptr_lsh_vg)
    , ranking_method(ts.ranking_method)
    , childrenHT(ts.childrenHT)
    , tIDsBasis(ts.tIDsBasis)
  {
    enc_arr = new encT[num_rows * b];
    std::copy(ts.enc_arr, ts.enc_arr + num_rows * b, enc_arr);
    scount_arr = new scT[num_rows * b];
    std::copy(ts.scount_arr, ts.scount_arr + num_rows * b, scount_arr);
    tlca_arr = new tT[num_rows * b];
    std::copy(ts.tlca_arr, ts.tlca_arr + num_rows * b, tlca_arr);
    ind_arr = new uint8_t[num_rows];
    std::copy(ts.ind_arr, ts.ind_arr + num_rows, ind_arr);
    table_size = num_rows * b * sizeof(encT) + num_rows * sizeof(uint8_t) +
                 num_rows * b * sizeof(scT) + num_rows * b * sizeof(tT);
  }
  HTs& operator=(const HTs& rhs)
  {
    std::copy(rhs.enc_arr, rhs.enc_arr + num_rows * b, enc_arr);
    std::copy(rhs.scount_arr, rhs.scount_arr + num_rows * b, scount_arr);
    std::copy(rhs.tlca_arr, rhs.tlca_arr + num_rows * b, tlca_arr);
    std::copy(rhs.ind_arr, rhs.ind_arr + num_rows, ind_arr);
    tID = rhs.tID;
    k = rhs.k;
    h = rhs.h;
    b = rhs.b;
    num_rows = rhs.num_rows;
    num_kmers = rhs.num_kmers;
    num_species = rhs.num_species;
    ptr_lsh_vg = rhs.ptr_lsh_vg;
    ranking_method = rhs.ranking_method;
    childrenHT = rhs.childrenHT;
    tIDsBasis = rhs.tIDsBasis;
    return *this;
  }

  ~HTs(void)
  {
    delete[] enc_arr;
    delete[] scount_arr;
    delete[] tlca_arr;
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
  void updateLCA();
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
  vvec<tT> tlca_vvec;
  RankingMethod ranking_method;
  std::vector<HTd<encT>> childrenHT;
  std::set<tT> tIDsBasis;

  HTd(tT tID,
      uint8_t k,
      uint8_t h,
      uint32_t num_rows,
      maskLSH* ptr_lsh_vg,
      RankingMethod ranking_method)
    : tID(tID)
    , k(k)
    , h(h)
    , num_rows(num_rows) // i.e., tbatch_size
    , num_kmers(0)
    , num_species(0)
    , ptr_lsh_vg(ptr_lsh_vg)
    , ranking_method(ranking_method)
  {
    enc_vvec.resize(num_rows);
    scount_vvec.resize(num_rows);
    tlca_vvec.resize(num_rows);
    table_size = num_rows * sizeof(std::vector<encT>) + num_rows * sizeof(std::vector<scT>) +
                 num_rows * sizeof(std::vector<tT>);
  }

  void clearRows();
  void updateSize();
  void sortColumns();
  bool areColumnsSorted();
  void makeUnique(bool update_size = true);
  void unionRows(HTd<encT>& sibling, bool update_size = true);
  void mergeRows(HTd<encT>& sibling, bool update_size = true);
  void convertHTs(HTs<encT>* new_table);
  void trimColumns(uint8_t b_max);
  void pruneColumns(uint8_t b_max);
  void shrinkHT(uint64_t num_rm, uint8_t b_max);
  void initBasis(tT tID);
  void updateLCA();
  std::unordered_map<uint8_t, uint64_t> histRowSizes();
};

#include "lca.h"
#include "scount.h"
#endif
