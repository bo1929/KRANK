#ifndef _TABLE_H
#define _TABLE_H

#include "assess.h"
#include "common.h"
#include "taxonomy.h"
#include "encode.h"
#include "io.h"
#include "lsh.h"

enum RankingMethod
{
  random_kmer,
  representative_kmer,
};

template<typename encT>
struct HTs
{
  tT trID;
  uint8_t k, h, b;
  uint32_t num_rows;
  uint64_t num_kmers;
  uint64_t num_basis;
  size_t table_size;
  maskLSH *ptr_lsh_vg;
  encT *enc_arr;
  scT *scount_arr;
  tT *tlca_arr;
  uint8_t *ind_arr;
  RankingMethod ranking_method;
  TaxonomyRecord<tT> *tax_record;
  std::vector<HTs<encT>> childrenHT;
  std::set<tT> trIDsBasis;

  HTs(tT trID,
      uint8_t k,
      uint8_t h,
      uint8_t b,
      uint32_t num_rows,
      maskLSH *ptr_lsh_vg,
      RankingMethod ranking_method,
      TaxonomyRecord<tT> *tax_record)
    : trID(trID)
    , k(k)
    , h(h)
    , b(b)
    , num_rows(num_rows) // i.e., tbatch_size
    , num_kmers(0)
    , num_basis(0)
    , ptr_lsh_vg(ptr_lsh_vg)
    , ranking_method(ranking_method)
    , tax_record(tax_record)
  {
    enc_arr = new encT[num_rows * b];
    std::fill(enc_arr, enc_arr + num_rows * b, 0);
    scount_arr = new scT[num_rows * b];
    std::fill(scount_arr, scount_arr + num_rows * b, 0);
    tlca_arr = new tT[num_rows * b];
    std::fill(tlca_arr, tlca_arr + num_rows * b, 0);
    ind_arr = new uint8_t[num_rows];
    std::fill(ind_arr, ind_arr + num_rows, 0);
    table_size =
      num_rows * b * sizeof(encT) + num_rows * sizeof(uint8_t) + num_rows * b * sizeof(scT) + num_rows * b * sizeof(tT);
  }
  HTs(const HTs &ts)
    : trID(ts.trID)
    , k(ts.k)
    , h(ts.h)
    , b(ts.b)
    , num_rows(ts.num_rows) // i.e., tbatch_size
    , num_kmers(ts.num_kmers)
    , num_basis(ts.num_basis)
    , ptr_lsh_vg(ts.ptr_lsh_vg)
    , ranking_method(ts.ranking_method)
    , childrenHT(ts.childrenHT)
    , trIDsBasis(ts.trIDsBasis)
  {
    enc_arr = new encT[num_rows * b];
    std::copy(ts.enc_arr, ts.enc_arr + num_rows * b, enc_arr);
    scount_arr = new scT[num_rows * b];
    std::copy(ts.scount_arr, ts.scount_arr + num_rows * b, scount_arr);
    tlca_arr = new tT[num_rows * b];
    std::copy(ts.tlca_arr, ts.tlca_arr + num_rows * b, tlca_arr);
    ind_arr = new uint8_t[num_rows];
    std::copy(ts.ind_arr, ts.ind_arr + num_rows, ind_arr);
    table_size =
      num_rows * b * sizeof(encT) + num_rows * sizeof(uint8_t) + num_rows * b * sizeof(scT) + num_rows * b * sizeof(tT);
  }
  HTs &operator=(const HTs &rhs)
  {
    std::copy(rhs.enc_arr, rhs.enc_arr + num_rows * b, enc_arr);
    std::copy(rhs.scount_arr, rhs.scount_arr + num_rows * b, scount_arr);
    std::copy(rhs.tlca_arr, rhs.tlca_arr + num_rows * b, tlca_arr);
    std::copy(rhs.ind_arr, rhs.ind_arr + num_rows, ind_arr);
    trID = rhs.trID;
    k = rhs.k;
    h = rhs.h;
    b = rhs.b;
    num_rows = rhs.num_rows;
    num_kmers = rhs.num_kmers;
    num_basis = rhs.num_basis;
    ptr_lsh_vg = rhs.ptr_lsh_vg;
    ranking_method = rhs.ranking_method;
    childrenHT = rhs.childrenHT;
    trIDsBasis = rhs.trIDsBasis;
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
  void updateLCA();
  void shrinkHT(uint64_t num_rm);
  void makeUnique(bool update_size = true);
  void unionRows(HTs<encT> &child, bool update_size = true);
  void getScores(std::vector<float> &scores_vec, uint32_t rix);
  void mapBefore(std::unordered_map<encT, scT> &scount_map, std::unordered_map<encT, tT> &tlca_map, uint32_t rix);
  void updateAfter(std::unordered_map<encT, scT> &scount_map, std::unordered_map<encT, tT> &tlca_map, uint32_t rix);
  void accumulateCounts(std::unordered_map<encT, scT> &scount_map, HTs<encT> &child, uint32_t rix);
  void getRowOrdering(std::vector<unsigned int> &row_order, bool reverse);
  std::map<uint8_t, uint64_t> histRowSizes();
};

template<typename encT>
struct HTd
{
  tT trID;
  uint8_t k, h;
  uint32_t num_rows;
  uint64_t num_kmers;
  uint64_t num_basis;
  size_t table_size;
  maskLSH *ptr_lsh_vg;
  vvec<encT> enc_vvec;
  vvec<scT> scount_vvec;
  vvec<tT> tlca_vvec;
  RankingMethod ranking_method;
  TaxonomyRecord<tT> *tax_record;
  std::vector<HTd<encT>> childrenHT;
  std::set<tT> trIDsBasis;

  HTd(tT trID,
      uint8_t k,
      uint8_t h,
      uint32_t num_rows,
      maskLSH *ptr_lsh_vg,
      RankingMethod ranking_method,
      TaxonomyRecord<tT> *tax_record)
    : trID(trID)
    , k(k)
    , h(h)
    , num_rows(num_rows) // i.e., tbatch_size
    , num_kmers(0)
    , num_basis(0)
    , ptr_lsh_vg(ptr_lsh_vg)
    , ranking_method(ranking_method)
    , tax_record(tax_record)
  {
    enc_vvec.resize(num_rows);
    scount_vvec.resize(num_rows);
    tlca_vvec.resize(num_rows);
    table_size =
      num_rows * sizeof(std::vector<encT>) + num_rows * sizeof(std::vector<scT>) + num_rows * sizeof(std::vector<tT>);
  }

  void clearRows();
  void updateSize();
  void sortColumns();
  bool areColumnsSorted();
  void mapBefore(std::unordered_map<encT, scT> &scount_map, std::unordered_map<encT, tT> &tlca_map, uint32_t rix);
  void updateAfter(std::unordered_map<encT, scT> &scount_map, std::unordered_map<encT, tT> &tlca_map, uint32_t rix);
  void accumulateCounts(std::unordered_map<encT, scT> &scount_map, HTd<encT> &child, uint32_t rix);
  void getScores(std::vector<float> &scores_vec, uint32_t rix);
  void filterLSR(std::vector<uint8_t> &depth_vec, uint8_t slr_depth);
  void getRowOrdering(std::vector<unsigned int> &row_order, bool reverse);
  void selectCoverage(std::vector<size_t> &ixs, uint32_t rix, size_t b_max);
  void unionRows(HTd<encT> &child, bool update_size = true);
  void shrinkHT(uint64_t num_rm, size_t b_max);
  void makeUnique(bool update_size = true);
  void convertHTs(HTs<encT> *new_table);
  void trimColumns(size_t b_max);
  void pruneColumns(size_t b_max);
  void initBasis(tT trID);
  void updateLCA();
  std::map<size_t, uint64_t> histRowSizes();
};

#endif
