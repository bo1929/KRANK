#ifndef _TABLE_H
#define _TABLE_H

#include "common.h"
#include "encode.h"
#include "io.h"
#include "lsh.h"
#include <cstddef>
#include <cstdint>
#include <sys/types.h>

template<typename encT>
struct tableC
{
  uint8_t k, h;
  uint32_t l_rix;
  uint32_t curr_vix;
  uint64_t num_rows;
  maskLSH* ptr_lsh_vg;
  std::vector<std::pair<uint32_t, encT>> rix_enc_vec;

  tableC(uint8_t k, uint8_t h, maskLSH* ptr_lsh_vg)
    : k(h)
    , h(h)
    , l_rix(0)
    , curr_vix(0)
    , ptr_lsh_vg(ptr_lsh_vg)
  {
    num_rows = pow(2, 2 * h);
  }
  std::unordered_map<uint8_t, uint64_t> histNumCols();
  bool saveVec(const char* filepath);
  bool loadVec(const char* filepath);
  uint64_t fillVec(const char* filepath, unsigned int batch_size, unsigned int num_threads);
  uint64_t getBatch(std::vector<std::vector<encT>>& batch_table, uint32_t batch_size);
};

template<typename encT>
struct tableS
{
  uint32_t f_rix;
  uint32_t curr_rix;
  const char* filepath;
  std::ifstream vec_ifs;
  std::streampos curr_pos;

  tableS(const char* filepath)
    : filepath(filepath)
    , curr_rix(0)
    , f_rix(0)
  {
  }
  bool openStream();
  uint64_t getBatch(std::vector<std::vector<encT>>& batch_table, uint32_t batch_size);
};

template<typename encT>
struct tableD
{
  uint8_t k, h;
  uint64_t num_rows;
  maskLSH* ptr_lsh_vg;
  std::vector<std::vector<encT>> dynamic_table;

  tableD(uint8_t k, uint8_t h, uint32_t batch_size, maskLSH* ptr_lsh_vg)
    : k(h)
    , h(h)
    , num_rows(batch_size)
    , ptr_lsh_vg(ptr_lsh_vg)
  {
    dynamic_table.resize(batch_size);
  }

  void makeUnique();
  void pruneColumns(uint8_t b);
};

/*template<typename encT> */
/* struct tableF */
/* { */
/*   uint8_t k, h, b; */
/*   size_t total_size; */
/*   uint64_t num_kmers; */
/*   std::unique_ptr<encT[]> enc_arr; */
/*   std::vector<bool> row_ind; */

/*   tableF(uint8_t k, uint8_t h, uint8_t b) */
/*     : k(k) */
/*     , h(h) */
/*     , b(b) */
/*   { */
/*     total_size = pow(2, 2 * h) * b; */
/*   } */
/* }; */

#endif
