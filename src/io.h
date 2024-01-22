#ifndef _IO_H
#define _IO_H

#include "common.h"
#include "encode.h"
#include "lsh.h"

extern "C"
{
#include "kseq.h"
}

struct sseq_t
{
  std::string nseq;
  std::string name;
  unsigned int len;
};

KSEQ_INIT(gzFile, gzread)

namespace IO {
  bool ensureDirectory(const char *dirpath);

  kseq_t *getReader(const char *path);

  uint64_t adjustBatchSize(uint64_t batch_size, uint8_t num_threads);

  void readBatch(std::vector<sseq_t> &seqRead, kseq_t *kseq, uint64_t batch_size);

  FILE *open_file(const char *filepath, bool &is_ok, const char *mode);

  std::ifstream open_ifstream(const char *filepath, bool is_ok);
} // namespace IO

template<typename encT>
struct StreamIM
{
  uint8_t k, w, h;
  uint32_t l_rix;
  uint64_t tnum_kmers;
  uint32_t curr_vix;
  std::vector<std::string> filepath_v;
  maskLSH *ptr_lsh_vg;
  std::vector<uint8_t> *ptr_npositions;
  std::vector<std::pair<uint32_t, encT>> lsh_enc_vec;

  StreamIM(std::vector<std::string> filepath_v,
           uint8_t k,
           uint8_t w,
           uint8_t h,
           maskLSH *ptr_lsh_vg,
           std::vector<uint8_t> *ptr_npositions)
    : l_rix(0)
    , curr_vix(0)
    , tnum_kmers(0)
    , filepath_v(filepath_v)
    , k(k)
    , h(h)
    , w(w)
    , ptr_lsh_vg(ptr_lsh_vg)
    , ptr_npositions(ptr_npositions)
  {}
  bool save(const char *filepath);
  bool load(const char *filepath);
  void clearStream();
  void resetStream();
  uint64_t readInput(uint64_t rbatch_size);
  uint64_t extractInput(uint64_t rbatch_size);
  uint64_t getBatch(vvec<encT> &batch_table, uint32_t tbatch_size);
  std::map<uint8_t, uint64_t> histRowSizes();
};

template<typename encT>
struct StreamOD
{
  uint32_t f_rix;
  uint32_t curr_rix;
  std::string filepath;
  std::ifstream vec_ifs;
  std::streampos curr_pos;
  bool is_open;

  StreamOD(std::string filepath)
    : filepath(filepath)
    , curr_rix(0)
    , f_rix(0)
    , is_open(false)
  {}
  void openStream();
  void closeStream();
  uint64_t getBatch(vvec<encT> &batch_table, uint32_t tbatch_size, bool contd = false);
  void load(std::vector<std::pair<uint32_t, encT>> &lsh_enc_vec, uint32_t bix, uint32_t eix);
};

#define DEFAULT_BATCH_SIZE 1048576
#define GENOME_BATCH_SIZE 50

#endif
