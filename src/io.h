#ifndef _IO_H
#define _IO_H

#include "lsh.h"
#include "common.h"
#include "encode.h"
#include "filesystem.hpp"

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

  size_t writeData(void *ptr, size_t size, size_t nmemb, FILE *stream);

  std::string downloadURL(std::string url);

  bool checkFASTAQ(const char *filepath);
}

template<typename encT>
struct inputHandler
{
  uint8_t k, w, h;
  uint32_t l_rix;
  uint64_t tnum_kmers;
  uint32_t curr_vix;
  std::vector<std::string> filepath_v;
  maskLSH *ptr_lsh_vg;
  std::vector<uint8_t> *ptr_npositions;
  std::vector<std::pair<uint32_t, encT>> lsh_enc_vec;
  std::unordered_map<encT, uint64_t> rcounts;

  inputHandler(std::vector<std::string> filepath_v,
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
  float readInput(uint64_t rbatch_size);
  void extractMers(std::vector<std::pair<uint32_t, encT>> &lsh_enc_vec_f, std::string filepath, uint64_t rbatch_size);
  void readMers(std::vector<std::pair<uint32_t, encT>> &lsh_enc_vec_f, std::string filepath, uint64_t rbatch_size);
  float extractInput(uint64_t rbatch_size);
  bool saveInput(const char *dirpath, tT trID_key, uint16_t total_batches, uint32_t tbatch_size);
  bool loadInput(const char *dirpath, tT trID_key, uint16_t total_batches);
  bool checkInput(const char *dirpath, tT trID_key, uint16_t total_batches);
  uint64_t getBatch(vvec<encT> &td, uint32_t tbatch_size);
  std::map<uint8_t, uint64_t> histRowSizes();
  void clearInput();
  void resetInput();
};

template<typename encT>
struct inputStream
{
  tT trID_key;
  std::string dirpath;

  inputStream(std::string dirpath, tT trID_key)
    : dirpath(dirpath)
    , trID_key(trID_key)
  {}
  void loadBatch(std::vector<std::pair<uint32_t, encT>> &lsh_enc_vec, unsigned int curr_batch);
  void loadCounts(std::unordered_map<encT, uint32_t> &rcounts);
  uint64_t retrieveBatch(vvec<encT> &td, uint32_t tbatch_size, unsigned int curr_batch, bool shared_table = false);
  void removeBatch(unsigned int curr_batch);
};

inline bool exists_test(const char *filepath)
{
  if (FILE *file = fopen(filepath, "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }
}

#define DEFAULT_BATCH_SIZE 1048576
#define GENOME_BATCH_SIZE 16

#endif
