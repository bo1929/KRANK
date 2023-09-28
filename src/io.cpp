#include "io.h"
#include <cstdint>
#include <cstring>
#include <fstream>
#include <utility>

bool
IO::ensureDirectory(const char* dirpath)
{
  bool is_ok;
  struct stat info;
  if (stat(dirpath, &info) != 0) {
    std::cerr << "Cannot access to given directory path " << dirpath << std::endl;
    is_ok = false;
  } else if (info.st_mode & S_IFDIR) {
    is_ok = true;
  } else {
    std::cerr << "Given path is not a directory " << dirpath << std::endl;
    is_ok = false;
  }
  return is_ok;
}

kseq_t*
IO::getReader(const char* fpath)
{
  gzFile fp;
  fp = gzopen(fpath, "r");
  if (fp == nullptr) {
    std::cerr << "Failed to open file at" << fpath << std::endl;
    exit(1);
  }

  kseq_t* kseq;
  kseq = kseq_init(fp);

  return kseq;
}

uint64_t
IO::adjustBatchSize(uint64_t batch_size, uint8_t num_threads)
{
  return (batch_size / num_threads) * num_threads;
}

std::vector<sseq_t>
IO::readBatch(kseq_t* kseq, uint64_t batch_size)
{
  int l;
  unsigned int i = 0;
  std::vector<sseq_t> seqRead;
  seqRead.reserve(batch_size);
  while ((i < batch_size) && (l = kseq_read(kseq)) >= 0) {
    sseq_t tmp_seq;
    tmp_seq.nseq.assign(kseq->seq.s);
    tmp_seq.name.assign(kseq->name.s);
    seqRead.push_back(tmp_seq);
    i++;
  }
  if ((i < batch_size) && (!seqRead.empty())) {
    kseq_destroy(kseq);
    gzclose(kseq->f->f);
  }
  return seqRead;
}

FILE*
IO::open_file(const char* filepath, bool& is_ok, const char* mode)
{
  FILE* f;
  f = std::fopen(filepath, mode);
  if (!f) {
    std::cerr << "File opening failed!" << filepath << std::endl;
    is_ok = false;
  }
  return f;
}

std::ifstream
IO::open_ifstream(const char* filepath, bool is_ok)
{
  std::ifstream ifs;
  ifs.open(filepath, std::ios::binary | std::ios::in);
  if (!ifs.is_open()) {
    std::fprintf(stderr, "File opening failed! %s", filepath);
    is_ok = false;
  }
  return ifs;
}

void
IO::read_encinput(std::string disk_path, std::vector<std::pair<uint32_t, encT>>& lsh_enc_vec)
{
  bool is_ok = false;
  std::ifstream vec_ifs = open_ifstream(disk_path.c_str(), is_ok);
  while (!vec_ifs.eof() && vec_ifs.good()) {
    std::pair<uint32_t, encT> lsh_enc;
    vec_ifs.read((char*)&lsh_enc, sizeof(std::pair<uint32_t, encT>));
    if (!vec_ifs.eof())
      lsh_enc_vec.push_back(lsh_enc);
  }
  if (vec_ifs.fail() && !vec_ifs.eof()) {
    std::puts("I/O rror when reading LSH-value and encoding pairs.\n");
  } else if (vec_ifs.eof()) {
    is_ok = true;
  }
  vec_ifs.close();
}
