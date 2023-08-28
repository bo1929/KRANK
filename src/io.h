#ifndef _IO_H
#define _IO_H

#include "common.h"
#include <zlib.h>

extern "C"
{
#include "kseq.h"
}

struct sseq_t
{
  std::string nseq;
  std::string name;
};

KSEQ_INIT(gzFile, gzread)

namespace IO {
bool
ensureDirectory(char* dirpath);

kseq_t*
getReader(const char* path);

uint64_t
adjustBatchSize(uint64_t batch_size, uint8_t num_threads);

std::vector<sseq_t>
readBatch(kseq_t* kseq, uint64_t batch_size);

FILE*
open_file(const char* filepath, bool is_ok, const char* mode);

std::ifstream
open_ifstream(const char* filepath, bool is_ok);
}

#define DEFAULT_BATCH_SIZE 1048576

#endif
