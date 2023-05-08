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

kseq_t*
getReader(char* path);

unsigned int
adjustBatchSize(unsigned int batch_size, unsigned int num_threads);

std::vector<sseq_t>
readBatch(kseq_t* kseq, unsigned int batch_size);

#define DEFAULT_BATCH_SIZE 65536

#endif
