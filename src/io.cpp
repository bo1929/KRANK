#include "io.h"
#include <cstring>

kseq_t*
getReader(char* fpath)
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

unsigned int
adjustBatchSize(unsigned int batch_size, unsigned int num_threads)
{
  return (batch_size / num_threads) * num_threads;
}

std::vector<sseq_t>
readBatch(kseq_t* kseq, unsigned int batch_size)
{
  int l;
  int i = 0;
  std::vector<sseq_t> seqRead;
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
