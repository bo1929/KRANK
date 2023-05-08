#include "encode.h"

void
kmerEncodingBPCompute(const char* seq, uint64_t& enc_bp)
{
  for (int i = 0; i < (unsigned)strlen(seq); i++) {
    enc_bp = enc_bp << 2;
    if (seq[i] == 'T') {
      enc_bp += 3;
    } else if (seq[i] == 'G') {
      enc_bp += 2;
    } else if (seq[i] == 'C') {
      enc_bp += 1;
    } else {
      enc_bp += 0;
    }
  }
}

void
kmerEncodingCompute(const char* seq, uint64_t& enc_lr, uint64_t& enc_bp)
{
  for (int i = 0; i < (unsigned)strlen(seq); i++) {
    enc_lr = enc_lr << 1;
    enc_bp = enc_bp << 2;
    if (seq[i] == 'T') {
      enc_lr += 4294967297;
      enc_bp += 3;
    } else if (seq[i] == 'G') {
      enc_lr += 4294967296;
      enc_bp += 2;
    } else if (seq[i] == 'C') {
      enc_lr += 1;
      enc_bp += 1;
    } else {
      enc_lr += 0;
      enc_bp += 0;
    }
  }
}

void
kmerEncodingComputeC(const char* seq, uint64_t& enc_lr, uint64_t& enc_bp)
{
  for (int i = 0; i < (unsigned)strlen(seq); i++) {
    enc_lr = enc_lr << 1;
    enc_bp = enc_bp << 2;

    if (seq[i] == 'A') {
      enc_lr += 4294967297;
      enc_bp += 3;
    } else if (seq[i] == 'C') {
      enc_lr += 4294967296;
      enc_bp += 2;
    } else if (seq[i] == 'G') {
      enc_lr += 1;
      enc_bp += 1;
    } else {
      enc_lr += 0;
      enc_bp += 0;
    }
  }
}

void
kmerEncodingUpdate(const char* seq, uint64_t& enc_lr, uint64_t& enc_bp)
{
  enc_lr = enc_lr << 1;
  enc_bp = enc_bp << 2;

  uint64_t mask = 4294967297;
  enc_lr = enc_lr & ~mask;

  if (seq[0] == 'T') {
    enc_lr += 4294967297;
    enc_bp += 3;
  } else if (seq[0] == 'G') {
    enc_lr += 4294967296;
    enc_bp += 2;
  } else if (seq[0] == 'C') {
    enc_lr += 1;
    enc_bp += 1;
  } else {
    enc_lr += 0;
    enc_bp += 0;
  }
}

void
kmerEncodingUpdateC(const char* seq, uint64_t& enc_lr, uint64_t& enc_bp)
{
  enc_lr = enc_lr << 1;
  enc_bp = enc_bp << 2;

  uint64_t mask = 4294967297;
  enc_lr = enc_lr & ~mask;

  if (seq[0] == 'A') {
    enc_lr += 4294967297;
    enc_bp += 3;
  } else if (seq[0] == 'C') {
    enc_lr += 4294967296;
    enc_bp += 2;
  } else if (seq[0] == 'G') {
    enc_lr += 1;
    enc_bp += 1;
  } else {
    enc_lr += 0;
    enc_bp += 0;
  }
}

void
retrieveEncodings(char* fpath, uint64_t*& enc_arr, uint32_t num_kmers, unsigned int batch_size, unsigned int num_threads)
{
  // Allocate memory for the encoding array.
  try {
    enc_arr = new uint64_t[num_kmers];
  } catch (std::bad_alloc& ba) {
    std::cerr << "Failed to allocate memory for the encoding array." << ba.what() << std::endl;
  }
  kseq_t* reader = getReader(fpath);
  batch_size = adjustBatchSize(batch_size, num_threads);
  std::vector<sseq_t> seqBatch = readBatch(reader, batch_size);
  uint32_t total_ix = 0;

  while (!seqBatch.empty()) {
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
    for (uint32_t ix = 0; ix < (uint32_t)seqBatch.size(); ++ix) {
      const char* kmer_seq;
      uint64_t enc_bp;
      uint32_t fix = total_ix + ix;
      assert(fix < UINT32_MAX);

      kmer_seq = seqBatch[ix].nseq.c_str();
      kmerEncodingBPCompute(kmer_seq, enc_bp);
      enc_arr[fix] = enc_bp;
    }
    total_ix += (uint32_t)seqBatch.size();
    if (seqBatch.size() == batch_size)
      seqBatch = readBatch(reader, batch_size);
    else
      seqBatch.clear();
  }
}
