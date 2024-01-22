#include "encode.h"

bool kmerEncodingBPCompute(const char *seq, uint64_t &enc_bp)
{
  enc_bp = 0;
  for (unsigned int i = 0; i < static_cast<unsigned int>(strlen(seq)); i++) {
    enc_bp = enc_bp << 2;
    if (seq[i] == 'T') {
      enc_bp += 3;
    } else if (seq[i] == 'G') {
      enc_bp += 2;
    } else if (seq[i] == 'C') {
      enc_bp += 1;
    } else if (seq[i] == 'A') {
      enc_bp += 0;
    } else {
      return false;
    }
  }
  return true;
}

bool kmerEncodingCompute(const char *seq, uint64_t &enc_lr, uint64_t &enc_bp)
{
  enc_lr = 0;
  enc_bp = 0;
  for (unsigned int i = 0; i < static_cast<unsigned int>(strlen(seq)); i++) {
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
    } else if (seq[i] == 'A') {
      enc_lr += 0;
      enc_bp += 0;
    } else {
      return false;
    }
  }
  return true;
}

bool kmerEncodingComputeC(const char *seq, uint64_t &enc_lr, uint64_t &enc_bp)
{
  enc_lr = 0;
  enc_bp = 0;
  for (unsigned int i = 0; i < static_cast<unsigned int>(strlen(seq)); i++) {
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
    } else if (seq[i] == 'T') {
      enc_lr += 0;
      enc_bp += 0;
    } else {
      return false;
    }
  }
  return true;
}

bool kmerEncodingUpdate(const char *seq, uint64_t &enc_lr, uint64_t &enc_bp)
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
  } else if (seq[0] == 'A') {
    enc_lr += 0;
    enc_bp += 0;
  } else {
    return false;
  }
  return true;
}

bool kmerEncodingUpdateC(const char *seq, uint64_t &enc_lr, uint64_t &enc_bp)
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
  } else if (seq[0] == 'T') {
    enc_lr += 0;
    enc_bp += 0;
  } else {
    return false;
  }
  return true;
}
