#ifndef _ENCODE_H
#define _ENCODE_H

#include "common.h"
#include "io.h"

void
kmerEncodingBPCompute(const char* seq, uint64_t& enc_bp);

void
kmerEncodingCompute(const char* seq, uint64_t& enc_lr, uint64_t& enc_bp);

void
kmerEncodingComputeC(const char* seq, uint64_t& enc_lr, uint64_t& enc_bp);

void
kmerEncodingUpdate(const char* seq, uint64_t& enc_lr, uint64_t& enc_bp);

void
kmerEncodingUpdateC(const char* seq, uint64_t& enc_lr, uint64_t& enc_bp);

void
retrieveEncodings(const char* fpath,
                  uint64_t*& enc_arr,
                  uint32_t num_kmers,
                  unsigned int batch_size);

uint8_t
computeHammingDistance64(uint64_t x, uint64_t y);

uint8_t
computeHammingDistance32(uint32_t x, uint32_t y);

inline void
drop64Encoding32(std::vector<uint8_t>& npositions,
                 uint64_t enc64_bp,
                 uint64_t enc64_lr,
                 uint32_t& enc32_bp,
                 uint32_t& enc32_lr)
{
  assert(npositions.size() <= 16);
  enc32_bp = 0;
  enc32_lr = 0;
  for (unsigned int i = npositions.size() - 1; i >= 1; --i) {
    enc32_lr += static_cast<uint32_t>((enc64_lr >> npositions[i]) & 1);
    enc32_lr += static_cast<uint32_t>((enc64_lr >> (npositions[i] + 32)) & 1) << 16;
    enc32_lr = enc32_lr << 1;
    enc32_bp += static_cast<uint32_t>((enc64_bp >> (npositions[i] * 2)) & 3);
    enc32_bp = enc32_bp << 2;
  }
  enc32_lr += static_cast<uint32_t>((enc64_lr >> npositions.front()) & 1);
  enc32_lr += static_cast<uint32_t>((enc64_lr >> (npositions.front() + 32)) & 1) << 16;
  enc32_bp += static_cast<uint32_t>((enc64_bp >> (npositions.front() * 2)) & 3);
}

inline uint8_t
computeHammingDistance64(const uint64_t x, const uint64_t y)
{
  uint64_t z1 = x ^ y;
  uint32_t z2 = z1 >> 32;
  uint32_t zc = z1 | z2;
  return static_cast<uint8_t>(__builtin_popcount(zc));
}

inline uint8_t
computeHammingDistance32(const uint32_t x, const uint32_t y)
{
  uint32_t z1 = x ^ y;
  uint16_t z2 = z1 >> 16;
  uint16_t zc = z1 | z2;
  return static_cast<uint8_t>(__builtin_popcount(zc));
}

inline u_int64_t
revcomp64bp(const u_int64_t& x, size_t sizeKmer)
{
  uint64_t res = ~x;
  res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
  res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
  res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
  res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
  res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
  return (res >> (2 * (32 - sizeKmer)));
}

inline uint64_t
rmoddp64(uint64_t x)
{
  x = x & 0x5555555555555555;
  x = (x | (x >> 1)) & 0x3333333333333333;
  x = (x | (x >> 2)) & 0x0f0f0f0f0f0f0f0f;
  x = (x | (x >> 4)) & 0x00ff00ff00ff00ff;
  x = (x | (x >> 8)) & 0x0000ffff0000ffff;
  x = (x | (x >> 16)) & 0x00000000ffffffff;
  return x;
}

inline u_int64_t
conv64bp2lr(const u_int64_t& x, size_t sizeKmer)
{
  return (rmoddp64(x >> 1) << 32) | rmoddp64(x);
}

static inline uint64_t
murmur64(uint64_t h)
{
  h ^= (h >> 33);
  h *= 0xff51afd7ed558ccdL;
  h ^= (h >> 33);
  h *= 0xc4ceb9fe1a85ec53L;
  h ^= (h >> 33);
  return h;
}

#endif
