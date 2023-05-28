#include "common.h"
#include <cstdint>

uint8_t
computeHammingDistance(uint64_t x, uint64_t y)
{
  uint64_t z1 = x ^ y;
  uint32_t z2 = z1 >> 32;
  uint32_t zc = z1 | z2;
  return (uint8_t)__builtin_popcount(zc);
}
