#include "lsh.h"
#include <cstdint>

maskLSH
generateMaskLSH(uint8_t k, uint8_t h)
{
  uint8_t n;
  std::vector<uint8_t> positions;
  for (uint8_t m = 0; m < h; m++) {
    n = rand() % k;
    if (count(positions.begin(), positions.end(), n)) {
      m -= 1;
    } else {
      positions.push_back(n);
    }
  }
  sort(positions.begin(), positions.end(), std::greater<uint8_t>());

  // VERBOSE START
#ifdef DEBUG
  std::cout << "LSH positions: ";
  for (unsigned int j = 0; j < (h - 1); j++) {
    std::cout << (unsigned int)positions[j] << ", ";
  }
  std::cout << (unsigned int)positions[h - 1] << std::endl;
  char* test_kmer = new char[k];
  for (uint8_t i = 0; i < k; i++) {
    uint8_t random_char = rand() % 4;
    if (random_char == 3) {
      test_kmer[i] = 'T';
    } else if (random_char == 2) {
      test_kmer[i] = 'G';
    } else if (random_char == 1) {
      test_kmer[i] = 'C';
    } else {
      test_kmer[i] = 'A';
    }
  }
  std::cout << test_kmer << std::endl;
  std::cout << "Grabbed bits: ";
  for (unsigned int j = 0; j < h; j++) {
    std::cout << test_kmer[31 - positions[j]];
  }
  std::cout << std::endl;
  std::cout << "Binary LSH: ";
  for (unsigned int j = 0; j < h; j++) {
    if (test_kmer[31 - positions[j]] == 'T') {
      std::cout << "11";
    } else if (test_kmer[31 - positions[j]] == 'G') {
      std::cout << "10";
    } else if (test_kmer[31 - positions[j]] == 'C') {
      std::cout << "01";
    } else {
      std::cout << "00";
    }
  }
  std::cout << std::endl;
#endif
  // VERBOSE END

  std::vector<int8_t> v;
  std::vector<int8_t> g;
  int8_t lp = 31;
  int8_t jp = 0;

  for (int8_t j = 0; j < h; j++) {
    if (j == 0) {
      v.push_back((lp - positions[j]) * 2);
      lp = positions[j];
      jp += 2;
    } else if ((positions[j - 1] - positions[j]) != 1) {
      v.push_back((lp - positions[j]) * 2);
      lp = positions[j];
      g.push_back(jp);
      jp = 2;
    } else {
      jp += 2;
    }
  }
  g.push_back(jp);

  // Complete above one-dimensional vector.
  v.push_back(-1);
  g.push_back(-1);

  maskLSH lsh_vg(v.size() < g.size() ? v.size() : g.size());
  for (unsigned int i = 0; i < lsh_vg.size(); i++) {
    lsh_vg[i] = std::make_pair(v[i], g[i]);
  }
  return lsh_vg;
}

uint32_t
computeValueLSH(uint64_t enc_bp, maskLSH lsh_vg)
{
  uint64_t res = 0;
  unsigned int i = 0;
  while (lsh_vg[i].first != -1) {
    enc_bp = enc_bp << lsh_vg[i].first;
    asm("shld %b3, %2, %0" : "=rm"(res) : "0"(res), "r"(enc_bp), "ic"(lsh_vg[i].second) : "cc");
    i++;
  }
  return (uint32_t)res;
}
