#include "lsh.h"

maskLSH
generateMaskLSH(std::vector<uint8_t>& positions)
{
  sort(positions.begin(), positions.end(), std::greater<uint8_t>());
  std::vector<int8_t> v;
  std::vector<int8_t> g;
  int8_t lp = 31;
  int8_t jp = 0;

  for (int8_t j = 0; j < positions.size(); j++) {
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
  return static_cast<uint32_t>(res);
}
