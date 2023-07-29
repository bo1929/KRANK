#include "assess.h"

inline unsigned int
computeFastModulo(const unsigned int input, const unsigned int ceil)
{
  return input >= ceil ? input % ceil : input;
}

template<typename T>
inline void
pruneListPseudorandom(std::list<T>& l, uint8_t b_max)
{
  while (l.size() > b_max) {
    l.erase(std::next(l.begin(), computeFastModulo(l.back(), l.size())));
  }
}

template<typename T>
inline void
pruneVectorPseudorandom(std::vector<T>& v, uint8_t b_max)
{
  while (v.size() > b_max) {
    v[computeFastModulo(v.back(), v.size())] = v.back();
    v.pop_back();
  }
  std::sort(v.begin(), v.end());
}

template void
pruneListPseudorandom(std::list<uint32_t>& l, uint8_t b_max);

template void
pruneListPseudorandom(std::list<uint64_t>& l, uint8_t b_max);

template void
pruneVectorPseudorandom(std::vector<uint32_t>& v, uint8_t b_max);

template void
pruneVectorPseudorandom(std::vector<uint64_t>& v, uint8_t b_max);
