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
  while (v.size() - b_max) {
    v[computeFastModulo(v.back(), v.size())] = v.back();
    v.pop_back();
  }
  std::sort(v.begin(), v.end());
  // Somewhat slower alternative
  /* std::list<T> l; */
  /* std::copy(v.begin(), v.end(), std::back_inserter(l)); */
  /* pruneListPseudorandom(l, b_max); */
  /* v.resize(b_max); */
  /* std::copy(std::begin(l), std::end(l), std::begin(v)); */
}

template void
pruneListPseudorandom(std::list<uint32_t>& l, uint8_t b_max);

template void
pruneListPseudorandom(std::list<uint64_t>& l, uint8_t b_max);

template void
pruneVectorPseudorandom(std::vector<uint32_t>& v, uint8_t b_max);

template void
pruneVectorPseudorandom(std::vector<uint64_t>& v, uint8_t b_max);
