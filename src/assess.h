#ifndef _ASSESS_H
#define _ASSESS_H
#include "common.h"

inline unsigned int
computeFastModulo(const unsigned int input, const unsigned int ceil);

template<typename T>
void
pruneVectorPseudorandom(std::vector<T>& v, uint8_t b_max);

template<typename T>
void
pruneListPseudorandom(std::list<T>& l, uint8_t b_max);

#endif
