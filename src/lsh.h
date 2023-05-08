#ifndef _LSH_H
#define _LSH_H

#include "common.h"

typedef std::vector<std::pair<int8_t, int8_t>> vectorLSH;
typedef std::vector<std::vector<std::pair<int8_t, int8_t>>> vectorMLSH;

vectorLSH
generateMaskLSH(uint8_t k, uint8_t h);

uint64_t
computeValueLSH(uint64_t enc_bp, vectorLSH vg);

#endif
