#ifndef _LSH_H
#define _LSH_H

#include "common.h"

typedef std::vector<std::pair<int8_t, int8_t>> maskLSH;

maskLSH
generateMaskLSH(uint8_t k, uint8_t h);

uint32_t
computeValueLSH(uint64_t enc_bp, maskLSH vg);

#endif
