#ifndef _INDEX_H
#define _INDEX_H

#include "common.h"
#include "encode.h"
#include "io.h"
#include "lsh.h"

uint32_t
fillIndexTable(char* fpath, std::vector<std::vector<uint32_t>>& index_table, uint8_t k, uint8_t h, unsigned int batch_size, unsigned int num_threads);

std::vector<uint64_t>
computeHistNumCols(std::vector<std::vector<uint32_t>> tableIx);

#endif
